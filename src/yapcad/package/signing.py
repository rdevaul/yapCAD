"""Package signing and verification for yapCAD packages.

This module provides cryptographic signing and verification of .ycpkg packages
using GPG or SSH keys. It enables individual contributors to sign packages
without requiring a centralized authority.

Usage:
    from yapcad.package.signing import sign_package, verify_package
    
    # Sign with GPG
    sign_package("my_design.ycpkg", method="gpg", key_id="ABCD1234")
    
    # Sign with SSH
    sign_package("my_design.ycpkg", method="ssh", key_path="~/.ssh/id_ed25519")
    
    # Verify
    result = verify_package("my_design.ycpkg")
    if result.is_valid:
        print(f"Verified: {result.trusted_signers} trusted signatures")
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import tempfile
from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml


class SignatureMethod(Enum):
    """Supported signature methods."""
    GPG = "gpg"
    SSH = "ssh"


class VerificationStatus(Enum):
    """Verification result status."""
    VALID = "VALID"              # Signature valid and signer trusted
    VALID_UNTRUSTED = "VALID_UNTRUSTED"  # Signature valid but signer not trusted
    INVALID = "INVALID"          # Signature doesn't match
    ERROR = "ERROR"              # Verification failed (missing key, etc.)


@dataclass
class SignatureInfo:
    """Information about a single signature."""
    signer: str
    method: SignatureMethod
    key_id: str
    timestamp: str
    signature_file: str
    status: Optional[VerificationStatus] = None
    error_message: Optional[str] = None


@dataclass
class VerificationResult:
    """Result of package verification."""
    package_path: str
    manifest_hash: str
    signatures: List[SignatureInfo] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    
    @property
    def is_valid(self) -> bool:
        """True if at least one trusted signature is valid."""
        return any(s.status == VerificationStatus.VALID for s in self.signatures)
    
    @property
    def trusted_signers(self) -> int:
        """Number of valid, trusted signatures."""
        return sum(1 for s in self.signatures if s.status == VerificationStatus.VALID)
    
    @property
    def untrusted_signers(self) -> int:
        """Number of valid but untrusted signatures."""
        return sum(1 for s in self.signatures if s.status == VerificationStatus.VALID_UNTRUSTED)
    
    @property
    def overall_status(self) -> VerificationStatus:
        """Overall verification status."""
        if self.errors:
            return VerificationStatus.ERROR
        if self.trusted_signers > 0:
            return VerificationStatus.VALID
        if self.untrusted_signers > 0:
            return VerificationStatus.VALID_UNTRUSTED
        if any(s.status == VerificationStatus.INVALID for s in self.signatures):
            return VerificationStatus.INVALID
        return VerificationStatus.ERROR


class SigningError(Exception):
    """Error during signing operation."""
    pass


class VerificationError(Exception):
    """Error during verification operation."""
    pass


def _compute_canonical_manifest_hash(manifest_path: Path) -> str:
    """Compute canonical hash of manifest for signing.
    
    The manifest is canonicalized by:
    1. Loading as YAML
    2. Removing the 'signatures' section
    3. Sorting all keys recursively
    4. Serializing to JSON
    5. Computing SHA-256
    """
    with open(manifest_path, 'r', encoding='utf-8') as f:
        data = yaml.safe_load(f)
    
    # Remove signatures section (we're signing the content, not the signatures)
    if 'signatures' in data:
        data = {k: v for k, v in data.items() if k != 'signatures'}
    
    # Canonicalize: sort keys, compact JSON
    canonical = json.dumps(data, sort_keys=True, separators=(',', ':'))
    
    # Compute hash
    hash_bytes = hashlib.sha256(canonical.encode('utf-8')).digest()
    return hash_bytes.hex()


def _create_signable_content(manifest_hash: str) -> bytes:
    """Create the content to be signed.
    
    Format:
        ycpkg-manifest-v1\n
        <sha256_hex>\n
    """
    content = f"ycpkg-manifest-v1\n{manifest_hash}\n"
    return content.encode('utf-8')


def _sign_with_gpg(content: bytes, key_id: Optional[str] = None) -> tuple[str, str, str]:
    """Sign content with GPG.
    
    Returns: (signature_armored, signer_id, key_id)
    """
    cmd = ["gpg", "--armor", "--detach-sign"]
    if key_id:
        cmd.extend(["--local-user", key_id])
    
    try:
        result = subprocess.run(
            cmd,
            input=content,
            capture_output=True,
            check=True
        )
        signature = result.stdout.decode('utf-8')
    except subprocess.CalledProcessError as e:
        raise SigningError(f"GPG signing failed: {e.stderr.decode('utf-8')}")
    except FileNotFoundError:
        raise SigningError("GPG not found. Install GnuPG to use GPG signing.")
    
    # Get signer info
    try:
        list_result = subprocess.run(
            ["gpg", "--list-keys", "--with-colons", key_id or ""],
            capture_output=True,
            check=True
        )
        lines = list_result.stdout.decode('utf-8').split('\n')
        signer_id = ""
        actual_key_id = key_id or ""
        for line in lines:
            parts = line.split(':')
            if parts[0] == 'uid' and not signer_id:
                signer_id = parts[9]
            elif parts[0] == 'pub':
                actual_key_id = parts[4]
    except Exception:
        signer_id = key_id or "unknown"
        actual_key_id = key_id or "unknown"
    
    return signature, signer_id, actual_key_id


def _sign_with_ssh(content: bytes, key_path: Path, namespace: str = "ycpkg") -> tuple[str, str, str]:
    """Sign content with SSH key.
    
    Returns: (signature, signer_id, key_fingerprint)
    """
    key_path = Path(key_path).expanduser()
    if not key_path.exists():
        raise SigningError(f"SSH key not found: {key_path}")
    
    # Write content to temp file (ssh-keygen needs file input)
    with tempfile.NamedTemporaryFile(mode='wb', delete=False) as f:
        f.write(content)
        content_file = f.name
    
    try:
        # Sign
        sig_file = content_file + ".sig"
        result = subprocess.run(
            ["ssh-keygen", "-Y", "sign", "-f", str(key_path), "-n", namespace, content_file],
            capture_output=True,
            check=True
        )
        
        # Read signature
        with open(sig_file, 'r') as f:
            signature = f.read()
        
        # Get key fingerprint
        fp_result = subprocess.run(
            ["ssh-keygen", "-lf", str(key_path)],
            capture_output=True,
            check=True
        )
        fp_line = fp_result.stdout.decode('utf-8').strip()
        # Format: "256 SHA256:xxxx comment (type)"
        parts = fp_line.split()
        fingerprint = parts[1] if len(parts) > 1 else "unknown"
        signer_id = parts[2] if len(parts) > 2 else str(key_path.name)
        
    except subprocess.CalledProcessError as e:
        raise SigningError(f"SSH signing failed: {e.stderr.decode('utf-8')}")
    except FileNotFoundError:
        raise SigningError("ssh-keygen not found. OpenSSH 8.0+ required for SSH signing.")
    finally:
        # Cleanup
        Path(content_file).unlink(missing_ok=True)
        Path(content_file + ".sig").unlink(missing_ok=True)
    
    return signature, signer_id, fingerprint


def _verify_gpg_signature(content: bytes, signature: str, keyring: Optional[Path] = None) -> tuple[VerificationStatus, str, Optional[str]]:
    """Verify GPG signature.
    
    Returns: (status, signer_id, error_message)
    """
    with tempfile.NamedTemporaryFile(mode='wb', delete=False, suffix='.content') as f:
        f.write(content)
        content_file = f.name
    
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.sig') as f:
        f.write(signature)
        sig_file = f.name
    
    try:
        cmd = ["gpg", "--verify", sig_file, content_file]
        if keyring:
            cmd.extend(["--keyring", str(keyring)])
        
        result = subprocess.run(cmd, capture_output=True)
        stderr = result.stderr.decode('utf-8')
        
        if result.returncode == 0:
            # Extract signer from output
            signer = "unknown"
            for line in stderr.split('\n'):
                if 'Good signature from' in line:
                    # Extract quoted name
                    start = line.find('"') + 1
                    end = line.rfind('"')
                    if start > 0 and end > start:
                        signer = line[start:end]
            return VerificationStatus.VALID, signer, None
        elif "Can't check signature: No public key" in stderr:
            return VerificationStatus.VALID_UNTRUSTED, "unknown", "Public key not found"
        else:
            return VerificationStatus.INVALID, "unknown", stderr
            
    except FileNotFoundError:
        return VerificationStatus.ERROR, "unknown", "GPG not found"
    finally:
        Path(content_file).unlink(missing_ok=True)
        Path(sig_file).unlink(missing_ok=True)


def _verify_ssh_signature(
    content: bytes, 
    signature: str, 
    allowed_signers: Optional[Path] = None,
    namespace: str = "ycpkg"
) -> tuple[VerificationStatus, str, Optional[str]]:
    """Verify SSH signature.
    
    Returns: (status, signer_id, error_message)
    """
    with tempfile.NamedTemporaryFile(mode='wb', delete=False, suffix='.content') as f:
        f.write(content)
        content_file = f.name
    
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.sig') as f:
        f.write(signature)
        sig_file = f.name
    
    try:
        if allowed_signers and allowed_signers.exists():
            # Verify with allowed signers (trusted verification)
            # Need to figure out principal - for now use wildcard
            cmd = [
                "ssh-keygen", "-Y", "verify",
                "-f", str(allowed_signers),
                "-I", "*",  # Match any principal
                "-n", namespace,
                "-s", sig_file
            ]
            result = subprocess.run(cmd, input=content, capture_output=True)
            
            if result.returncode == 0:
                return VerificationStatus.VALID, "verified", None
            else:
                stderr = result.stderr.decode('utf-8')
                if "Could not verify" in stderr:
                    return VerificationStatus.INVALID, "unknown", stderr
                else:
                    return VerificationStatus.VALID_UNTRUSTED, "unknown", "Signer not in allowed_signers"
        else:
            # Without allowed_signers, we can only check signature format validity
            # This is a limitation - can't verify without knowing the public key
            return VerificationStatus.VALID_UNTRUSTED, "unknown", "No allowed_signers file provided"
            
    except FileNotFoundError:
        return VerificationStatus.ERROR, "unknown", "ssh-keygen not found"
    finally:
        Path(content_file).unlink(missing_ok=True)
        Path(sig_file).unlink(missing_ok=True)


def sign_package(
    package_path: str | Path,
    method: str | SignatureMethod = SignatureMethod.GPG,
    key_id: Optional[str] = None,
    key_path: Optional[str | Path] = None,
    add: bool = False,
    signer_name: Optional[str] = None,
) -> SignatureInfo:
    """Sign a yapCAD package.
    
    Args:
        package_path: Path to .ycpkg package directory
        method: "gpg" or "ssh"
        key_id: GPG key ID (for GPG signing)
        key_path: Path to SSH private key (for SSH signing)
        add: If True, add signature without removing existing ones
        signer_name: Override signer name in signature entry
        
    Returns:
        SignatureInfo for the new signature
    """
    pkg_path = Path(package_path)
    manifest_path = pkg_path / "manifest.yaml"
    
    if not manifest_path.exists():
        raise SigningError(f"Package manifest not found: {manifest_path}")
    
    # Convert method to enum
    if isinstance(method, str):
        method = SignatureMethod(method.lower())
    
    # Compute manifest hash
    manifest_hash = _compute_canonical_manifest_hash(manifest_path)
    content = _create_signable_content(manifest_hash)
    
    # Create signatures directory
    sig_dir = pkg_path / "signatures"
    sig_dir.mkdir(exist_ok=True)
    
    # Sign based on method
    timestamp = datetime.now(timezone.utc).isoformat()
    
    if method == SignatureMethod.GPG:
        signature, signer_id, actual_key_id = _sign_with_gpg(content, key_id)
        sig_file = sig_dir / "manifest.sig.asc"
        sig_file.write_text(signature)
        
        sig_info = SignatureInfo(
            signer=signer_name or signer_id,
            method=method,
            key_id=actual_key_id,
            timestamp=timestamp,
            signature_file=f"signatures/{sig_file.name}",
        )
    else:
        if not key_path:
            raise SigningError("SSH signing requires key_path")
        signature, signer_id, fingerprint = _sign_with_ssh(content, Path(key_path))
        
        # Generate unique filename for SSH signatures
        safe_name = signer_id.replace('@', '_').replace(' ', '_').replace('/', '_')[:20]
        sig_file = sig_dir / f"manifest.sig.{safe_name}"
        sig_file.write_text(signature)
        
        sig_info = SignatureInfo(
            signer=signer_name or signer_id,
            method=method,
            key_id=fingerprint,
            timestamp=timestamp,
            signature_file=f"signatures/{sig_file.name}",
        )
    
    # Update manifest
    with open(manifest_path, 'r', encoding='utf-8') as f:
        manifest = yaml.safe_load(f)
    
    sig_entry = {
        "signer": sig_info.signer,
        "method": sig_info.method.value,
        "key_id": sig_info.key_id,
        "timestamp": sig_info.timestamp,
        "signature_file": sig_info.signature_file,
    }
    
    if add and "signatures" in manifest:
        manifest["signatures"].append(sig_entry)
    else:
        manifest["signatures"] = [sig_entry]
    
    with open(manifest_path, 'w', encoding='utf-8') as f:
        yaml.safe_dump(manifest, f, default_flow_style=False, sort_keys=False)
    
    return sig_info


def verify_package(
    package_path: str | Path,
    allowed_signers: Optional[str | Path] = None,
    gpg_keyring: Optional[str | Path] = None,
) -> VerificationResult:
    """Verify signatures on a yapCAD package.
    
    Args:
        package_path: Path to .ycpkg package directory
        allowed_signers: Path to SSH allowed_signers file
        gpg_keyring: Path to custom GPG keyring
        
    Returns:
        VerificationResult with status of all signatures
    """
    pkg_path = Path(package_path)
    manifest_path = pkg_path / "manifest.yaml"
    
    if not manifest_path.exists():
        return VerificationResult(
            package_path=str(pkg_path),
            manifest_hash="",
            errors=["Package manifest not found"]
        )
    
    # Compute expected manifest hash
    manifest_hash = _compute_canonical_manifest_hash(manifest_path)
    content = _create_signable_content(manifest_hash)
    
    # Load manifest
    with open(manifest_path, 'r', encoding='utf-8') as f:
        manifest = yaml.safe_load(f)
    
    result = VerificationResult(
        package_path=str(pkg_path),
        manifest_hash=f"sha256:{manifest_hash}",
    )
    
    signatures = manifest.get("signatures", [])
    if not signatures:
        result.warnings.append("Package has no signatures")
        return result
    
    # Verify each signature
    for sig_entry in signatures:
        method_str = sig_entry.get("method", "gpg")
        try:
            method = SignatureMethod(method_str)
        except ValueError:
            result.errors.append(f"Unknown signature method: {method_str}")
            continue
        
        sig_file_path = pkg_path / sig_entry.get("signature_file", "")
        if not sig_file_path.exists():
            sig_info = SignatureInfo(
                signer=sig_entry.get("signer", "unknown"),
                method=method,
                key_id=sig_entry.get("key_id", "unknown"),
                timestamp=sig_entry.get("timestamp", ""),
                signature_file=sig_entry.get("signature_file", ""),
                status=VerificationStatus.ERROR,
                error_message="Signature file not found"
            )
            result.signatures.append(sig_info)
            continue
        
        signature = sig_file_path.read_text()
        
        if method == SignatureMethod.GPG:
            status, signer, error = _verify_gpg_signature(
                content, signature, 
                keyring=Path(gpg_keyring) if gpg_keyring else None
            )
        else:
            status, signer, error = _verify_ssh_signature(
                content, signature,
                allowed_signers=Path(allowed_signers) if allowed_signers else None
            )
        
        sig_info = SignatureInfo(
            signer=sig_entry.get("signer", signer),
            method=method,
            key_id=sig_entry.get("key_id", "unknown"),
            timestamp=sig_entry.get("timestamp", ""),
            signature_file=sig_entry.get("signature_file", ""),
            status=status,
            error_message=error
        )
        result.signatures.append(sig_info)
    
    return result


def list_signatures(package_path: str | Path) -> List[SignatureInfo]:
    """List all signatures on a package.
    
    Args:
        package_path: Path to .ycpkg package directory
        
    Returns:
        List of SignatureInfo for each signature
    """
    pkg_path = Path(package_path)
    manifest_path = pkg_path / "manifest.yaml"
    
    if not manifest_path.exists():
        return []
    
    with open(manifest_path, 'r', encoding='utf-8') as f:
        manifest = yaml.safe_load(f)
    
    signatures = manifest.get("signatures", [])
    result = []
    
    for sig_entry in signatures:
        try:
            method = SignatureMethod(sig_entry.get("method", "gpg"))
        except ValueError:
            method = SignatureMethod.GPG
            
        sig_info = SignatureInfo(
            signer=sig_entry.get("signer", "unknown"),
            method=method,
            key_id=sig_entry.get("key_id", "unknown"),
            timestamp=sig_entry.get("timestamp", ""),
            signature_file=sig_entry.get("signature_file", ""),
        )
        result.append(sig_info)
    
    return result


__all__ = [
    "SignatureMethod",
    "VerificationStatus", 
    "SignatureInfo",
    "VerificationResult",
    "SigningError",
    "VerificationError",
    "sign_package",
    "verify_package",
    "list_signatures",
]

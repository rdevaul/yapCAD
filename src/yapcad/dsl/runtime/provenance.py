"""
Provenance tracking for DSL execution.

Captures invocation metadata for audit trails and reproducibility.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, Optional
from datetime import datetime, timezone
import hashlib
import json


@dataclass
class Provenance:
    """
    Provenance information for a DSL command execution.

    This metadata is attached to emitted geometry to enable:
    - Reproducibility: Re-run with same parameters
    - Audit trail: Track what generated each piece of geometry
    - Version control: Detect if source changed
    """
    # Command identification
    module_name: str
    command_name: str
    version: str = "0.0.0"

    # Parameters used
    parameters: Dict[str, Any] = field(default_factory=dict)

    # Source signature for change detection
    source_signature: str = ""

    # Execution timestamp
    timestamp: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())

    # Optional extra metadata
    extra: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "invocation": {
                "module": self.module_name,
                "command": self.command_name,
                "version": self.version,
                "parameters": self.parameters,
                "sourceSignature": self.source_signature,
                "timestamp": self.timestamp,
            },
            "extra": self.extra,
        }

    def to_json(self) -> str:
        """Convert to JSON string."""
        return json.dumps(self.to_dict(), indent=2)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Provenance":
        """Create from dictionary."""
        inv = data.get("invocation", {})
        return cls(
            module_name=inv.get("module", ""),
            command_name=inv.get("command", ""),
            version=inv.get("version", "0.0.0"),
            parameters=inv.get("parameters", {}),
            source_signature=inv.get("sourceSignature", ""),
            timestamp=inv.get("timestamp", ""),
            extra=data.get("extra", {}),
        )


def compute_source_signature(source: str) -> str:
    """
    Compute a signature for source code.

    Uses SHA-256 hash of the normalized source.
    """
    # Normalize: strip whitespace, normalize line endings
    normalized = source.strip().replace('\r\n', '\n').replace('\r', '\n')
    return f"sha256:{hashlib.sha256(normalized.encode()).hexdigest()}"


def create_provenance(
    module_name: str,
    command_name: str,
    parameters: Dict[str, Any],
    source: str = "",
    version: str = "0.0.0",
    **extra: Any,
) -> Provenance:
    """
    Create a Provenance record for a command execution.

    Args:
        module_name: Name of the DSL module
        command_name: Name of the command being executed
        parameters: Parameter values passed to the command
        source: Original source code (for signature)
        version: Version of the module
        **extra: Additional metadata to include

    Returns:
        A Provenance record
    """
    # Serialize parameters to JSON-compatible types
    serialized_params = {}
    for key, value in parameters.items():
        if hasattr(value, 'tolist'):  # numpy arrays
            serialized_params[key] = value.tolist()
        elif hasattr(value, '__iter__') and not isinstance(value, (str, dict)):
            serialized_params[key] = list(value)
        else:
            serialized_params[key] = value

    return Provenance(
        module_name=module_name,
        command_name=command_name,
        version=version,
        parameters=serialized_params,
        source_signature=compute_source_signature(source) if source else "",
        extra=dict(extra) if extra else {},
    )


def verify_source_signature(source: str, signature: str) -> bool:
    """
    Verify that source code matches a signature.

    Args:
        source: The source code to verify
        signature: The expected signature (e.g., "sha256:abc123...")

    Returns:
        True if the source matches the signature
    """
    computed = compute_source_signature(source)
    return computed == signature

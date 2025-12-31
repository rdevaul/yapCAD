Package Signing Specification
=============================

**Version:** ``ycpkg-signing-v0.1``
**Status:** Draft – Provisional signing for yapCAD 1.0

1. Overview
-----------

This document specifies a provisional package signing system for yapCAD
packages (``.ycpkg``). The system leverages existing PKI infrastructure
(GPG or SSH keys) to allow individual contributors to sign packages
without requiring a centralized authority.

**Goals:**

- Enable cryptographic verification of package integrity and authorship
- Work with existing developer tooling (GPG, SSH)
- Support offline verification (no network required)
- Provide clear trust semantics for automation pipelines

**Non-Goals (deferred to 1.1+):**

- Multi-signature approval workflows
- Delegation and authority chains
- Centralized key servers or PKI infrastructure
- Revocation lists

2. Signing Methods
------------------

Two signing methods are supported:

2.1 GPG Signatures
~~~~~~~~~~~~~~~~~~

Uses GnuPG for signing. Keys are identified by fingerprint or key ID.
This is the preferred method for organizations with existing GPG infrastructure.

**Requirements:**

- ``gpg`` command-line tool installed
- Signing key in user's GPG keyring
- Public key available for verification

**Signature format:** Detached ASCII-armored (``.asc``)

2.2 SSH Signatures
~~~~~~~~~~~~~~~~~~

Uses SSH keys for signing (OpenSSH 8.0+). More accessible for developers
who already have SSH keys for Git/GitHub.

**Requirements:**

- ``ssh-keygen`` command-line tool (OpenSSH 8.0+)
- SSH private key (ed25519, rsa, ecdsa)
- Public key available for verification

**Signature format:** SSH signature format (``ssh-keygen -Y sign``)

3. What Gets Signed
-------------------

The signature covers a **canonical manifest hash**, not individual files.
This provides:

- Single signature covers entire package
- Manifest already contains hashes of all files
- Efficient verification (hash manifest, verify signature)

3.1 Canonical Manifest Format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before signing, the manifest is canonicalized:

1. Load ``manifest.yaml`` as YAML
2. Remove the ``signatures`` section (if present)
3. Sort all keys alphabetically (recursive)
4. Serialize to JSON with sorted keys, no extra whitespace
5. Compute SHA-256 hash of UTF-8 encoded JSON bytes

This produces a deterministic hash regardless of YAML formatting.

3.2 Signed Content
~~~~~~~~~~~~~~~~~~

The signature is computed over:

.. code-block:: text

   ycpkg-manifest-v1\n
   <sha256_hex_lowercase>\n

The prefix ensures signatures can't be repurposed from other contexts.

4. Manifest Schema
------------------

Signatures are stored in the manifest:

.. code-block:: yaml

   # In manifest.yaml
   signatures:
     - signer: "Alice Developer <alice@example.com>"
       method: gpg
       key_id: "ABCD1234EFGH5678"
       timestamp: "2025-12-30T15:30:00Z"
       signature_file: "signatures/manifest.sig.asc"
       
     - signer: "Bob Engineer <bob@example.com>"
       method: ssh
       key_fingerprint: "SHA256:abcd1234..."
       public_key_file: "signatures/bob.pub"
       timestamp: "2025-12-30T16:00:00Z"
       signature_file: "signatures/manifest.sig.bob"

4.1 Signature Entry Fields
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Required:**

- ``method``: ``gpg`` or ``ssh``
- ``timestamp``: ISO 8601 timestamp of signing
- ``signature_file``: Path to signature file (relative to package root)

**For GPG signatures:**

- ``signer``: User ID from signing key
- ``key_id``: GPG key ID (short or long form)

**For SSH signatures:**

- ``signer``: Identifier (name, email)
- ``key_fingerprint``: SSH key fingerprint (SHA256 format)
- ``public_key_file``: Path to public key file (optional, for self-contained packages)

5. Package Layout
-----------------

.. code-block:: text

   my_design.ycpkg/
   ├── manifest.yaml
   ├── geometry/
   │   └── primary.json
   ├── signatures/                    # Signature files
   │   ├── manifest.sig.asc           # GPG signature
   │   ├── manifest.sig.bob           # SSH signature (Bob)
   │   └── bob.pub                    # Bob's public key (optional)
   └── ...

6. CLI Operations
-----------------

6.1 Signing a Package
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Sign with GPG (uses default key)
   python -m yapcad.package sign my_design.ycpkg/
   
   # Sign with specific GPG key
   python -m yapcad.package sign my_design.ycpkg/ --gpg-key ABCD1234
   
   # Sign with SSH key
   python -m yapcad.package sign my_design.ycpkg/ --ssh-key ~/.ssh/id_ed25519
   
   # Add signature (don't replace existing)
   python -m yapcad.package sign my_design.ycpkg/ --add

6.2 Verifying Signatures
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Verify all signatures
   python -m yapcad.package verify my_design.ycpkg/
   
   # Verify specific signer
   python -m yapcad.package verify my_design.ycpkg/ --signer alice@example.com
   
   # Verify with allowed signers file (SSH)
   python -m yapcad.package verify my_design.ycpkg/ --allowed-signers ~/.ssh/allowed_signers
   
   # Verify with GPG keyring
   python -m yapcad.package verify my_design.ycpkg/ --gpg-keyring custom.gpg

6.3 Listing Signatures
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   python -m yapcad.package signatures my_design.ycpkg/
   
   # Output:
   # Signatures for my_design.ycpkg:
   #   1. Alice Developer <alice@example.com>
   #      Method: gpg, Key: ABCD1234
   #      Signed: 2025-12-30T15:30:00Z
   #      Status: VALID
   #   2. Bob Engineer <bob@example.com>
   #      Method: ssh, Key: SHA256:abcd1234...
   #      Signed: 2025-12-30T16:00:00Z
   #      Status: VALID (key in package)

7. Verification Semantics
-------------------------

7.1 Verification Levels
~~~~~~~~~~~~~~~~~~~~~~~

``VALID``
  Signature cryptographically valid AND signer key is trusted.
  
``VALID_UNTRUSTED``
  Signature cryptographically valid BUT signer key not in trust store.
  The package hasn't been tampered with, but the signer is unknown.
  
``INVALID``
  Signature doesn't match manifest. Package may have been modified.
  
``ERROR``
  Verification failed (missing key, corrupted signature, etc.)

7.2 Trust Model
~~~~~~~~~~~~~~~

For GPG signatures:
  Key must be in GPG keyring with at least marginal trust.
  
For SSH signatures:
  Key must be in ``allowed_signers`` file (same format as git).

7.3 Verification Output
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: json

   {
     "package": "my_design.ycpkg",
     "manifest_hash": "sha256:abcd1234...",
     "signatures": [
       {
         "signer": "Alice Developer <alice@example.com>",
         "method": "gpg",
         "status": "VALID",
         "key_id": "ABCD1234EFGH5678",
         "timestamp": "2025-12-30T15:30:00Z"
       }
     ],
     "overall_status": "VALID",
     "trusted_signers": 1,
     "untrusted_signers": 0
   }

8. API Surface
--------------

.. code-block:: python

   from yapcad.package.signing import (
       sign_package,
       verify_package,
       SignatureMethod,
       VerificationResult,
       SigningError,
   )
   
   # Sign a package
   sign_package(
       package_path="my_design.ycpkg",
       method=SignatureMethod.GPG,
       key_id="ABCD1234",
   )
   
   # Verify signatures
   result: VerificationResult = verify_package(
       package_path="my_design.ycpkg",
       allowed_signers=Path("~/.ssh/allowed_signers"),
   )
   
   if result.is_valid:
       print(f"Package verified: {result.trusted_signers} trusted signatures")
   else:
       print(f"Verification failed: {result.errors}")

9. Security Considerations
--------------------------

9.1 What Signing Provides
~~~~~~~~~~~~~~~~~~~~~~~~~

- **Integrity:** Detect modifications to any file in the package
- **Attribution:** Identify who signed the package
- **Non-repudiation:** Signer cannot deny signing (if key is secured)

9.2 What Signing Does NOT Provide
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **Authorization:** Signing doesn't grant permissions
- **Trust transitivity:** Signing by A doesn't imply trust in A's dependencies
- **Timestamp verification:** Timestamps are self-reported, not authenticated

9.3 Recommendations
~~~~~~~~~~~~~~~~~~~

1. Use hardware-backed keys (YubiKey, etc.) for high-value signatures
2. Verify signer identity out-of-band before trusting a key
3. Rotate signing keys periodically
4. Keep private keys secure (encrypted, limited access)

10. Interoperability
--------------------

10.1 Git Integration
~~~~~~~~~~~~~~~~~~~~

SSH signatures use the same format as Git commit signing:

.. code-block:: bash

   # Same allowed_signers file works for both
   git config --global gpg.ssh.allowedSignersFile ~/.ssh/allowed_signers

10.2 GitHub Verification
~~~~~~~~~~~~~~~~~~~~~~~~

Packages signed with SSH keys registered on GitHub can be verified using
GitHub's ``/users/{username}/keys`` API (future enhancement).

11. Future Extensions (1.1+)
----------------------------

- **Multi-signature requirements:** Require N of M signatures for approval
- **Role-based signing:** Designer, reviewer, approver roles
- **Delegation:** Allow signing authority delegation
- **Revocation:** Handle compromised keys
- **Timestamping:** RFC 3161 trusted timestamps
- **Key discovery:** Automatic key retrieval from keyservers

12. Example Workflow
--------------------

.. code-block:: bash

   # Designer creates and signs package
   python -m yapcad.dsl run design.dsl MAKE_PART --package part.ycpkg
   python -m yapcad.package sign part.ycpkg --ssh-key ~/.ssh/id_ed25519
   
   # Reviewer verifies, runs tests, adds signature
   python -m yapcad.package verify part.ycpkg
   python tools/ycpkg_analyze.py part.ycpkg --all
   python -m yapcad.package sign part.ycpkg --add --ssh-key ~/.ssh/reviewer_key
   
   # CI pipeline verifies before deployment
   python -m yapcad.package verify part.ycpkg --allowed-signers trusted_signers.txt
   # Exit code 0 = valid, non-zero = invalid

13. See Also
------------

- ``docs/ycpkg_spec.rst`` - Package specification
- ``docs/yapCADone.rst`` - yapCAD 1.0 roadmap
- OpenSSH ssh-keygen(1) - SSH signing documentation
- GnuPG gpg(1) - GPG signing documentation

# security

kind: let

source:
```prose
security = session: reviewer
  prompt: "Review for security vulnerabilities"
```

---

## Security Review: PS-TEROS Codebase

**Review Date:** 2026-01-19
**Reviewer:** Claude Code Security Analysis
**Scope:** Full codebase security audit focusing on input validation, injection vulnerabilities, unsafe file operations, and credential exposure risks

---

## Executive Summary

The PS-TEROS codebase demonstrates **generally good security practices** for a scientific computing application. The code operates within the controlled environment of AiiDA workflow management and computational chemistry clusters, which provides inherent security boundaries. However, several **medium-risk issues** were identified that should be addressed, particularly around input validation and subprocess execution.

**Overall Risk Level:** MEDIUM

**Key Findings:**
- ✅ No hardcoded credentials or API keys detected
- ⚠️ Subprocess execution requires careful review (4 instances)
- ⚠️ Path traversal validation needed in file loading functions
- ⚠️ Temporary file handling could be improved
- ✅ No SQL injection vectors (uses AiiDA ORM)
- ✅ No pickle/marshal deserialization vulnerabilities

---

## Critical Findings (Priority 1)

### None Identified

No critical security vulnerabilities were found that require immediate remediation.

---

## High-Risk Findings (Priority 2)

### None Identified

No high-risk vulnerabilities were found.

---

## Medium-Risk Findings (Priority 3)

### 1. Subprocess Execution Without Input Sanitization

**Location:** `teros/core/fukui/tasks.py`

**Lines:** 444-481, 476-481, 699-710, 830-842

**Description:**
The code executes Python subprocess calls to run FukuiGrid wrapper scripts. While the scripts are internal, the delta_n values and file names are constructed from user inputs without explicit validation.

**Vulnerable Code:**
```python
# Line 468-473
args = [
    sys.executable,
    WRAPPER_SCRIPT,
    ftype,
    delta_n_csv,
] + file_names

result = subprocess.run(
    args,
    cwd=str(tmppath),
    capture_output=True,
    text=True,
)
```

**Risk:**
- If an attacker could control `ftype` or `delta_n_csv` values, they could potentially inject shell commands
- The current implementation appears to construct file names like `CHGCAR_{delta_n:.2f}`, which provides some protection through formatting
- Risk is mitigated by AiiDA's controlled execution environment

**Recommendation:**
```python
# Add explicit validation
def validate_fukui_type(ftype: str) -> str:
    """Validate Fukui type is one of allowed values."""
    if ftype not in ('plus', 'minus'):
        raise ValueError(f"Invalid Fukui type: {ftype}. Must be 'plus' or 'minus'")
    return ftype

def validate_delta_n_list(delta_n_list: list) -> list:
    """Validate delta_n values are numeric and within reasonable range."""
    for dn in delta_n_list:
        if not isinstance(dn, (int, float)):
            raise TypeError(f"delta_n must be numeric, got {type(dn)}")
        if not (0.0 <= dn <= 1.0):
            raise ValueError(f"delta_n out of range: {dn}")
    return delta_n_list

# Apply before subprocess.run()
ftype = validate_fukui_type(fukui_type.value)
delta_n_list = validate_delta_n_list(delta_n_values.get_list())
```

**Impact:** MEDIUM - Could allow command injection if AiiDA security boundaries are bypassed

---

### 2. Path Traversal in File Loading

**Location:** Multiple files
- `teros/core/workgraph.py:45`
- `teros/core/surface_energy/workgraph.py:50`
- `teros/core/convergence/workgraph.py:350`
- `teros/core/surface_hydroxylation/tasks.py:137`

**Description:**
The `load_structure_from_file()` function and similar file reading operations accept filepath strings without validating against path traversal attacks (e.g., `../../etc/passwd`).

**Vulnerable Code:**
```python
def load_structure_from_file(filepath: str) -> orm.StructureData:
    """Load structure from a file using ASE."""
    atoms = read(filepath)  # No path validation
    return orm.StructureData(ase=atoms)
```

**Risk:**
- An attacker with access to workflow submission could read arbitrary files on the system
- ASE's `read()` function can handle many file formats, potentially exposing different parsers to malicious input
- Risk is partially mitigated by filesystem permissions and AiiDA's execution context

**Recommendation:**
```python
from pathlib import Path

def validate_filepath(filepath: str, allowed_dirs: list = None) -> Path:
    """
    Validate filepath is safe and optionally within allowed directories.

    Args:
        filepath: Path to validate
        allowed_dirs: Optional list of allowed parent directories

    Returns:
        Validated Path object

    Raises:
        ValueError: If path contains traversal attempts or is outside allowed dirs
    """
    path = Path(filepath).resolve()

    # Check for path traversal
    if '..' in filepath or filepath.startswith('/'):
        raise ValueError(f"Path traversal detected in: {filepath}")

    # Validate against allowed directories if specified
    if allowed_dirs:
        if not any(path.is_relative_to(Path(d).resolve()) for d in allowed_dirs):
            raise ValueError(f"Path {filepath} is outside allowed directories")

    # Check file exists and is readable
    if not path.exists():
        raise FileNotFoundError(f"File not found: {filepath}")
    if not path.is_file():
        raise ValueError(f"Not a file: {filepath}")

    return path

def load_structure_from_file(filepath: str, allowed_dirs: list = None) -> orm.StructureData:
    """Load structure from a validated file path."""
    validated_path = validate_filepath(filepath, allowed_dirs)
    atoms = read(str(validated_path))
    return orm.StructureData(ase=atoms)
```

**Impact:** MEDIUM - Could expose sensitive files if combined with other vulnerabilities

---

### 3. Temporary File Cleanup Not Guaranteed

**Location:** `teros/core/fukui/tasks.py`

**Lines:** 286-308, 458-503, 673-732, 816-860

**Description:**
Multiple functions use `tempfile.TemporaryDirectory()` context managers, which is good practice. However, the `extract_locpot_from_retrieved()` function manually creates temporary files without guaranteed cleanup on exceptions.

**Vulnerable Code:**
```python
# Line 767-774
with tempfile.NamedTemporaryFile(delete=False, suffix='_LOCPOT') as f:
    f.write(locpot_content)
    temp_path = f.name

try:
    return orm.SinglefileData(temp_path, filename='LOCPOT')
finally:
    Path(temp_path).unlink(missing_ok=True)  # Only cleans up on success
```

**Risk:**
- If `orm.SinglefileData()` raises an exception, the temporary file may not be cleaned up
- Over time, this could lead to disk space exhaustion
- Sensitive data (charge densities, potentials) could remain on disk

**Recommendation:**
```python
def extract_locpot_from_retrieved(retrieved: orm.FolderData) -> orm.SinglefileData:
    """Extract LOCPOT file with guaranteed cleanup."""
    try:
        locpot_content = retrieved.get_object_content('LOCPOT', mode='rb')
    except (FileNotFoundError, OSError) as e:
        raise FileNotFoundError(
            f"LOCPOT not found in retrieved folder. "
            f"Ensure 'LOCPOT' is in ADDITIONAL_RETRIEVE_LIST. Error: {e}"
        )

    # Use context manager for guaranteed cleanup
    with tempfile.TemporaryDirectory() as tmpdir:
        temp_path = Path(tmpdir) / 'LOCPOT'
        temp_path.write_bytes(locpot_content)
        # SinglefileData copies the file, so tmpdir cleanup is safe
        return orm.SinglefileData(str(temp_path), filename='LOCPOT')
```

**Impact:** MEDIUM - Information disclosure and resource exhaustion

---

## Low-Risk Findings (Priority 4)

### 4. Environment Variable Usage Without Validation

**Location:** `examples/fukui/02_Electrodes/run_fukui_potential.py:97`

**Description:**
The example script reads `MP_API_KEY` from environment variables without validation. While this is in example code (not core library), it demonstrates a pattern that could be misused.

**Code:**
```python
api_key = os.environ.get('MP_API_KEY')
if api_key:
    mpr = MPRester(api_key)
```

**Recommendation:**
```python
def get_validated_api_key(env_var: str = 'MP_API_KEY', min_length: int = 16) -> str:
    """Get and validate API key from environment."""
    api_key = os.environ.get(env_var)
    if not api_key:
        return None

    # Basic validation
    if len(api_key) < min_length:
        raise ValueError(f"API key too short (minimum {min_length} characters)")
    if not api_key.strip():
        raise ValueError("API key is empty or whitespace")

    return api_key.strip()
```

**Impact:** LOW - Example code only; proper API key handling should be documented

---

### 5. Module Patching System

**Location:** `teros/core/vasp_parallelization/patches/`

**Description:**
The module includes scripts to patch installed packages (`aiida-vasp`, `aiida-workgraph`) by directly modifying files in site-packages. While this is done for legitimate bug fixes, it represents a security concern:

**Files:**
- `patch_aiida_vasp.py`
- `patch_aiida_workgraph.py`

**Risk:**
- Modifications to installed packages could be exploited if an attacker gains write access
- Patches are applied with user privileges but modify shared library code
- No cryptographic verification of patch integrity

**Recommendation:**
1. Document all patches clearly (✅ Already done in CLAUDE.md)
2. Add integrity checks:
```python
import hashlib

EXPECTED_HASH = {
    'vasp.py': 'sha256:abc123...',  # Hash of original file
}

def verify_file_hash(filepath: Path, expected_hash: str) -> bool:
    """Verify file matches expected hash before patching."""
    content = filepath.read_bytes()
    actual_hash = hashlib.sha256(content).hexdigest()
    return f"sha256:{actual_hash}" == expected_hash
```
3. Consider submitting patches upstream to aiida-vasp and aiida-workgraph
4. Add warning messages when applying patches about security implications

**Impact:** LOW - Requires local file system access; well-documented; legitimate use case

---

### 6. No Input Size Limits

**Location:** Multiple calcfunctions processing VASP output

**Description:**
Functions that parse OUTCAR, CHGCAR, and other VASP output files do not impose size limits. Very large files could cause memory exhaustion.

**Example:** `teros/core/fukui/tasks.py:100-104`
```python
chgcar_content = retrieved.get_object_content('CHGCAR')
# No size check before writing
output_path.write_text(chgcar_content)
```

**Recommendation:**
```python
MAX_CHGCAR_SIZE = 1024 * 1024 * 500  # 500 MB

def get_validated_content(folder: orm.FolderData, filename: str,
                         max_size: int = MAX_CHGCAR_SIZE) -> bytes:
    """Get file content with size validation."""
    # Check size first if possible
    try:
        size = folder.get_object_size(filename)
        if size > max_size:
            raise ValueError(
                f"File {filename} too large: {size} bytes (max {max_size})"
            )
    except AttributeError:
        pass  # get_object_size may not be available

    content = folder.get_object_content(filename, mode='rb')

    if len(content) > max_size:
        raise ValueError(
            f"File {filename} exceeds size limit: {len(content)} bytes"
        )

    return content
```

**Impact:** LOW - DoS via memory exhaustion; requires ability to submit jobs

---

## Positive Security Practices Observed

### 1. ✅ No Hardcoded Credentials
- No API keys, passwords, or tokens found in source code
- Example scripts properly use environment variables for API keys
- `.gitignore` includes `.env` files

### 2. ✅ Safe Database Usage
- All database operations go through AiiDA's ORM
- No raw SQL queries found
- No SQL injection vectors identified

### 3. ✅ No Unsafe Deserialization
- No `pickle.loads()` or `marshal.loads()` calls
- JSON used for data serialization (`thermodynamics_data.json`)
- AiiDA nodes used for structured data

### 4. ✅ Proper Use of Context Managers
- Most temporary file operations use `with tempfile.TemporaryDirectory()`
- Files are properly closed after use
- Resources are generally cleaned up

### 5. ✅ Type Annotations
- Extensive use of type hints throughout codebase
- Helps prevent type confusion vulnerabilities
- Improves code maintainability

### 6. ✅ Input Validation in Some Areas
- AiiDA nodes provide type validation (orm.Float, orm.Int, orm.Str)
- Some functions validate ranges and domains
- Error messages are informative without leaking sensitive info

---

## Recommendations Summary

### Immediate Actions (Next Sprint)

1. **Add input validation to subprocess calls** in `fukui/tasks.py`
   - Validate `fukui_type` against whitelist
   - Validate `delta_n` values are numeric and in range
   - Validate file names match expected pattern

2. **Implement path validation** in all file loading functions
   - Add `validate_filepath()` utility function
   - Enforce allowed directories where applicable
   - Check for path traversal attempts

3. **Fix temporary file cleanup** in `extract_locpot_from_retrieved()`
   - Use `TemporaryDirectory` context manager
   - Ensure cleanup on all code paths

### Short-term Actions (1-2 Months)

4. **Add file size limits** to all file reading operations
   - Prevent memory exhaustion attacks
   - Set reasonable limits based on typical file sizes

5. **Document security considerations**
   - Add SECURITY.md to repository
   - Document trusted execution environment assumptions
   - Provide guidelines for users running on shared systems

6. **Submit patches upstream**
   - Work with aiida-vasp and aiida-workgraph maintainers
   - Reduce need for local package modifications

### Long-term Actions (3-6 Months)

7. **Security audit of dependencies**
   - Regular dependency scanning with tools like `safety` or `pip-audit`
   - Keep AiiDA, pymatgen, and other dependencies updated

8. **Add integration tests for security scenarios**
   - Test path traversal prevention
   - Test subprocess call validation
   - Test file size limits

9. **Consider sandboxing**
   - Evaluate running user-submitted structures in containers
   - Document cluster security best practices

---

## Security Assumptions

The security analysis assumes:

1. **Trusted Users:** Users have legitimate access to the AiiDA system
2. **Controlled Environment:** Calculations run on managed HPC clusters
3. **Filesystem Permissions:** Standard Unix file permissions are enforced
4. **Network Security:** Cluster nodes are on protected networks
5. **AiiDA Security:** AiiDA's security model is trusted and up-to-date

These assumptions are **reasonable for academic/research environments** but may not hold for:
- Public-facing web services
- Multi-tenant cloud environments
- Untrusted user scenarios

---

## Conclusion

The PS-TEROS codebase demonstrates **good security hygiene** for a scientific computing application. The identified issues are **manageable** and can be addressed incrementally without major refactoring. The code benefits from operating within AiiDA's controlled environment, which provides inherent security boundaries.

**Key Strengths:**
- No critical vulnerabilities found
- Good use of modern Python security practices
- Well-structured code with clear boundaries
- Comprehensive documentation and provenance tracking

**Priority Actions:**
1. Add subprocess input validation (Medium priority)
2. Implement path traversal protection (Medium priority)
3. Fix temporary file cleanup (Medium priority)
4. Document security model (Low priority)

**Overall Assessment:** The codebase is **suitable for production use in controlled HPC environments** with the recommended security enhancements applied.

---

## Appendix: Security Checklist

| Category | Status | Notes |
|----------|--------|-------|
| SQL Injection | ✅ PASS | Uses AiiDA ORM exclusively |
| Command Injection | ⚠️ REVIEW | Subprocess calls need validation |
| Path Traversal | ⚠️ REVIEW | File loading needs validation |
| XSS | ✅ N/A | No web interface |
| CSRF | ✅ N/A | No web interface |
| Authentication | ✅ PASS | Handled by AiiDA |
| Authorization | ✅ PASS | Handled by cluster/AiiDA |
| Crypto | ✅ PASS | No custom crypto; uses system SSL |
| Secrets Management | ✅ PASS | No hardcoded secrets |
| Logging | ✅ PASS | Appropriate logging without secrets |
| Error Handling | ✅ PASS | Informative without data leakage |
| Input Validation | ⚠️ REVIEW | Needs improvement in some areas |
| Output Encoding | ✅ N/A | No user-facing output |
| File Operations | ⚠️ REVIEW | Temp file cleanup needed |
| Dependencies | ✅ PASS | Using well-maintained libraries |

Legend: ✅ PASS | ⚠️ REVIEW | ❌ FAIL | ➖ N/A

---

*End of Security Review*

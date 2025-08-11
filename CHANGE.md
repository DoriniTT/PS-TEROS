# Changelog

## [v0.1.2] - Enthalpy of Formation Bug Fix & Oxygen Chemical Potential Limits

### Bug Fixes
- Fixed calculation of enthalpy of formation to correct numerical errors in previous versions.
- Corrected the limits for the oxygen chemical potential to ensure physically meaningful bounds are enforced.

---

## [v0.1.1] - Added DEFECT_TYPES File Generation

### New Feature
- Added functionality to analyze stoichiometric deviations in surface slab terminations
- Generated DEFECT_TYPES file documents excess/deficit elements in each termination
- Enables subsequent calculations with charge compensations for specific terminations

### Implementation
- Modified `get_slabs` function in `functions/slabs.py`:
  - Added `analyze_defects` helper function to compare slab/bulk compositions
  - Calculates element excess relative to bulk stoichiometry for each termination
  - Outputs results to DEFECT_TYPES with columns: termination, [elements]

### File Format
DEFECT_TYPES contains:
- termination: Slab identifier (s_0, s_1, etc.)
- One column per element in bulk structure

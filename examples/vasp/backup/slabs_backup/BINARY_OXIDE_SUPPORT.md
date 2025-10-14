# Binary Oxide Support in PS-TEROS

## Summary

PS-TEROS now fully supports **binary oxide systems** (e.g., Ag₂O, CuO, ZnO) in addition to ternary oxides (e.g., Ag₃PO₄, Fe₂WO₆).

## Changes Made

### 1. Core Workgraph Module (`teros/core/workgraph.py`)

**Modified Functions:**
- `core_workgraph()` - Main @task.graph function
- `build_core_workgraph()` - Convenience wrapper
- `build_core_workgraph_with_map()` - Alias wrapper

**Key Changes:**
- Made nonmetal parameters optional (default to None):
  - `nonmetal_name: str = None`
  - `nonmetal_parameters: dict = None`
  - `nonmetal_options: dict = None`
  - `nonmetal_potential_mapping: dict = None`
  
- Reordered function parameters:
  - Required parameters first: `metal_name`, `oxygen_name`
  - Optional nonmetal parameters after oxygen parameters
  
- Added conditional logic for nonmetal relaxation:
  ```python
  if nonmetal_name is not None:
      # Run nonmetal relaxation (for ternary oxides)
  else:
      # Use dummy values (for binary oxides)
  ```

- Updated docstrings to indicate binary oxide support

### 2. Example File (`examples/slabs/slabs_input_relax_ag2o.py`)

**Adapted for Ag₂O Binary Oxide:**
- Changed bulk structure: `ag3po4.cif` → `ag2o.cif`
- Removed all nonmetal (P) references
- Removed parameters:
  - `nonmetal_name`
  - `nonmetal_parameters`
  - `nonmetal_options`
  - `nonmetal_potential_mapping`
  
- Updated potential mappings:
  - Bulk: `{"Ag": "Ag", "O": "O"}`
  - Slab: `{"Ag": "Ag", "O": "O"}`
  
- Updated slab file names:
  - From: `slab_term_0.cif`
  - To: `ag2o_slab_term_0.cif`

### 3. Core Calculation Modules (No Changes Needed!)

The following modules already supported binary oxides:

**`teros/core/hf.py`:**
- `calculate_formation_enthalpy()` automatically detects binary vs ternary oxides
- For binary: ignores nonmetal contribution even if dummy values are provided
- Returns `oxide_type = 'binary'` in output Dict

**`teros/core/thermodynamics.py`:**
- `identify_oxide_type()` identifies binary (1 metal) vs ternary (2 metals)
- `calculate_surface_energy_binary()` computes γ(Δμ_O) for binary oxides
- `calculate_surface_energy_ternary()` computes γ(Δμ_M, Δμ_O) for ternary oxides

## Usage

### For Binary Oxides (e.g., Ag₂O)

```python
from teros.core.workgraph import build_core_workgraph_with_map

wg = build_core_workgraph_with_map(
    structures_dir="path/to/structures",
    bulk_name="ag2o.cif",
    metal_name="Ag.cif",
    oxygen_name="O2.cif",
    # NO nonmetal_name parameter!
    code_label="VASP-VTST-6.4.3@bohr",
    potential_family="PBE",
    bulk_potential_mapping={"Ag": "Ag", "O": "O"},
    metal_potential_mapping={"Ag": "Ag"},
    oxygen_potential_mapping={"O": "O"},
    # NO nonmetal_potential_mapping!
    bulk_parameters={...},
    metal_parameters={...},
    oxygen_parameters={...},
    # NO nonmetal_parameters!
    input_slabs=input_slabs,
    relax_slabs=True,
    compute_thermodynamics=True,
    name="Ag2O_Slabs",
)
```

### For Ternary Oxides (e.g., Ag₃PO₄)

```python
wg = build_core_workgraph_with_map(
    structures_dir="path/to/structures",
    bulk_name="ag3po4.cif",
    metal_name="Ag.cif",
    nonmetal_name="P.cif",  # Required for ternary
    oxygen_name="O2.cif",
    code_label="VASP-VTST-6.4.3@bohr",
    potential_family="PBE",
    bulk_potential_mapping={"Ag": "Ag", "P": "P", "O": "O"},
    metal_potential_mapping={"Ag": "Ag"},
    nonmetal_potential_mapping={"P": "P"},  # Required for ternary
    oxygen_potential_mapping={"O": "O"},
    bulk_parameters={...},
    metal_parameters={...},
    nonmetal_parameters={...},  # Required for ternary
    oxygen_parameters={...},
    input_slabs=input_slabs,
    relax_slabs=True,
    compute_thermodynamics=True,
    name="Ag3PO4_Slabs",
)
```

## Examples

### Binary Oxide Example
- **File:** `examples/slabs/slabs_input_relax_ag2o.py`
- **System:** Ag₂O
- **References:** Ag metal + O₂ molecule

### Ternary Oxide Example
- **File:** `examples/slabs/slabs_input_relax_ag3po4.py`
- **System:** Ag₃PO₄
- **References:** Ag metal + P + O₂ molecule

## Required Structure Files

### For Binary Oxides (Ag₂O example):
- `ag2o.cif` - Bulk Ag₂O structure
- `Ag.cif` - Ag metal reference
- `O2.cif` - O₂ molecule reference
- `ag2o_slab_term_*.cif` - Pre-generated slab structures (optional)

### For Ternary Oxides (Ag₃PO₄ example):
- `ag3po4.cif` - Bulk Ag₃PO₄ structure
- `Ag.cif` - Ag metal reference
- `P.cif` - P reference
- `O2.cif` - O₂ molecule reference
- `slab_term_*.cif` - Pre-generated slab structures (optional)

## Output

For both binary and ternary oxides, the workflow outputs:

1. **Formation Enthalpy:**
   - `formation_enthalpy` (Dict with ΔH_f and oxide_type)

2. **Relaxed Structures:**
   - `bulk_structure`, `metal_structure`, `oxygen_structure`
   - `nonmetal_structure` (only for ternary oxides)

3. **Slab Outputs (if relax_slabs=True):**
   - `slab_structures` - Input/generated slabs
   - `relaxed_slabs` - Relaxed slab structures
   - `slab_energies` - Total energies of relaxed slabs

4. **Surface Energies (if compute_thermodynamics=True):**
   - Binary: γ(Δμ_O) - 1D function of oxygen chemical potential
   - Ternary: γ(Δμ_M, Δμ_O) - 2D function of metal and oxygen chemical potentials

## Testing

After making changes:
```bash
# Clean Python cache
find . -type d -name __pycache__ -exec rm -rf {} +
find . -name "*.pyc" -delete

# Restart AiiDA daemon
verdi daemon restart

# Run binary oxide example
source ~/envs/psteros/bin/activate && python examples/slabs/slabs_input_relax_ag2o.py
```

## Notes

- The nonmetal relaxation task is only created when `nonmetal_name` is provided
- For binary oxides, dummy nonmetal values are passed to `calculate_formation_enthalpy()`, which automatically ignores them
- The `oxide_type` field in the output Dict indicates whether the system is 'binary' or 'ternary'
- Surface energy calculations automatically use the appropriate function (`calculate_surface_energy_binary` or `calculate_surface_energy_ternary`) based on the oxide type

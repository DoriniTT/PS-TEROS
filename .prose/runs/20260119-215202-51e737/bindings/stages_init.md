# stages_init

kind: let

source:
```prose
let stages_init = session: implementer
  context: analysis
  prompt: """Create the stages module structure..."""
```

---

## Summary

Created the `teros/core/stages/__init__.py` module with the following components:

### StageContext Dataclass

A comprehensive dataclass that holds shared state between workflow stages:

**Workflow Configuration:**
- `wg`: WorkGraph instance (optional, may be None before build)
- `resolved_flags`: dict of workflow flags
- `resolved_preset_name`: str

**Code Configuration:**
- `code_label`: str (main code, used as fallback)
- `bulk_code_label`, `slab_code_label`, `metal_code_label`, `oxygen_code_label`, `nonmetal_code_label`, `aimd_code_label`: Optional[str]

**Potential Configuration:**
- `potential_family`: str
- `bulk_potential_mapping`: dict
- `slab_potential_mapping`, `metal_potential_mapping`, `oxygen_potential_mapping`, `nonmetal_potential_mapping`, `aimd_potential_mapping`, `adsorption_potential_mapping`: Optional[dict]

**Parameter Configuration (INCAR):**
- `bulk_parameters`: dict
- `slab_parameters`, `metal_parameters`, `oxygen_parameters`, `nonmetal_parameters`, `aimd_parameters`, `adsorption_parameters`, `bands_parameters`, `slab_bands_parameters`: Optional[dict]

**Options Configuration (Scheduler):**
- `bulk_options`: dict
- `slab_options`, `metal_options`, `oxygen_options`, `nonmetal_options`, `aimd_options`, `adsorption_options`, `bands_options`, `slab_bands_options`: Optional[dict]

**K-points Configuration:**
- `kpoints_spacing`: float
- `slab_kpoints_spacing`, `aimd_kpoints_spacing`, `adsorption_kpoints_spacing`: Optional[float]

**Other Configuration:**
- `clean_workdir`: bool
- `logger`: logging.Logger

**Restart/Input Data:**
- `restart_folders`: Optional[dict] (RemoteData)
- `restart_slabs`: Optional[dict] (StructureData)
- `input_slabs`: Optional[dict] (StructureData)
- `use_input_slabs`: bool

### Helper Methods

The StageContext includes convenience methods for parameter fallback chains:

- `get_effective_code_label(calculation_type)`: Returns the code label with fallback to main code_label
- `get_effective_kpoints_spacing(calculation_type)`: Returns spacing with fallback to main kpoints_spacing
- `get_effective_potential_mapping(calculation_type)`: Returns mapping with fallback chain (specific -> slab -> bulk)
- `get_effective_parameters(calculation_type)`: Returns INCAR parameters with fallback chain
- `get_effective_options(calculation_type)`: Returns scheduler options with fallback chain

### Stage Function Placeholders

All 21 stage functions are defined as placeholders that raise `NotImplementedError`:

**Stage 1:** `resolve_workflow_preset`
**Stage 2:** `validate_required_inputs`, `log_workflow_configuration`
**Stage 3:** `extract_restart_data`, `prepare_input_slabs`
**Stage 4:** `setup_cp2k_files`, `check_cleavage_compatibility`
**Stage 5:** `build_core_workgraph_base`
**Stage 6:** `resolve_slab_parameters`, `add_restart_slab_tasks`, `add_normal_slab_tasks`, `add_thermodynamics_task`, `add_cleavage_task`
**Stage 7:** `add_slab_electronic_properties`
**Stage 8:** `add_bulk_electronic_properties`
**Stage 9:** `resolve_aimd_parameters`, `prepare_fixed_atoms`, `determine_initial_slabs_source`, `add_aimd_stages`
**Stage 10:** `resolve_adsorption_parameters`, `extract_builder_params`, `add_adsorption_energy_task`

### Module Documentation

Comprehensive docstrings explain:
- The stage builder pattern
- Purpose of each stage
- Usage examples for StageContext
- Public API listing

### Verification

The module was tested and imports correctly:
```python
from teros.core.stages import StageContext
ctx = StageContext(code_label='test')
print(ctx.get_effective_code_label('bulk'))  # 'test'
```

### Files Created

- `/home/thiagotd/git/PS-TEROS/teros/core/stages/__init__.py` (560 lines)

### Next Steps

The placeholder functions should be replaced with actual implementations by importing from their respective modules as they are created:
- `preset_resolution.py`
- `input_validation.py`
- `restart_handling.py`
- `cp2k_setup.py`
- `core_build.py`
- `input_slabs.py`
- `slab_electronic.py`
- `bulk_electronic.py`
- `aimd_stage.py`
- `adsorption_stage.py`

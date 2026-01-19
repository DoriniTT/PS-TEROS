# Electronic Properties Stage Extraction

## Status: COMPLETED

## Created File

`/home/thiagotd/git/PS-TEROS/teros/core/stages/electronic_properties.py`

## Extracted Functions

### 1. `add_slab_electronic_properties_stage()`

**Source:** workgraph.py lines 1564-1612

**Parameters:**
- `wg: WorkGraph` - WorkGraph instance to add the task to
- `code_label: str` - AiiDA code label for VASP calculations
- `potential_family: str` - POTCAR family name
- `bulk_options: Dict[str, Any]` - Default scheduler options (fallback)
- `slab_potential_mapping: Optional[Dict[str, str]]` - Element to potential mapping for slabs
- `slab_bands_parameters: Optional[Dict[str, Any]]` - VASP parameters for slab DOS/bands
- `slab_bands_options: Optional[Dict[str, Any]]` - Scheduler options for slab bands
- `slab_band_settings: Optional[Dict[str, Any]]` - Band workflow settings
- `slab_electronic_properties: Dict[str, Dict[str, Any]]` - Per-slab parameter configs
- `max_concurrent_jobs: Optional[int]` - Maximum concurrent VASP calculations
- `use_restart_mode: bool` - Whether restart mode is being used
- `logger: logging.Logger` - Logger instance
- `scatter_task: Optional[Any]` - Scatter task (non-restart mode)
- `collector: Optional[Any]` - Collector task (restart mode)

**Functionality:**
- Adds `calculate_electronic_properties_slabs_scatter` task
- Handles restart vs non-restart sources for relaxed slabs
- Connects outputs: `slab_bands`, `slab_dos`, `slab_primitive_structures`, `slab_seekpath_parameters`
- Logs status with count of slabs enabled

### 2. `add_bulk_electronic_properties_stage()`

**Source:** workgraph.py lines 1614-1696

**Parameters:**
- `wg: WorkGraph` - WorkGraph instance to add the task to
- `code_label: str` - AiiDA code label for VASP calculations
- `potential_family: str` - POTCAR family name
- `bulk_potential_mapping: Dict[str, str]` - Element to potential mapping for bulk
- `bulk_options: Dict[str, Any]` - Default scheduler options
- `bands_parameters: Optional[Dict[str, Any]]` - VASP parameters for DOS/bands (scf, bands, dos keys)
- `bands_options: Optional[Dict[str, Any]]` - Scheduler options for bands
- `band_settings: Optional[Dict[str, Any]]` - Band workflow settings
- `clean_workdir: bool` - Whether to clean remote directories
- `logger: logging.Logger` - Logger instance

**Functionality:**
- Wraps `BandsWorkChain` as a task using `WorkflowFactory('vasp.v2.bands')`
- Gets bulk structure from `wg.tasks['VaspWorkChain'].outputs.structure`
- Builds inputs for SCF, bands, and DOS namespaces based on `bands_parameters`
- Connects outputs: `bulk_bands`, `bulk_dos`, `bulk_primitive_structure`, `bulk_seekpath_parameters`
- Logs status when enabled

## Design Patterns

Both functions follow the established stage pattern:
- Take `wg` (WorkGraph) as first parameter
- Modify workgraph in place (no return value)
- Use `wg.add_task()` to add tasks
- Connect outputs via `wg.outputs.xxx = task.outputs.xxx`
- Log status with `logger.info()`
- Full Google-style docstrings with Args, Returns, Example sections
- Type hints throughout

## Notes

- The slab function requires either `scatter_task` or `collector` based on `use_restart_mode`
- The bulk function assumes a task named 'VaspWorkChain' exists in the workgraph
- Both functions use lazy imports for AiiDA modules

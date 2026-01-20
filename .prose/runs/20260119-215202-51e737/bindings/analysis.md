# analysis

kind: let

source:
```prose
let analysis = session: architect
  prompt: """Analyze the build_core_workgraph() function..."""
```

---

# Refactoring Plan: build_core_workgraph() Decomposition

## Executive Summary

The `build_core_workgraph()` function (lines 638-1974, ~1350 lines) orchestrates 10 distinct stages of workflow construction. This plan proposes extracting each stage into a dedicated module under `teros/core/stages/`, maintaining 100% backward compatibility while improving testability and maintainability.

## Current Architecture Analysis

### Function Signature
- **104 parameters** (lines 638-742)
- **407 lines of docstring** (lines 743-1045)
- **~930 lines of implementation** (lines 1046-1974)

### Identified Stages (with line ranges)

| Stage | Lines | Purpose | Complexity |
|-------|-------|---------|------------|
| 1. Workflow preset resolution | 1046-1131 | Resolve presets, validate inputs | Medium |
| 2. Input validation | 1133-1168 | Validate required inputs, log config | Low |
| 3. Restart handling | 1170-1218 | Extract data from previous workgraph | Medium |
| 4. CP2K setup | 1222-1260 | Create basis/pseudo files for CP2K | Low |
| 5. Core workgraph build | 1262-1346 | Call core_workgraph.build() | Medium |
| 6. Input slabs handling | 1348-1562 | Relax/thermodynamics/cleavage for manual slabs | High |
| 7. Slab electronic properties | 1564-1612 | DOS/bands for slabs | Medium |
| 8. Bulk electronic properties | 1614-1696 | DOS/bands for bulk | Medium |
| 9. AIMD calculation | 1698-1862 | Multi-stage AIMD simulation | High |
| 10. Adsorption energy | 1864-1965 | Adsorption energy workflow | Medium |

## Proposed Module Structure

```
teros/core/stages/
    __init__.py                    # Public API exports
    _types.py                      # Shared type definitions & dataclasses
    preset_resolution.py           # Stage 1: Preset resolution
    input_validation.py            # Stage 2: Input validation
    restart_handling.py            # Stage 3: Restart handling
    cp2k_setup.py                  # Stage 4: CP2K-specific setup
    core_build.py                  # Stage 5: Core workgraph construction
    input_slabs.py                 # Stage 6: Manual/restart slab handling
    slab_electronic.py             # Stage 7: Slab electronic properties
    bulk_electronic.py             # Stage 8: Bulk electronic properties
    aimd_stage.py                  # Stage 9: AIMD calculation
    adsorption_stage.py            # Stage 10: Adsorption energy
```

## Detailed Function Signatures

### Stage 0: Shared Types (`_types.py`)

```python
"""Shared type definitions for workflow stages."""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any
from aiida import orm

@dataclass
class ResolvedFlags:
    """Resolved workflow flags after preset application."""
    relax_slabs: bool
    compute_thermodynamics: bool
    compute_cleavage: bool
    compute_relaxation_energy: bool
    compute_electronic_properties_bulk: bool
    compute_electronic_properties_slabs: bool
    run_aimd: bool
    run_adsorption_energy: bool

@dataclass
class RestartData:
    """Data extracted from a previous workgraph for restart."""
    folders: Optional[Dict[str, orm.RemoteData]] = None
    slabs: Optional[Dict[str, orm.StructureData]] = None
    
    @property
    def has_restart_data(self) -> bool:
        return self.folders is not None or self.slabs is not None

@dataclass
class CP2KFiles:
    """CP2K basis and pseudopotential files."""
    basis_file: Optional[orm.SinglefileData] = None
    pseudo_file: Optional[orm.SinglefileData] = None

@dataclass  
class SlabParameters:
    """Resolved slab calculation parameters."""
    parameters: dict
    options: dict
    potential_mapping: dict
    kpoints_spacing: float
    code: orm.Code

@dataclass
class StageContext:
    """Shared context passed between stages."""
    resolved_preset: str
    flags: ResolvedFlags
    restart_data: RestartData
    cp2k_files: CP2KFiles
    use_input_slabs: bool
    input_slabs: Optional[Dict[str, orm.StructureData]]
    # References to created tasks for inter-stage dependencies
    created_tasks: Dict[str, Any] = field(default_factory=dict)
```

### Stage 1: Preset Resolution (`preset_resolution.py`)

```python
"""Stage 1: Workflow preset resolution and validation."""

from typing import Tuple, Optional
from teros.core.stages._types import ResolvedFlags

def resolve_workflow_preset(
    workflow_preset: Optional[str],
    relax_slabs: Optional[bool],
    compute_thermodynamics: Optional[bool],
    compute_cleavage: Optional[bool],
    compute_relaxation_energy: Optional[bool],
    compute_electronic_properties_bulk: Optional[bool],
    compute_electronic_properties_slabs: Optional[bool],
    run_aimd: Optional[bool],
    run_adsorption_energy: Optional[bool],
    # Validation parameters
    metal_name: Optional[str] = None,
    oxygen_name: Optional[str] = None,
    nonmetal_name: Optional[str] = None,
    miller_indices: Optional[list] = None,
    bands_parameters: Optional[dict] = None,
    bands_options: Optional[dict] = None,
    band_settings: Optional[dict] = None,
    aimd_sequence: Optional[list] = None,
    aimd_parameters: Optional[dict] = None,
    aimd_options: Optional[dict] = None,
    aimd_potential_mapping: Optional[dict] = None,
    aimd_kpoints_spacing: Optional[float] = None,
    slab_bands_parameters: Optional[dict] = None,
    slab_bands_options: Optional[dict] = None,
    slab_band_settings: Optional[dict] = None,
    slab_electronic_properties: Optional[dict] = None,
    adsorption_structures: Optional[dict] = None,
    adsorption_formulas: Optional[dict] = None,
    input_slabs: Optional[dict] = None,
) -> Tuple[str, ResolvedFlags]:
    """
    Resolve workflow preset and apply user overrides.
    
    Performs:
    1. Check for deprecated old-style API usage
    2. Resolve preset and apply user overrides
    3. Validate preset requirements
    4. Check flag dependencies and emit warnings
    
    Args:
        workflow_preset: Name of workflow preset to use
        relax_slabs: Override preset default
        ... (all flag overrides)
        ... (validation parameters)
        
    Returns:
        Tuple of (resolved_preset_name, ResolvedFlags)
        
    Raises:
        ValueError: If preset validation fails or flag dependencies are invalid
    """
    ...
```

### Stage 2: Input Validation (`input_validation.py`)

```python
"""Stage 2: Input validation and configuration logging."""

from typing import Optional
import logging
from teros.core.stages._types import ResolvedFlags

logger = logging.getLogger(__name__)

def validate_required_inputs(
    resolved_preset_name: str,
    flags: ResolvedFlags,
    structures_dir: Optional[str],
    bulk_name: Optional[str],
    metal_name: Optional[str],
    oxygen_name: Optional[str],
    miller_indices: Optional[list],
    input_slabs: Optional[dict],
) -> None:
    """
    Validate that required inputs are provided based on workflow configuration.
    
    Performs:
    1. Determine if bulk structure is needed
    2. Validate structures_dir and bulk_name when needed
    3. Log workflow configuration
    
    Args:
        resolved_preset_name: Name of resolved preset
        flags: Resolved workflow flags
        structures_dir: Directory containing structure files
        bulk_name: Filename of bulk structure
        metal_name: Metal reference filename
        oxygen_name: Oxygen reference filename
        miller_indices: Miller indices for slab generation
        input_slabs: Pre-generated slab structures
        
    Raises:
        ValueError: If required inputs are missing
    """
    ...

def log_workflow_configuration(
    preset_name: str,
    flags: ResolvedFlags,
) -> None:
    """Log the workflow configuration in a formatted manner."""
    ...
```

### Stage 3: Restart Handling (`restart_handling.py`)

```python
"""Stage 3: Restart data extraction from previous workgraph."""

from typing import Optional
import logging
from aiida import orm
from teros.core.stages._types import RestartData

logger = logging.getLogger(__name__)

def extract_restart_data(
    restart_from_node: Optional[int],
) -> RestartData:
    """
    Extract restart folders and slab structures from a previous workgraph.
    
    Performs:
    1. Load previous workgraph node
    2. Extract RemoteData restart folders
    3. Extract relaxed or generated slab structures
    4. Return RestartData with extracted information
    
    Args:
        restart_from_node: PK of previous workgraph to restart from
        
    Returns:
        RestartData containing folders and slabs (or empty RestartData)
        
    Note:
        Returns empty RestartData if restart_from_node is None or extraction fails
    """
    ...

def prepare_input_slabs(
    input_slabs: Optional[dict],
    restart_data: RestartData,
) -> tuple[Optional[dict], bool]:
    """
    Prepare input slabs, potentially overriding with restart data.
    
    Args:
        input_slabs: User-provided input slabs
        restart_data: Extracted restart data
        
    Returns:
        Tuple of (final_input_slabs, use_input_slabs_flag)
    """
    ...
```

### Stage 4: CP2K Setup (`cp2k_setup.py`)

```python
"""Stage 4: CP2K-specific file setup."""

from typing import Optional
import logging
from aiida import orm
from teros.core.stages._types import CP2KFiles

logger = logging.getLogger(__name__)

def setup_cp2k_files(
    calculator: str,
    basis_content: Optional[str] = None,
    pseudo_content: Optional[str] = None,
) -> CP2KFiles:
    """
    Create CP2K basis and pseudopotential files if using CP2K calculator.
    
    Performs:
    1. Check if calculator is 'cp2k'
    2. Load default or provided basis/pseudo content
    3. Create SinglefileData nodes
    
    Args:
        calculator: Calculator type ('vasp' or 'cp2k')
        basis_content: Custom basis set content (optional)
        pseudo_content: Custom pseudopotential content (optional)
        
    Returns:
        CP2KFiles with basis_file and pseudo_file (None for VASP)
    """
    ...

def check_cleavage_compatibility(
    use_input_slabs: bool,
    compute_cleavage: bool,
) -> bool:
    """
    Check and potentially disable cleavage for manual slabs.
    
    Args:
        use_input_slabs: Whether using manual input slabs
        compute_cleavage: Current cleavage flag value
        
    Returns:
        Updated compute_cleavage flag (False if incompatible)
    """
    ...
```

### Stage 5: Core Workgraph Build (`core_build.py`)

```python
"""Stage 5: Core workgraph construction."""

from typing import Optional, Dict, Any
from aiida import orm
from aiida_workgraph import WorkGraph

def build_core_workgraph_base(
    # Structure parameters
    structures_dir: Optional[str],
    bulk_name: Optional[str],
    # Code parameters
    code_label: str,
    bulk_code_label: Optional[str],
    metal_code_label: Optional[str],
    nonmetal_code_label: Optional[str],
    oxygen_code_label: Optional[str],
    slab_code_label: Optional[str],
    # Potential parameters
    potential_family: str,
    bulk_potential_mapping: dict,
    kpoints_spacing: float,
    # Bulk parameters
    bulk_parameters: dict,
    bulk_options: dict,
    clean_workdir: bool,
    # Reference parameters
    metal_name: Optional[str],
    oxygen_name: Optional[str],
    metal_potential_mapping: Optional[dict],
    metal_parameters: Optional[dict],
    metal_options: Optional[dict],
    oxygen_potential_mapping: Optional[dict],
    oxygen_parameters: Optional[dict],
    oxygen_options: Optional[dict],
    nonmetal_name: Optional[str],
    nonmetal_potential_mapping: Optional[dict],
    nonmetal_parameters: Optional[dict],
    nonmetal_options: Optional[dict],
    # Slab generation parameters
    miller_indices: list,
    min_slab_thickness: float,
    min_vacuum_thickness: float,
    slab_parameters: Optional[dict],
    slab_options: Optional[dict],
    slab_potential_mapping: Optional[dict],
    slab_kpoints_spacing: Optional[float],
    slab_relax_builder_inputs: Optional[dict],
    structure_specific_relax_builder_inputs: Optional[dict],
    slab_scf_builder_inputs: Optional[dict],
    structure_specific_scf_builder_inputs: Optional[dict],
    # Slab generation options
    lll_reduce: bool,
    center_slab: bool,
    symmetrize: bool,
    primitive: bool,
    in_unit_planes: bool,
    max_normal_search: Optional[int],
    # Workflow flags
    relax_slabs: bool,
    compute_thermodynamics: bool,
    thermodynamics_sampling: int,
    compute_relaxation_energy: bool,
    use_input_slabs: bool,
    compute_cleavage: bool,
    # AIMD parameters (passed through)
    run_aimd: bool,
    aimd_sequence: Optional[list],
    aimd_parameters: Optional[dict],
    aimd_options: Optional[dict],
    aimd_potential_mapping: Optional[dict],
    aimd_kpoints_spacing: Optional[float],
    # Slab electronic parameters (passed through)
    compute_electronic_properties_slabs: bool,
    slab_electronic_properties: Optional[dict],
    slab_bands_parameters: Optional[dict],
    slab_bands_options: Optional[dict],
    slab_band_settings: Optional[dict],
    # Adsorption parameters (passed through)
    run_adsorption_energy: bool,
    adsorption_structures: Optional[dict],
    adsorption_formulas: Optional[dict],
    adsorption_parameters: Optional[dict],
    adsorption_options: Optional[dict],
    adsorption_potential_mapping: Optional[dict],
    adsorption_kpoints_spacing: Optional[float],
    relax_before_adsorption: bool,
    adsorption_relax_builder_inputs: Optional[dict],
    adsorption_scf_builder_inputs: Optional[dict],
    adsorption_fix_atoms: bool,
    adsorption_fix_type: Optional[str],
    adsorption_fix_thickness: float,
    adsorption_fix_elements: Optional[list],
    adsorption_structure_specific_relax_builder_inputs: Optional[dict],
    adsorption_structure_specific_scf_builder_inputs: Optional[dict],
    adsorption_structure_component_specific_scf_builder_inputs: Optional[dict],
    # Concurrency
    max_concurrent_jobs: Optional[int],
    max_concurrent_jobs_slabs: Optional[int],
) -> WorkGraph:
    """
    Build the core workgraph by calling core_workgraph.build().
    
    This is a thin wrapper that:
    1. Handles dummy values for input_slabs mode
    2. Calls the @task.graph decorated core_workgraph function
    3. Returns the constructed WorkGraph
    
    Args:
        ... (all core_workgraph parameters)
        
    Returns:
        WorkGraph instance from core_workgraph.build()
    """
    ...
```

### Stage 6: Input Slabs Handling (`input_slabs.py`)

```python
"""Stage 6: Manual/restart slab handling with thermodynamics and cleavage."""

from typing import Optional, Dict, Any
import logging
from aiida import orm
from aiida_workgraph import WorkGraph
from teros.core.stages._types import RestartData, SlabParameters

logger = logging.getLogger(__name__)

def resolve_slab_parameters(
    slab_parameters: Optional[dict],
    slab_options: Optional[dict],
    slab_potential_mapping: Optional[dict],
    slab_kpoints_spacing: Optional[float],
    bulk_parameters: dict,
    bulk_options: dict,
    bulk_potential_mapping: dict,
    kpoints_spacing: float,
    code_label: str,
    slab_code_label: Optional[str],
) -> SlabParameters:
    """
    Resolve slab parameters, falling back to bulk parameters when needed.
    
    Args:
        slab_*: Slab-specific parameters (may be None)
        bulk_*: Bulk parameters (used as fallback)
        code_label: Main code label
        slab_code_label: Optional slab-specific code label
        
    Returns:
        SlabParameters with resolved values
    """
    ...

def add_restart_slab_tasks(
    wg: WorkGraph,
    input_slabs: Dict[str, orm.StructureData],
    restart_folders: Dict[str, orm.RemoteData],
    slab_params: SlabParameters,
    potential_family: str,
    clean_workdir: bool,
) -> Dict[str, Any]:
    """
    Add VASP tasks with restart folders for each slab.
    
    Performs:
    1. Create VaspWorkChain task for each slab with restart_folder
    2. Create energy extraction tasks
    3. Create collector task for outputs
    4. Connect outputs to WorkGraph
    
    Args:
        wg: WorkGraph to add tasks to
        input_slabs: Dict of slab structures
        restart_folders: Dict of RemoteData for restart
        slab_params: Resolved slab parameters
        potential_family: Potential family name
        clean_workdir: Whether to clean work directory
        
    Returns:
        Dict with references to created outputs (relaxed_slabs_dict, slab_energies_dict)
    """
    ...

def add_normal_slab_tasks(
    wg: WorkGraph,
    input_slabs: Dict[str, orm.StructureData],
    slab_params: SlabParameters,
    potential_family: str,
    clean_workdir: bool,
    max_concurrent_jobs: Optional[int],
    slab_scf_builder_inputs: Optional[dict],
    structure_specific_scf_builder_inputs: Optional[dict],
    slab_relax_builder_inputs: Optional[dict],
    structure_specific_relax_builder_inputs: Optional[dict],
) -> Dict[str, Any]:
    """
    Add SCF + relax + relaxation energy tasks for input slabs.
    
    Performs:
    1. Add SCF task for unrelaxed energies
    2. Add relaxation task
    3. Add relaxation energy calculation
    4. Connect outputs to WorkGraph
    
    Args:
        ... (similar to add_restart_slab_tasks)
        
    Returns:
        Dict with references to created outputs
    """
    ...

def add_thermodynamics_task(
    wg: WorkGraph,
    relaxed_slabs_dict: Any,
    slab_energies_dict: Any,
    thermodynamics_sampling: int,
) -> None:
    """
    Add surface thermodynamics calculation task.
    
    Performs:
    1. Get bulk and formation enthalpy tasks from workgraph
    2. Add oxide type identification task
    3. Add surface energies scatter task
    4. Connect surface_energies output
    
    Args:
        wg: WorkGraph to add tasks to
        relaxed_slabs_dict: Socket or dict of relaxed slab structures
        slab_energies_dict: Socket or dict of slab energies
        thermodynamics_sampling: Sampling resolution
    """
    ...

def add_cleavage_task(
    wg: WorkGraph,
    relaxed_slabs_dict: Any,
    slab_energies_dict: Any,
) -> None:
    """
    Add cleavage energy calculation task.
    
    Args:
        wg: WorkGraph to add tasks to
        relaxed_slabs_dict: Socket or dict of relaxed slab structures
        slab_energies_dict: Socket or dict of slab energies
    """
    ...
```

### Stage 7: Slab Electronic Properties (`slab_electronic.py`)

```python
"""Stage 7: Slab electronic properties (DOS/bands) calculation."""

from typing import Optional, Dict, Any
import logging
from aiida import orm
from aiida_workgraph import WorkGraph
from teros.core.stages._types import RestartData

logger = logging.getLogger(__name__)

def add_slab_electronic_properties(
    wg: WorkGraph,
    slab_electronic_properties: Dict[str, dict],
    code_label: str,
    potential_family: str,
    bulk_options: dict,
    bulk_potential_mapping: dict,
    slab_parameters: Optional[dict],
    slab_options: Optional[dict],
    slab_potential_mapping: Optional[dict],
    slab_bands_parameters: Optional[dict],
    slab_bands_options: Optional[dict],
    slab_band_settings: Optional[dict],
    clean_workdir: bool,
    max_concurrent_jobs: Optional[int],
    restart_folders: Optional[Dict[str, orm.RemoteData]],
    use_collector: bool,
) -> None:
    """
    Add electronic properties calculation for selected slabs.
    
    Performs:
    1. Resolve default parameters
    2. Determine relaxed slabs source (collector or scatter task)
    3. Add electronic properties scatter task
    4. Connect outputs (slab_bands, slab_dos, etc.)
    
    Args:
        wg: WorkGraph to add tasks to
        slab_electronic_properties: Dict mapping slab labels to parameter overrides
        code_label: VASP code label
        potential_family: Potential family name
        ... (parameter fallbacks and options)
        restart_folders: If not None, indicates restart mode
        use_collector: Whether to use collector task outputs
    """
    ...
```

### Stage 8: Bulk Electronic Properties (`bulk_electronic.py`)

```python
"""Stage 8: Bulk electronic properties (DOS/bands) calculation."""

from typing import Optional, Dict
import logging
from aiida import orm
from aiida_workgraph import WorkGraph

logger = logging.getLogger(__name__)

def add_bulk_electronic_properties(
    wg: WorkGraph,
    code_label: str,
    potential_family: str,
    bulk_potential_mapping: dict,
    bulk_options: dict,
    bands_parameters: Optional[dict],
    bands_options: Optional[dict],
    band_settings: Optional[dict],
    clean_workdir: bool,
) -> None:
    """
    Add electronic properties calculation for relaxed bulk structure.
    
    Performs:
    1. Get BandsWorkChain and wrap as task
    2. Get bulk VASP task from workgraph
    3. Build SCF namespace inputs
    4. Build Bands namespace inputs (if provided)
    5. Build DOS namespace inputs (if provided)
    6. Add BandsWorkChain task
    7. Connect outputs (bulk_bands, bulk_dos, etc.)
    
    Args:
        wg: WorkGraph to add tasks to
        code_label: VASP code label
        potential_family: Potential family name
        bulk_potential_mapping: Element to potential mapping
        bulk_options: Scheduler options (fallback for bands)
        bands_parameters: Dict with 'scf', 'bands', 'dos' INCAR parameters
        bands_options: Scheduler options for bands (overrides bulk_options)
        band_settings: Band workflow settings
        clean_workdir: Whether to clean work directory
    """
    ...
```

### Stage 9: AIMD Calculation (`aimd_stage.py`)

```python
"""Stage 9: Ab initio molecular dynamics (AIMD) calculation."""

from typing import Optional, Dict, List, Any
import logging
from aiida import orm
from aiida_workgraph import WorkGraph
from teros.core.stages._types import CP2KFiles

logger = logging.getLogger(__name__)

def resolve_aimd_parameters(
    aimd_parameters: Optional[dict],
    aimd_options: Optional[dict],
    aimd_potential_mapping: Optional[dict],
    aimd_kpoints_spacing: Optional[float],
    slab_parameters: Optional[dict],
    slab_options: Optional[dict],
    slab_potential_mapping: Optional[dict],
    slab_kpoints_spacing: Optional[float],
    bulk_parameters: dict,
    bulk_options: dict,
    bulk_potential_mapping: dict,
    kpoints_spacing: float,
) -> Dict[str, Any]:
    """
    Resolve AIMD parameters with fallback chain.
    
    Fallback order: aimd_* -> slab_* -> bulk_*
    
    Returns:
        Dict with keys: parameters, options, potential_mapping, kpoints_spacing
    """
    ...

def prepare_fixed_atoms(
    input_slabs: Optional[Dict[str, orm.StructureData]],
    fix_atoms: bool,
    fix_type: Optional[str],
    fix_thickness: float,
    fix_elements: Optional[List[str]],
) -> Dict[str, List[int]]:
    """
    Prepare fixed atoms lists for each slab structure.
    
    Args:
        input_slabs: Dict of slab structures
        fix_atoms: Whether fixing is enabled
        fix_type: Type of constraint ('bottom', 'top', etc.)
        fix_thickness: Thickness for fixing
        fix_elements: Elements to fix
        
    Returns:
        Dict mapping slab labels to lists of fixed atom indices
    """
    ...

def determine_initial_slabs_source(
    wg: WorkGraph,
    relax_slabs: bool,
    input_slabs: Optional[Dict[str, orm.StructureData]],
    restart_folders: Optional[dict],
) -> Any:
    """
    Determine the source of initial structures for AIMD.
    
    Priority:
    1. relax_slabs_scatter task outputs (if relax_slabs and task exists)
    2. collect_slab_outputs_restart task outputs (if restart mode)
    3. input_slabs directly (if provided)
    4. generate_slab_structures task outputs (fallback)
    
    Returns:
        Socket or dict representing initial structures
    """
    ...

def add_aimd_stages(
    wg: WorkGraph,
    aimd_sequence: List[dict],
    calculator: str,
    code_label: str,
    aimd_code_label: Optional[str],
    aimd_params: Dict[str, Any],
    potential_family: str,
    cp2k_files: CP2KFiles,
    fixed_atoms_lists: Dict[str, List[int]],
    fix_atoms: bool,
    fix_type: Optional[str],
    fix_thickness: float,
    fix_elements: Optional[List[str]],
    fix_components: str,
    initial_slabs_source: Any,
    aimd_supercell: Optional[List[int]],
    clean_workdir: bool,
    max_concurrent_jobs: Optional[int],
) -> None:
    """
    Add sequential AIMD stages to the workgraph.
    
    Performs:
    1. Select appropriate scatter function (VASP or CP2K)
    2. Optionally add supercell creation task
    3. Add sequential AIMD stage tasks
    4. Wire outputs between stages
    
    Args:
        wg: WorkGraph to add tasks to
        aimd_sequence: List of stage configurations
        calculator: 'vasp' or 'cp2k'
        ... (all AIMD-specific parameters)
    """
    ...
```

### Stage 10: Adsorption Energy (`adsorption_stage.py`)

```python
"""Stage 10: Adsorption energy calculation."""

from typing import Optional, Dict
import logging
from aiida import orm
from aiida_workgraph import WorkGraph

logger = logging.getLogger(__name__)

def resolve_adsorption_parameters(
    adsorption_parameters: Optional[dict],
    adsorption_options: Optional[dict],
    adsorption_potential_mapping: Optional[dict],
    adsorption_kpoints_spacing: Optional[float],
    slab_parameters: Optional[dict],
    slab_options: Optional[dict],
    slab_potential_mapping: Optional[dict],
    slab_kpoints_spacing: Optional[float],
    bulk_parameters: dict,
    bulk_options: dict,
    bulk_potential_mapping: dict,
    kpoints_spacing: float,
) -> Dict[str, any]:
    """
    Resolve adsorption parameters with fallback chain.
    
    Fallback order: adsorption_* -> slab_* -> bulk_*
    
    Returns:
        Dict with keys: parameters, options, potential_mapping, kpoints_spacing
    """
    ...

def extract_builder_params(
    builder_inputs: Optional[dict],
    fallback_params: Optional[dict],
) -> Optional[dict]:
    """
    Extract INCAR parameters from builder_inputs or use fallback.
    
    Handles both new-style (builder_inputs with 'parameters'/'incar') 
    and old-style (direct parameters dict).
    
    Returns:
        Dict of INCAR parameters or None
    """
    ...

def add_adsorption_energy_task(
    wg: WorkGraph,
    adsorption_structures: Dict[str, orm.StructureData],
    adsorption_formulas: Dict[str, str],
    code_label: str,
    potential_family: str,
    ads_params: Dict[str, any],
    relax_before_adsorption: bool,
    relax_params: Optional[dict],
    scf_params: Optional[dict],
    adsorption_relax_builder_inputs: Optional[dict],
    adsorption_scf_builder_inputs: Optional[dict],
    adsorption_structure_specific_relax_builder_inputs: Optional[dict],
    adsorption_structure_specific_scf_builder_inputs: Optional[dict],
    adsorption_structure_component_specific_scf_builder_inputs: Optional[dict],
    adsorption_fix_atoms: bool,
    adsorption_fix_type: Optional[str],
    adsorption_fix_thickness: float,
    adsorption_fix_elements: Optional[list],
    clean_workdir: bool,
    max_concurrent_jobs: Optional[int],
) -> None:
    """
    Add adsorption energy calculation task.
    
    Performs:
    1. Add compute_adsorption_energies_scatter task
    2. Connect all outputs (relaxed_complete_structures, separated_structures,
       substrate_energies, molecule_energies, complete_energies, adsorption_energies)
    
    Args:
        wg: WorkGraph to add tasks to
        ... (all adsorption-specific parameters)
    """
    ...
```

### Public API (`__init__.py`)

```python
"""
PS-TEROS Workflow Stages Module

This module provides composable stage functions for building the core workgraph.
Each stage is independently testable and maintains clear responsibilities.

Public API:
    resolve_workflow_preset: Stage 1 - Preset resolution
    validate_required_inputs: Stage 2 - Input validation
    extract_restart_data: Stage 3 - Restart handling
    setup_cp2k_files: Stage 4 - CP2K setup
    build_core_workgraph_base: Stage 5 - Core workgraph construction
    add_input_slab_tasks: Stage 6 - Input slabs handling
    add_slab_electronic_properties: Stage 7 - Slab electronic properties
    add_bulk_electronic_properties: Stage 8 - Bulk electronic properties
    add_aimd_stages: Stage 9 - AIMD calculation
    add_adsorption_energy_task: Stage 10 - Adsorption energy
    
Type definitions:
    ResolvedFlags: Dataclass for resolved workflow flags
    RestartData: Dataclass for restart information
    CP2KFiles: Dataclass for CP2K files
    SlabParameters: Dataclass for slab parameters
    StageContext: Shared context between stages
"""

from ._types import (
    ResolvedFlags,
    RestartData,
    CP2KFiles,
    SlabParameters,
    StageContext,
)

from .preset_resolution import resolve_workflow_preset
from .input_validation import validate_required_inputs, log_workflow_configuration
from .restart_handling import extract_restart_data, prepare_input_slabs
from .cp2k_setup import setup_cp2k_files, check_cleavage_compatibility
from .core_build import build_core_workgraph_base
from .input_slabs import (
    resolve_slab_parameters,
    add_restart_slab_tasks,
    add_normal_slab_tasks,
    add_thermodynamics_task,
    add_cleavage_task,
)
from .slab_electronic import add_slab_electronic_properties
from .bulk_electronic import add_bulk_electronic_properties
from .aimd_stage import (
    resolve_aimd_parameters,
    prepare_fixed_atoms,
    determine_initial_slabs_source,
    add_aimd_stages,
)
from .adsorption_stage import (
    resolve_adsorption_parameters,
    extract_builder_params,
    add_adsorption_energy_task,
)

__all__ = [
    # Types
    'ResolvedFlags',
    'RestartData',
    'CP2KFiles',
    'SlabParameters',
    'StageContext',
    # Stage 1
    'resolve_workflow_preset',
    # Stage 2
    'validate_required_inputs',
    'log_workflow_configuration',
    # Stage 3
    'extract_restart_data',
    'prepare_input_slabs',
    # Stage 4
    'setup_cp2k_files',
    'check_cleavage_compatibility',
    # Stage 5
    'build_core_workgraph_base',
    # Stage 6
    'resolve_slab_parameters',
    'add_restart_slab_tasks',
    'add_normal_slab_tasks',
    'add_thermodynamics_task',
    'add_cleavage_task',
    # Stage 7
    'add_slab_electronic_properties',
    # Stage 8
    'add_bulk_electronic_properties',
    # Stage 9
    'resolve_aimd_parameters',
    'prepare_fixed_atoms',
    'determine_initial_slabs_source',
    'add_aimd_stages',
    # Stage 10
    'resolve_adsorption_parameters',
    'extract_builder_params',
    'add_adsorption_energy_task',
]
```

## Refactored build_core_workgraph()

After refactoring, the main function becomes a clean orchestrator (~150 lines):

```python
def build_core_workgraph(
    # ... same 104 parameters ...
) -> WorkGraph:
    """Build a centralized WorkGraph for oxide surface calculations.
    
    [Same docstring, but now refers to stages module for implementation details]
    """
    from teros.core.stages import (
        resolve_workflow_preset,
        validate_required_inputs,
        extract_restart_data,
        prepare_input_slabs,
        setup_cp2k_files,
        check_cleavage_compatibility,
        build_core_workgraph_base,
        resolve_slab_parameters,
        add_restart_slab_tasks,
        add_normal_slab_tasks,
        add_thermodynamics_task,
        add_cleavage_task,
        add_slab_electronic_properties,
        add_bulk_electronic_properties,
        add_aimd_stages,
        resolve_aimd_parameters,
        prepare_fixed_atoms,
        determine_initial_slabs_source,
        resolve_adsorption_parameters,
        extract_builder_params,
        add_adsorption_energy_task,
    )
    
    # Stage 1: Resolve workflow preset
    resolved_preset_name, flags = resolve_workflow_preset(
        workflow_preset, relax_slabs, compute_thermodynamics, ...
    )
    
    # Stage 2: Validate inputs
    validate_required_inputs(resolved_preset_name, flags, structures_dir, bulk_name, ...)
    
    # Stage 3: Handle restart
    restart_data = extract_restart_data(restart_from_node)
    input_slabs, use_input_slabs = prepare_input_slabs(input_slabs, restart_data)
    
    # Stage 4: CP2K setup
    cp2k_files = setup_cp2k_files(calculator, basis_content, pseudo_content)
    compute_cleavage = check_cleavage_compatibility(use_input_slabs, flags.compute_cleavage)
    
    # Stage 5: Build core workgraph
    wg = build_core_workgraph_base(
        structures_dir, bulk_name, code_label, ...
    )
    
    # Stage 6: Input slabs handling
    if use_input_slabs and flags.relax_slabs:
        slab_params = resolve_slab_parameters(...)
        if restart_data.folders:
            outputs = add_restart_slab_tasks(wg, input_slabs, restart_data.folders, ...)
        else:
            outputs = add_normal_slab_tasks(wg, input_slabs, slab_params, ...)
        
        if flags.compute_thermodynamics:
            add_thermodynamics_task(wg, outputs['relaxed_slabs'], outputs['energies'], ...)
        
        if compute_cleavage:
            add_cleavage_task(wg, outputs['relaxed_slabs'], outputs['energies'])
    
    # Stage 7: Slab electronic properties
    if flags.compute_electronic_properties_slabs and flags.relax_slabs and slab_electronic_properties:
        add_slab_electronic_properties(wg, slab_electronic_properties, ...)
    
    # Stage 8: Bulk electronic properties
    if flags.compute_electronic_properties_bulk:
        add_bulk_electronic_properties(wg, code_label, potential_family, ...)
    
    # Stage 9: AIMD
    if flags.run_aimd and aimd_sequence and (input_slabs or miller_indices):
        aimd_params = resolve_aimd_parameters(...)
        fixed_atoms = prepare_fixed_atoms(input_slabs, fix_atoms, ...)
        initial_source = determine_initial_slabs_source(wg, flags.relax_slabs, ...)
        add_aimd_stages(wg, aimd_sequence, calculator, ...)
    
    # Stage 10: Adsorption energy
    if flags.run_adsorption_energy and adsorption_structures and adsorption_formulas:
        ads_params = resolve_adsorption_parameters(...)
        relax_params = extract_builder_params(adsorption_relax_builder_inputs, ...)
        scf_params = extract_builder_params(adsorption_scf_builder_inputs, ...)
        add_adsorption_energy_task(wg, adsorption_structures, adsorption_formulas, ...)
    
    # Finalize
    wg.name = name
    if max_concurrent_jobs is not None:
        wg.max_number_jobs = max_concurrent_jobs
    
    return wg
```

## Implementation Strategy

### Phase 1: Create Module Structure (Week 1)
1. Create `teros/core/stages/` directory
2. Implement `_types.py` with dataclasses
3. Implement `__init__.py` with imports

### Phase 2: Extract Helper Functions (Week 2)
1. Extract parameter resolution functions (stages 6, 9, 10)
2. Extract validation functions (stages 1, 2)
3. Add unit tests for each helper

### Phase 3: Extract Stage Functions (Weeks 3-4)
1. Extract stages 1-4 (simpler, fewer dependencies)
2. Extract stages 5-6 (core build and input slabs)
3. Extract stages 7-10 (optional calculations)

### Phase 4: Integration and Testing (Week 5)
1. Refactor `build_core_workgraph()` to use stages
2. Run full test suite
3. Verify backward compatibility with existing examples

## Testing Strategy

### Unit Tests (per stage)
```python
# tests/core/stages/test_preset_resolution.py
def test_resolve_preset_default():
    """Test default preset resolution."""
    preset, flags = resolve_workflow_preset(None, None, None, ...)
    assert preset == 'surface_thermodynamics'
    assert flags.relax_slabs == True

def test_resolve_preset_override():
    """Test preset with user override."""
    preset, flags = resolve_workflow_preset(
        'surface_thermodynamics', 
        relax_slabs=False,  # Override
        ...
    )
    assert flags.relax_slabs == False

def test_resolve_preset_invalid():
    """Test invalid preset name raises error."""
    with pytest.raises(ValueError):
        resolve_workflow_preset('invalid_preset', ...)
```

### Integration Tests
```python
# tests/core/stages/test_integration.py
def test_full_workflow_equivalent():
    """Verify refactored function produces identical workgraph."""
    # Build with old function
    wg_old = build_core_workgraph_original(...)
    
    # Build with refactored function
    wg_new = build_core_workgraph(...)
    
    # Compare task names, connections, outputs
    assert set(wg_old.tasks.keys()) == set(wg_new.tasks.keys())
    assert set(wg_old.outputs.keys()) == set(wg_new.outputs.keys())
```

## Benefits of Refactoring

1. **Testability**: Each stage function can be unit tested in isolation
2. **Readability**: Main function becomes a clear orchestrator (~150 lines vs ~1350)
3. **Maintainability**: Changes to one stage don't affect others
4. **Reusability**: Stage functions can be used by other workflows
5. **Documentation**: Each stage has focused docstrings
6. **Type Safety**: Dataclasses provide clear contracts between stages

## Risk Mitigation

1. **Backward Compatibility**: Keep same function signature and behavior
2. **Gradual Migration**: Implement one stage at a time
3. **Comprehensive Testing**: Unit + integration tests at each step
4. **Feature Flags**: Can toggle between old/new implementations during transition

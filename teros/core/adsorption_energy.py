# teros/core/adsorption_energy.py
"""
Adsorption Energy Calculations Module

This module provides functions to calculate adsorption energies from
substrate+adsorbate structures using AiiDA-WorkGraph and VASP.
"""

from __future__ import annotations

import copy
import typing as t
from collections import Counter

import numpy as np
from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task, namespace, dynamic
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.core import Composition

from .slabs import extract_total_energy, get_settings
from .fixed_atoms import get_fixed_atoms_list, add_fixed_atoms_to_vasp_parameters
from .utils import deep_merge_dicts, extract_max_jobs_value


def _build_vasp_inputs(
    structure: orm.StructureData,
    code: orm.Code,
    builder_inputs: dict | None = None,
    parameters: dict | None = None,
    options: dict | None = None,
    potential_family: str | None = None,
    potential_mapping: dict | None = None,
    kpoints_spacing: float | None = None,
    clean_workdir: bool = True,
    force_scf: bool = False,
    workchain_type: str = 'vasp',
    # Atom fixing parameters
    fix_atoms: bool = False,
    fix_type: str | None = None,
    fix_thickness: float = 0.0,
    fix_elements: list[str] | None = None,
) -> dict:
    """
    Build VASP WorkChain inputs from either builder_inputs or old-style parameters.

    This helper provides backward compatibility by accepting either:
    1. New-style: builder_inputs dict (full control, recommended)
    2. Old-style: parameters/options/etc dicts (automatic conversion)

    **Priority Behavior:**
    If both builder_inputs and parameters are provided, builder_inputs takes priority
    and parameters is ignored. This allows explicit control over the input construction.

    **WorkChain Type:**
    - 'vasp': vasp.v2.vasp plugin (flat structure, SCF calculations)
    - 'relax': vasp.v2.relax plugin (nested structure with 'vasp' namespace and 'relax_settings')
      * When using builder_inputs with 'relax', provide full structure:
        {'relax_settings': orm.Dict(...), 'vasp': {...}}
      * When using parameters with 'relax', relax_settings is auto-generated from INCAR

    **Atom Fixing:**
    Use fix_atoms=True with fix_type/fix_thickness to constrain atoms during relaxation.
    This is useful for slab calculations where you want to fix bottom layers.

    Example showing priority:
        >>> # If both are provided, builder_inputs is used
        >>> inputs = _build_vasp_inputs(
        ...     structure=struct,
        ...     code=code,
        ...     builder_inputs={'parameters': {'incar': {'ENCUT': 600}}},
        ...     parameters={'ENCUT': 400},  # This is ignored
        ... )
        >>> inputs['parameters']['incar']['ENCUT']  # Returns 600

    Example with atom fixing:
        >>> # Fix bottom 7 Å of slab during relaxation
        >>> inputs = _build_vasp_inputs(
        ...     structure=slab,
        ...     code=code,
        ...     parameters={'NSW': 100, 'IBRION': 2},
        ...     workchain_type='relax',
        ...     fix_atoms=True,
        ...     fix_type='bottom',
        ...     fix_thickness=7.0,
        ... )

    Args:
        structure: Structure to calculate
        code: AiiDA code for VASP
        builder_inputs: New-style builder dict (takes priority if provided)
        parameters: Old-style INCAR dict (fallback)
        options: Old-style scheduler dict (fallback)
        potential_family: Pseudopotential family name
        potential_mapping: Element → potential mapping
        kpoints_spacing: K-points spacing in Angstrom^-1
        clean_workdir: Whether to clean remote working directory
        force_scf: If True, enforce NSW=0 and IBRION=-1 for single-point calculation
        workchain_type: Type of workchain ('vasp' or 'relax')
        fix_atoms: If True, apply atom fixing constraints
        fix_type: Type of fixing ('bottom', 'top', 'center')
        fix_thickness: Thickness in Angstrom for fixing region
        fix_elements: Optional list of elements to fix (None = all elements)

    Returns:
        Complete input dictionary for VASP WorkChain

    Raises:
        ValueError: If neither builder_inputs nor parameters is provided
    """
    if builder_inputs is not None:
        # Use new-style builder inputs
        # Make a deep copy to avoid modifying original
        inputs = copy.deepcopy(builder_inputs)

        # Override structure (always set by workflow)
        inputs['structure'] = structure

        # Set code in appropriate location based on workchain type
        if workchain_type == 'relax' and 'vasp' in inputs:
            # For relax workchain with vasp namespace, code goes inside vasp
            inputs['vasp']['code'] = code
        else:
            # For vasp workchain or when no vasp namespace, code goes at top level
            inputs['code'] = code

    elif parameters is not None:
        # Construct from old-style parameters
        inputs = {
            'structure': structure,
            'code': code,
            'parameters': {'incar': dict(parameters)},
            'options': dict(options) if options is not None else {},
            'potential_family': potential_family,
            'potential_mapping': dict(potential_mapping) if potential_mapping is not None else {},
            'clean_workdir': clean_workdir,
            'settings': orm.Dict(dict={
                'parser_settings': {
                    'add_trajectory': True,
                    'add_structure': True,
                    'add_kpoints': True,
                }
            }),
        }

        if kpoints_spacing is not None:
            inputs['kpoints_spacing'] = kpoints_spacing

    else:
        raise ValueError(
            "Must provide either 'builder_inputs' or 'parameters'. "
            "Use builder_inputs for full control (recommended) or parameters for backward compatibility."
        )

    # Apply atom fixing constraints if requested (before restructuring for relax)
    if fix_atoms and fix_type is not None and fix_thickness > 0.0:
        # Get list of atoms to fix
        fixed_atoms_list = get_fixed_atoms_list(
            structure=structure,
            fix_type=fix_type,
            fix_thickness=fix_thickness,
            fix_elements=fix_elements,
        )

        if fixed_atoms_list:
            # Determine where parameters are located
            # For relax with builder_inputs, params are in vasp namespace
            if 'vasp' in inputs and 'parameters' in inputs['vasp']:
                params_location = inputs['vasp']['parameters']
            elif 'parameters' in inputs:
                params_location = inputs['parameters']
            else:
                # Create at top level by default
                inputs['parameters'] = {'incar': {}}
                params_location = inputs['parameters']

            # Ensure incar dict exists
            if 'incar' not in params_location:
                params_location['incar'] = {}

            # Add constraints to structure and update parameters
            params_location['incar'], constrained_structure = add_fixed_atoms_to_vasp_parameters(
                base_parameters=params_location['incar'],
                structure=inputs['structure'],
                fixed_atoms_list=fixed_atoms_list,
            )

            # Replace structure with constrained version
            inputs['structure'] = constrained_structure

    # Force SCF mode if requested
    if force_scf:
        # Determine where parameters are located
        if 'vasp' in inputs and 'parameters' in inputs['vasp']:
            params_location = inputs['vasp']['parameters']
        elif 'parameters' in inputs:
            params_location = inputs['parameters']
        else:
            inputs['parameters'] = {'incar': {}}
            params_location = inputs['parameters']

        # Ensure incar dict exists
        if 'incar' not in params_location:
            params_location['incar'] = {}

        # Override NSW and IBRION for single-point calculation
        params_location['incar']['NSW'] = 0
        params_location['incar']['IBRION'] = -1

    # Restructure for VaspRelaxWorkChain if needed
    if workchain_type == 'relax':
        # Check if user already provided relax_settings (full builder control)
        if 'relax_settings' in inputs:
            # User provided full builder structure
            # Convert relax_settings to orm.Dict if it's a plain dict
            if isinstance(inputs['relax_settings'], dict) and not isinstance(inputs['relax_settings'], orm.Dict):
                inputs['relax_settings'] = orm.Dict(dict=inputs['relax_settings'])

            # Ensure structure is at top level
            if 'structure' not in inputs:
                raise ValueError(
                    "When providing relax_settings in builder_inputs, "
                    "structure must be included at top level"
                )
            # Structure is already set correctly above
        else:
            # Auto-generate relax_settings from INCAR parameters (backward compatibility)
            incar = inputs.get('parameters', {}).get('incar', {})

            # Build relax_settings from INCAR parameters
            relax_settings = {
                'algo': 'cg',  # Default: conjugate gradient
                'force_cutoff': abs(incar.get('EDIFFG', 0.05)),  # Default 0.05 eV/Å
                'steps': incar.get('NSW', 100),  # Default 100 steps
                'positions': True,  # Always relax positions
                'shape': incar.get('ISIF', 2) in [3, 4, 5, 6, 7],  # ISIF 3-7 include shape
                'volume': incar.get('ISIF', 2) in [3, 7],  # ISIF 3,7 include volume
                'convergence_on': True,  # Enable convergence checks
                'convergence_absolute': False,
                'convergence_max_iterations': 5,
                'convergence_positions': 0.01,  # 0.01 Å
                'convergence_volume': 0.01,  # 0.01 Å³
                'convergence_shape_lengths': 0.1,  # 0.1 Å
                'convergence_shape_angles': 0.1,  # 0.1 degrees
                'perform': True,  # Actually perform the relaxation
            }

            # Restructure: move everything except 'structure' into 'vasp' namespace
            structure_to_use = inputs.pop('structure')
            vasp_inputs = inputs  # All remaining keys go into vasp namespace

            # Build final nested structure
            inputs = {
                'structure': structure_to_use,
                'relax_settings': orm.Dict(dict=relax_settings),
                'vasp': vasp_inputs,
            }

    return inputs


def parse_formula(formula_str: str) -> dict[str, int]:
    """
    Parse chemical formula string into element counts.

    Args:
        formula_str: Chemical formula (e.g., "OOH", "H2O")

    Returns:
        Dictionary mapping element symbols to counts

    Example:
        >>> parse_formula("OOH")
        {'O': 2, 'H': 1}
    """
    comp = Composition(formula_str)
    return {str(el): int(count) for el, count in comp.items()}


def _separate_adsorbate_structure_impl(
    structure: orm.StructureData,
    adsorbate_formula: orm.Str,
    bond_cutoff: float = 1.8,
) -> dict:
    """
    Internal implementation for separating substrate+adsorbate structure.

    Uses CrystalNN connectivity analysis first, with a distance-based fallback
    for bonds that CrystalNN doesn't detect (e.g., O-O peroxide bonds in OOH).

    Args:
        structure: Complete structure (substrate + adsorbate)
        adsorbate_formula: Chemical formula of adsorbate (e.g., "OOH", "OH")
        bond_cutoff: Maximum bond distance in Angstrom for fallback detection (default: 1.8)

    Returns:
        Dictionary with three StructureData nodes:
        - substrate: Structure with adsorbate removed, original cell preserved
        - molecule: Adsorbate in original cell at original position
        - complete: Original structure (for provenance)

    Raises:
        ValueError: If formula is invalid, adsorbate not found, or multiple matches
    """
    from pymatgen.io.ase import AseAtomsAdaptor
    import networkx as nx

    # Validation: check formula is not empty
    formula_str = adsorbate_formula.value.strip()
    if not formula_str:
        raise ValueError("Adsorbate formula cannot be empty")

    # Parse adsorbate formula
    try:
        target_composition = Counter(parse_formula(formula_str))
    except Exception as e:
        raise ValueError(f"Invalid chemical formula '{formula_str}': {e}")

    # Get pymatgen structure
    pmg_structure = structure.get_pymatgen()

    # Validate structure has enough atoms
    total_target_atoms = sum(target_composition.values())
    if len(pmg_structure) < total_target_atoms:
        raise ValueError(
            f"Structure has only {len(pmg_structure)} atoms, "
            f"cannot contain adsorbate {formula_str} ({total_target_atoms} atoms)"
        )

    # Strategy: Find adsorbate by building subgraph of candidate atoms only
    # This works even if adsorbate is bonded to substrate

    # Step 1: Find all atoms matching elements in adsorbate formula
    target_elements = set(target_composition.keys())
    candidate_indices = []
    for i, site in enumerate(pmg_structure):
        if str(site.specie) in target_elements:
            candidate_indices.append(i)

    if len(candidate_indices) < total_target_atoms:
        raise ValueError(
            f"Not enough atoms to form adsorbate {formula_str}. "
            f"Need {total_target_atoms} atoms from {target_elements}, "
            f"but structure only has {len(candidate_indices)} matching atoms."
        )

    # Step 2: Build connectivity graph using CrystalNN
    sg = StructureGraph.from_local_env_strategy(pmg_structure, CrystalNN())

    # Step 3: Build subgraph containing only edges between candidate atoms
    candidate_graph = nx.Graph()
    candidate_graph.add_nodes_from(candidate_indices)

    for i in candidate_indices:
        # Get all bonded neighbors from CrystalNN graph
        neighbors = sg.get_connected_sites(i)
        for neighbor_info in neighbors:
            j = neighbor_info.index
            if j in candidate_indices and i < j:
                candidate_graph.add_edge(i, j)

    # Step 4: Find connected components
    components = list(nx.connected_components(candidate_graph))

    # Step 5: Find cluster matching adsorbate formula
    matching_clusters = []
    for component in components:
        node_indices = sorted(list(component))
        cluster_composition = Counter()
        for idx in node_indices:
            element = str(pmg_structure[idx].specie)
            cluster_composition[element] += 1
        if cluster_composition == target_composition:
            matching_clusters.append(node_indices)

    # FALLBACK: If CrystalNN didn't find the adsorbate, try distance-based connectivity
    # This handles cases like O-O peroxide bonds that CrystalNN doesn't detect
    if len(matching_clusters) == 0:
        # Rebuild graph using distance-based bonds
        candidate_graph_distance = nx.Graph()
        candidate_graph_distance.add_nodes_from(candidate_indices)

        for i_idx, i in enumerate(candidate_indices):
            pos_i = np.array(pmg_structure[i].coords)
            for j in candidate_indices[i_idx + 1:]:
                pos_j = np.array(pmg_structure[j].coords)
                dist = np.linalg.norm(pos_i - pos_j)
                if dist <= bond_cutoff:
                    candidate_graph_distance.add_edge(i, j)

        # Find components in distance-based graph
        components_distance = list(nx.connected_components(candidate_graph_distance))

        for component in components_distance:
            node_indices = sorted(list(component))
            cluster_composition = Counter()
            for idx in node_indices:
                element = str(pmg_structure[idx].specie)
                cluster_composition[element] += 1
            if cluster_composition == target_composition:
                matching_clusters.append(node_indices)

    # Validation: check we found exactly one match
    if len(matching_clusters) == 0:
        available_formulas = []
        for component in components:
            node_indices = list(component)
            comp = Counter()
            for idx in node_indices:
                element = str(pmg_structure[idx].specie)
                comp[element] += 1
            formula = "".join(f"{el}{count}" for el, count in sorted(comp.items()))
            available_formulas.append(formula)

        raise ValueError(
            f"Could not find adsorbate '{formula_str}' among candidate atom clusters. "
            f"Found clusters with formulas: {available_formulas}. "
            f"Tried CrystalNN and distance-based (cutoff={bond_cutoff} A) connectivity."
        )

    if len(matching_clusters) > 1:
        raise ValueError(
            f"Found {len(matching_clusters)} clusters matching '{formula_str}'. "
            f"Structure is ambiguous - please provide structure with unique adsorbate."
        )

    # Get the matched adsorbate indices
    adsorbate_indices = set(matching_clusters[0])

    # Validation: check substrate is not too small
    substrate_indices = [i for i in range(len(pmg_structure)) if i not in adsorbate_indices]
    if len(substrate_indices) < 3:
        raise ValueError(
            f"Substrate only has {len(substrate_indices)} atoms after removing adsorbate. "
            f"Structure may be invalid."
        )

    # Create three structures using ASE for manipulation
    adaptor = AseAtomsAdaptor()
    ase_structure = adaptor.get_atoms(pmg_structure)

    # Convert indices to sorted lists
    adsorbate_indices_list = sorted(list(adsorbate_indices), reverse=True)
    substrate_indices_list = sorted(substrate_indices, reverse=True)

    # 1. Substrate: remove adsorbate atoms, keep cell
    substrate_ase = ase_structure.copy()
    for idx in adsorbate_indices_list:
        del substrate_ase[idx]

    # 2. Molecule: keep only adsorbate atoms, keep cell and position
    molecule_ase = ase_structure.copy()
    for idx in substrate_indices_list:
        del molecule_ase[idx]

    # 3. Complete: return as-is for provenance
    complete_ase = ase_structure.copy()

    # Convert back to AiiDA StructureData
    return {
        'substrate': orm.StructureData(ase=substrate_ase),
        'molecule': orm.StructureData(ase=molecule_ase),
        'complete': orm.StructureData(ase=complete_ase),
    }


@task.calcfunction
def separate_adsorbate_structure(
    structure: orm.StructureData,
    adsorbate_formula: orm.Str,
) -> t.Annotated[dict, namespace(
    substrate=orm.StructureData,
    molecule=orm.StructureData,
    complete=orm.StructureData
)]:
    """
    Separate substrate+adsorbate structure into three components.

    Uses connectivity analysis (pymatgen StructureGraph) to identify bonded
    clusters and match the adsorbate by chemical formula.

    Args:
        structure: Complete structure (substrate + adsorbate)
        adsorbate_formula: Chemical formula of adsorbate (e.g., "OOH", "OH")

    Returns:
        Dictionary with three StructureData nodes:
        - substrate: Structure with adsorbate removed, original cell preserved
        - molecule: Adsorbate in original cell at original position
        - complete: Original structure (for provenance)

    Raises:
        ValueError: If formula is invalid, adsorbate not found, or multiple matches
    """
    return _separate_adsorbate_structure_impl(structure, adsorbate_formula)


def _calculate_adsorption_energy_impl(
    E_complete: orm.Float,
    E_substrate: orm.Float,
    E_molecule: orm.Float,
) -> orm.Float:
    """
    Internal implementation for calculating adsorption energy.

    Args:
        E_complete: Total energy of substrate+adsorbate system (eV)
        E_substrate: Total energy of bare substrate (eV)
        E_molecule: Total energy of isolated adsorbate molecule (eV)

    Returns:
        Adsorption energy in eV
    """
    E_ads = E_complete.value - E_substrate.value - E_molecule.value
    return orm.Float(E_ads)


@task.calcfunction
def calculate_adsorption_energy(
    E_complete: orm.Float,
    E_substrate: orm.Float,
    E_molecule: orm.Float,
) -> orm.Float:
    """
    Calculate adsorption energy from component energies.

    Uses the formula:
        E_ads = E_complete - E_substrate - E_molecule

    A negative adsorption energy indicates exothermic (favorable) adsorption.
    A positive adsorption energy indicates endothermic (unfavorable) adsorption.

    Args:
        E_complete: Total energy of substrate+adsorbate system (eV)
        E_substrate: Total energy of bare substrate (eV)
        E_molecule: Total energy of isolated adsorbate molecule (eV)

    Returns:
        Adsorption energy in eV
    """
    return _calculate_adsorption_energy_impl(E_complete, E_substrate, E_molecule)


@task.graph
def compute_adsorption_energies_scatter(
    structures: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    adsorbate_formulas: t.Annotated[dict[str, str], dict],
    code: orm.Code,
    potential_family: str,
    potential_mapping: t.Mapping[str, str],

    # Relaxation control
    relax_before_adsorption: bool = False,
    relax_parameters: t.Mapping[str, t.Any] | None = None,

    # SCF parameters
    scf_parameters: t.Mapping[str, t.Any] | None = None,

    # Common settings
    options: t.Mapping[str, t.Any] | None = None,
    kpoints_spacing: float | None = None,
    clean_workdir: bool = True,

    # Builder inputs (NEW - complete control over VASP parameters)
    relax_builder_inputs: dict | None = None,
    scf_builder_inputs: dict | None = None,
    structure_specific_relax_builder_inputs: dict | None = None,
    structure_specific_scf_builder_inputs: dict | None = None,
    structure_component_specific_scf_builder_inputs: dict | None = None,

    # Atom fixing control (for relaxation)
    fix_atoms: bool = False,
    fix_type: str | None = None,
    fix_thickness: float = 0.0,
    fix_elements: list[str] | None = None,

    # Concurrency control
    max_number_jobs: int = None,
) -> t.Annotated[dict, namespace(
    relaxed_complete_structures=dynamic(orm.StructureData),
    separated_structures=dynamic(dict),
    substrate_energies=dynamic(orm.Float),
    molecule_energies=dynamic(orm.Float),
    complete_energies=dynamic(orm.Float),
    adsorption_energies=dynamic(orm.Float),
)]:
    """
    Scatter-gather workflow for calculating adsorption energies using vasp.v2.vasp.

    This workflow supports two modes:
    1. Direct SCF: separate → SCF (3N jobs)
    2. Relax + SCF: relax → separate → SCF (N+3N jobs)

    Workflow phases:
    - Phase 1 (optional): Relax complete structures using vasp.v2.vasp (NSW > 0)
    - Phase 2: Separate relaxed (or original) structures into substrate/molecule/complete
    - Phase 3: SCF calculations using vasp.v2.vasp (NSW=0) for all three components

    Args:
        structures: Dynamic namespace of complete structures (substrate + adsorbate)
        adsorbate_formulas: Dictionary mapping structure keys to adsorbate formulas
                           Example: {'system1': 'OOH', 'system2': 'OH'}
        code: AiiDA code for VASP
        potential_family: Pseudopotential family name
        potential_mapping: Element to potential mapping

        relax_before_adsorption: If True, relax complete structures before separation
        relax_parameters: INCAR parameters for relaxation (must include NSW > 0, IBRION, ISIF)
        scf_parameters: INCAR parameters for SCF (NSW will be overridden to 0)

        options: Scheduler options dict
        kpoints_spacing: K-points spacing in Angstrom^-1
        clean_workdir: Whether to clean remote working directories

        relax_builder_inputs: Complete VASP builder dict for relaxation (NEW)
                             Takes priority over relax_parameters if provided
                             Format: {'parameters': {'incar': {...}}, 'options': {...}, 'kpoints_spacing': 0.2, ...}
        scf_builder_inputs: Complete VASP builder dict for SCF calculations (NEW)
                           Takes priority over scf_parameters if provided
                           Format: {'parameters': {'incar': {...}}, 'options': {...}, 'kpoints_spacing': 0.2, ...}
        structure_specific_relax_builder_inputs: Per-structure overrides for relaxation (NEW)
                                                 Dict mapping structure indices to builder dicts
                                                 Uses deep merge - only specified parameters are overridden
                                                 Format: {0: {'parameters': {'incar': {'ALGO': 'Normal'}}}, ...}
        structure_specific_scf_builder_inputs: Per-structure overrides for SCF calculations (NEW)
                                               Dict mapping structure indices to builder dicts
                                               Applies to all three SCF calculations (substrate, molecule, complete)
                                               Uses deep merge - only specified parameters are overridden
                                               Format: {1: {'kpoints_spacing': 0.15}, ...}
        structure_component_specific_scf_builder_inputs: Per-structure, per-component overrides (NEW)
                                                         Dict mapping structure indices to component-specific dicts
                                                         Allows modifying only substrate/molecule/complete for specific structures
                                                         Uses deep merge - only specified parameters are overridden
                                                         Format: {2: {'substrate': {'parameters': {'incar': {'ALGO': 'Fast'}}}}, ...}
                                                         Components: 'substrate', 'molecule', 'complete'

        fix_atoms: If True, apply atom fixing constraints during relaxation
        fix_type: Type of fixing ('bottom', 'top', 'center')
        fix_thickness: Thickness in Angstrom for fixing region
        fix_elements: Optional list of elements to fix (None = all elements)

        max_number_jobs: Maximum number of concurrent VASP calculations (None = unlimited)

    Returns:
        Dictionary with namespaces:
        - relaxed_complete_structures: Relaxed structures (empty if relax=False)
        - separated_structures: Dict of separated systems for each input
        - substrate_energies: Energies of bare substrates
        - molecule_energies: Energies of isolated molecules
        - complete_energies: Energies of complete systems
        - adsorption_energies: Final adsorption energies (E_ads = E_complete - E_substrate - E_molecule)

    Example:
        >>> relax_params = {'NSW': 200, 'IBRION': 2, 'ISIF': 2, 'EDIFFG': -0.05, 'ENCUT': 520, ...}
        >>> scf_params = {'ENCUT': 520, 'EDIFF': 1e-6, ...}  # NSW=0 added automatically
        >>> options = {'resources': {'num_machines': 1, 'num_cores_per_machine': 40}, 'queue_name': 'par40'}
        >>> results = compute_adsorption_energies_scatter(
        ...     structures={'oh_lamno3': structure},
        ...     adsorbate_formulas={'oh_lamno3': 'OH'},
        ...     code=code,
        ...     potential_family='PBE',
        ...     potential_mapping={'La': 'La', 'Ni': 'Ni', 'O': 'O', 'H': 'H'},
        ...     relax_before_adsorption=True,
        ...     relax_parameters=relax_params,
        ...     scf_parameters=scf_params,
        ...     options=options,
        ...     kpoints_spacing=0.6,
        ... )
    """
    from aiida_workgraph import get_current_graph

    # Set max_number_jobs on this workgraph to control concurrency
    if max_number_jobs is not None:
        wg = get_current_graph()
        max_jobs_value = extract_max_jobs_value(max_number_jobs)
        wg.max_number_jobs = max_jobs_value

    # Validate inputs
    if structures.keys() != adsorbate_formulas.keys():
        raise ValueError(
            f"Mismatch between structures and adsorbate_formulas keys. "
            f"Structures: {list(structures.keys())}, "
            f"Formulas: {list(adsorbate_formulas.keys())}"
        )

    # Get vasp.v2.vasp workflow
    try:
        vasp_wc = WorkflowFactory('vasp.v2.vasp')
    except Exception as e:
        raise RuntimeError(
            f"Failed to load vasp.v2.vasp workflow plugin. "
            f"Ensure aiida-vasp is installed and registered. Error: {e}"
        ) from e

    vasp_task_cls = task(vasp_wc)

    # ===== PHASE 1: OPTIONAL RELAXATION =====
    # Relax complete structures before separation if requested
    relaxed_structures: dict[str, orm.StructureData] = {}

    if relax_before_adsorption:
        # Validate that we have either builder_inputs or parameters
        if relax_builder_inputs is None and relax_parameters is None:
            raise ValueError(
                "When relax_before_adsorption=True, either relax_builder_inputs or "
                "relax_parameters must be provided"
            )

        # Create index mapping: key -> index
        key_to_index = {key: idx for idx, key in enumerate(structures.keys())}

        # Build default relax builder inputs
        if relax_builder_inputs is not None:
            # New style: use builder_inputs directly
            default_relax_builder = copy.deepcopy(relax_builder_inputs)
        else:
            # Old style: construct from parameters
            default_relax_builder = {
                'parameters': {'incar': dict(relax_parameters)},
                'options': dict(options) if options is not None else {},
                'potential_family': potential_family,
                'potential_mapping': dict(potential_mapping),
                'clean_workdir': clean_workdir,
                'settings': orm.Dict(dict=get_settings()),
            }
            if kpoints_spacing is not None:
                default_relax_builder['kpoints_spacing'] = kpoints_spacing

        # Convert structure_specific keys from int to str for lookup
        structure_specific_relax = {}
        if structure_specific_relax_builder_inputs is not None:
            for idx_key, override in structure_specific_relax_builder_inputs.items():
                # Accept both int and str keys
                str_key = str(idx_key)
                structure_specific_relax[str_key] = override

        # Process each structure
        for key, structure in structures.items():
            # Get structure index
            idx = key_to_index[key]
            idx_str = str(idx)

            # Deep merge structure-specific overrides if they exist
            if idx_str in structure_specific_relax:
                relax_inputs = deep_merge_dicts(
                    default_relax_builder,
                    structure_specific_relax[idx_str]
                )
            else:
                relax_inputs = copy.deepcopy(default_relax_builder)

            # Always set structure and code (override any user-provided values)
            relax_inputs['structure'] = structure
            relax_inputs['code'] = code

            # Apply atom fixing if requested (modifies relax_inputs in-place)
            if fix_atoms and fix_type is not None and fix_thickness > 0.0:
                fixed_atoms_list = get_fixed_atoms_list(
                    structure=relax_inputs['structure'],
                    fix_type=fix_type,
                    fix_thickness=fix_thickness,
                    fix_elements=fix_elements,
                )

                if fixed_atoms_list:
                    # Get parameters location
                    if 'parameters' in relax_inputs and 'incar' in relax_inputs['parameters']:
                        incar_params = relax_inputs['parameters']['incar']
                    else:
                        relax_inputs['parameters'] = {'incar': {}}
                        incar_params = relax_inputs['parameters']['incar']

                    # Apply fixing
                    updated_params, constrained_structure = add_fixed_atoms_to_vasp_parameters(
                        base_parameters=incar_params,
                        structure=relax_inputs['structure'],
                        fixed_atoms_list=fixed_atoms_list,
                    )

                    # Update inputs
                    relax_inputs['parameters']['incar'] = updated_params
                    relax_inputs['structure'] = constrained_structure

            # Run relaxation using vasp.v2.vasp
            relax_task = vasp_task_cls(**relax_inputs)

            # Extract relaxed structure
            relaxed_structures[key] = relax_task.structure

        # Use relaxed structures for separation
        structures_to_separate = relaxed_structures
    else:
        # No relaxation: use original structures
        structures_to_separate = structures

    # Output dictionaries
    separated_structures_dict: dict[str, dict] = {}
    substrate_energies: dict[str, orm.Float] = {}
    molecule_energies: dict[str, orm.Float] = {}
    complete_energies: dict[str, orm.Float] = {}
    adsorption_energies: dict[str, orm.Float] = {}

    # ===== PHASE 2: STRUCTURE SEPARATION =====
    # Separate relaxed (if Phase 1 ran) or original structures
    for key, structure in structures_to_separate.items():
        adsorbate_str = adsorbate_formulas[key]

        separated = separate_adsorbate_structure(
            structure=structure,
            adsorbate_formula=orm.Str(adsorbate_str)
        )

        separated_structures_dict[key] = {
            'substrate': separated.substrate,
            'molecule': separated.molecule,
            'complete': separated.complete,
        }

    # ===== PHASE 3: SCF CALCULATIONS =====
    # Single-point calculations for substrate, molecule, and complete systems

    # Build default SCF builder inputs
    if scf_builder_inputs is not None:
        # New style: use builder_inputs directly
        default_scf_builder = copy.deepcopy(scf_builder_inputs)
    elif scf_parameters is not None:
        # Old style: construct from scf_parameters
        default_scf_builder = {
            'parameters': {'incar': dict(scf_parameters)},
            'options': dict(options) if options is not None else {},
            'potential_family': potential_family,
            'potential_mapping': dict(potential_mapping),
            'clean_workdir': clean_workdir,
            'settings': orm.Dict(dict=get_settings()),
        }
        if kpoints_spacing is not None:
            default_scf_builder['kpoints_spacing'] = kpoints_spacing
    elif relax_parameters is not None:
        # Fallback: use relax_parameters as base for SCF
        default_scf_builder = {
            'parameters': {'incar': dict(relax_parameters)},
            'options': dict(options) if options is not None else {},
            'potential_family': potential_family,
            'potential_mapping': dict(potential_mapping),
            'clean_workdir': clean_workdir,
            'settings': orm.Dict(dict=get_settings()),
        }
        if kpoints_spacing is not None:
            default_scf_builder['kpoints_spacing'] = kpoints_spacing
    else:
        # No parameters provided at all - use minimal defaults
        default_scf_builder = {
            'parameters': {'incar': {}},
            'options': dict(options) if options is not None else {},
            'potential_family': potential_family,
            'potential_mapping': dict(potential_mapping),
            'clean_workdir': clean_workdir,
            'settings': orm.Dict(dict=get_settings()),
        }
        if kpoints_spacing is not None:
            default_scf_builder['kpoints_spacing'] = kpoints_spacing

    # Convert structure_specific keys from int to str for lookup
    structure_specific_scf = {}
    if structure_specific_scf_builder_inputs is not None:
        for idx_key, override in structure_specific_scf_builder_inputs.items():
            # Accept both int and str keys
            str_key = str(idx_key)
            structure_specific_scf[str_key] = override

    # Convert structure_component_specific keys from int to str for lookup
    structure_component_specific_scf = {}
    if structure_component_specific_scf_builder_inputs is not None:
        for idx_key, component_overrides in structure_component_specific_scf_builder_inputs.items():
            # Accept both int and str keys
            str_key = str(idx_key)
            structure_component_specific_scf[str_key] = component_overrides

    # Recreate key_to_index mapping for SCF phase
    key_to_index = {key: idx for idx, key in enumerate(structures.keys())}

    # Process each structure's SCF calculations
    for key, separated in separated_structures_dict.items():
        # Get structure index
        idx = key_to_index[key]
        idx_str = str(idx)

        # Deep merge structure-specific overrides if they exist
        # These overrides apply to all three SCF calculations for this structure
        if idx_str in structure_specific_scf:
            base_scf_inputs = deep_merge_dicts(
                default_scf_builder,
                structure_specific_scf[idx_str]
            )
        else:
            base_scf_inputs = copy.deepcopy(default_scf_builder)

        # Run SCF calculations for all three components
        scf_components = [
            ('substrate', separated['substrate'], substrate_energies),
            ('molecule', separated['molecule'], molecule_energies),
            ('complete', separated['complete'], complete_energies),
        ]

        for component_name, struct, energy_dict in scf_components:
            # Make a copy of base_scf_inputs for this component
            scf_inputs = copy.deepcopy(base_scf_inputs)

            # Apply component-specific overrides if they exist
            # Priority: base_scf_inputs < structure_specific < component_specific
            if idx_str in structure_component_specific_scf:
                component_overrides = structure_component_specific_scf[idx_str]
                if component_name in component_overrides:
                    scf_inputs = deep_merge_dicts(
                        scf_inputs,
                        component_overrides[component_name]
                    )

            # Always set structure and code
            scf_inputs['structure'] = struct
            scf_inputs['code'] = code

            # Force single-point calculation (NSW=0, IBRION=-1)
            if 'parameters' not in scf_inputs:
                scf_inputs['parameters'] = {'incar': {}}
            if 'incar' not in scf_inputs['parameters']:
                scf_inputs['parameters']['incar'] = {}

            scf_inputs['parameters']['incar']['NSW'] = 0
            scf_inputs['parameters']['incar']['IBRION'] = -1

            # Run SCF using vasp.v2.vasp
            scf_calc = vasp_task_cls(**scf_inputs)
            energy_dict[key] = extract_total_energy(energies=scf_calc.misc).result

    # ===== PHASE 4: ADSORPTION ENERGY CALCULATION =====
    # Calculate E_ads = E_complete - E_substrate - E_molecule
    for key in structures.keys():
        E_ads = calculate_adsorption_energy(
            E_complete=complete_energies[key],
            E_substrate=substrate_energies[key],
            E_molecule=molecule_energies[key],
        ).result

        adsorption_energies[key] = E_ads

    # Return all results
    return {
        'relaxed_complete_structures': relaxed_structures,
        'separated_structures': separated_structures_dict,
        'substrate_energies': substrate_energies,
        'molecule_energies': molecule_energies,
        'complete_energies': complete_energies,
        'adsorption_energies': adsorption_energies,
    }


@task.workfunction
def extract_from_completed_workgraph(
    workgraph_pk: orm.Int,
) -> t.Annotated[dict, namespace(
    structures=dynamic(orm.StructureData),
    energies=dynamic(orm.Float),
    retrieved=dynamic(orm.FolderData),
)]:
    """
    Extract structures, energies, and retrieved files from a completed WorkGraph.

    This function traverses the WorkGraph outputs namespace to extract the complete
    system data needed for adsorption energy calculations when the complete system
    SCF has already been run.

    For SCF calculations (NSW=0), structures may be in inputs rather than outputs
    since they don't change during the calculation. This function checks both
    locations.

    Args:
        workgraph_pk: PK of the completed WorkGraph node

    Returns:
        Dictionary with namespaces:
        - structures: Dict of StructureData keyed by structure label
        - energies: Dict of Float energies keyed by structure label
        - retrieved: Dict of FolderData keyed by structure label

    Example:
        >>> data = extract_from_completed_workgraph(orm.Int(32945))
        >>> print(data['energies'])  # {'ooh_on_O_000': <Float>, ...}
    """
    wg = orm.load_node(workgraph_pk.value)

    structures = {}
    energies = {}
    retrieved = {}

    # Handle namespaced outputs (multiple structures)
    if hasattr(wg.outputs, 'energies'):
        # Multiple structures with namespaced outputs
        for key in wg.outputs.energies:
            energies[key] = wg.outputs.energies[key]

        # Try to get structures from outputs first
        if hasattr(wg.outputs, 'structures'):
            for key in wg.outputs.structures:
                structures[key] = wg.outputs.structures[key]
        else:
            # Fall back to graph_inputs.structures (common for SCF calculations)
            if hasattr(wg.inputs, 'graph_inputs') and hasattr(wg.inputs.graph_inputs, 'structures'):
                for key in wg.inputs.graph_inputs.structures:
                    structures[key] = wg.inputs.graph_inputs.structures[key]

        if hasattr(wg.outputs, 'retrieved'):
            for key in wg.outputs.retrieved:
                retrieved[key] = wg.outputs.retrieved[key]

    elif hasattr(wg.outputs, 'energy'):
        # Single structure - wrap in dict with default key
        energies['structure'] = wg.outputs.energy

        # Try outputs first, then inputs
        if hasattr(wg.outputs, 'structure'):
            structures['structure'] = wg.outputs.structure
        elif hasattr(wg.inputs, 'graph_inputs') and hasattr(wg.inputs.graph_inputs, 'structure'):
            structures['structure'] = wg.inputs.graph_inputs.structure

        if hasattr(wg.outputs, 'retrieved'):
            retrieved['structure'] = wg.outputs.retrieved

    else:
        raise ValueError(
            f"WorkGraph {workgraph_pk.value} does not have expected output structure. "
            f"Expected 'energy' or 'energies' namespace in outputs."
        )

    # Validate that we have structures for all energies
    if not structures:
        raise ValueError(
            f"Could not find structures in WorkGraph {workgraph_pk.value}. "
            f"Checked outputs.structures and inputs.graph_inputs.structures."
        )

    return {
        'structures': structures,
        'energies': energies,
        'retrieved': retrieved,
    }


@task.graph
def compute_adsorption_energies_from_completed(
    # Input Option 1: From completed WorkGraph PK
    completed_workgraph_pk: int = None,

    # Input Option 2: Direct inputs (takes priority if provided)
    complete_structures: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)] = None,
    complete_energies: t.Annotated[dict[str, orm.Float], dynamic(orm.Float)] = None,
    complete_retrieved: t.Annotated[dict[str, orm.FolderData], dynamic(orm.FolderData)] = None,

    # Adsorbate specification
    adsorbate_formulas: t.Annotated[dict[str, str], dict] = None,

    # VASP code configuration
    code: orm.Code = None,

    # Builder inputs for SCF calculations
    scf_builder_inputs: dict = None,
    structure_specific_scf_builder_inputs: dict = None,
    structure_component_specific_scf_builder_inputs: dict = None,

    # Legacy parameter support
    potential_family: str = None,
    potential_mapping: t.Mapping[str, str] = None,
    scf_parameters: t.Mapping[str, t.Any] = None,
    options: t.Mapping[str, t.Any] = None,
    kpoints_spacing: float = None,
    clean_workdir: bool = True,

    # Bader analysis support
    retrieve_bader_files: bool = True,

    # Concurrency control
    max_number_jobs: int = None,
) -> t.Annotated[dict, namespace(
    separated_structures=dynamic(dict),
    substrate_energies=dynamic(orm.Float),
    molecule_energies=dynamic(orm.Float),
    complete_energies=dynamic(orm.Float),
    adsorption_energies=dynamic(orm.Float),
    substrate_retrieved=dynamic(orm.FolderData),
    molecule_retrieved=dynamic(orm.FolderData),
    complete_retrieved=dynamic(orm.FolderData),
)]:
    """
    Calculate adsorption energies when the complete system SCF is already done.

    This workflow is designed for cases where you have already run SCF calculations
    on complete (substrate+adsorbate) structures and want to calculate adsorption
    energies by running SCF only on the separated components (substrate, molecule).

    Workflow phases:
    1. Extract: Get structures, energies, and retrieved files from completed WorkGraph
    2. Separate: Split complete structures into substrate and molecule components
    3. SCF: Run parallel SCF calculations on substrate and molecule only (skip complete)
    4. Calculate: Compute E_ads = E_complete - E_substrate - E_molecule

    Args:
        completed_workgraph_pk: PK of completed WorkGraph with complete system SCF results
        complete_structures: Direct input of complete structures (alternative to PK)
        complete_energies: Direct input of complete system energies (alternative to PK)
        complete_retrieved: Direct input of retrieved files from complete SCF (alternative to PK)

        adsorbate_formulas: Dictionary mapping structure keys to adsorbate formulas
                           Example: {'ooh_on_O_000': 'OOH', 'ooh_on_Ni_000': 'OOH'}

        code: AiiDA code for VASP
        scf_builder_inputs: Complete VASP builder dict for SCF calculations
                           Format: {'parameters': {'incar': {...}}, 'options': {...}, ...}
        structure_specific_scf_builder_inputs: Per-structure overrides
        structure_component_specific_scf_builder_inputs: Per-structure, per-component overrides

        potential_family: Pseudopotential family name (legacy support)
        potential_mapping: Element to potential mapping (legacy support)
        scf_parameters: INCAR parameters for SCF (legacy support)
        options: Scheduler options dict (legacy support)
        kpoints_spacing: K-points spacing in Angstrom^-1 (legacy support)
        clean_workdir: Whether to clean remote working directories

        retrieve_bader_files: If True, add CHGCAR/AECCAR files to retrieve list
        max_number_jobs: Maximum number of concurrent VASP calculations

    Returns:
        Dictionary with namespaces:
        - separated_structures: Dict of {key: {substrate, molecule, complete}}
        - substrate_energies: Energies of bare substrates
        - molecule_energies: Energies of isolated molecules
        - complete_energies: Energies from input (passed through)
        - adsorption_energies: Final E_ads values
        - substrate_retrieved: FolderData from substrate SCF
        - molecule_retrieved: FolderData from molecule SCF
        - complete_retrieved: FolderData from input (passed through)

    Example:
        >>> wg = compute_adsorption_energies_from_completed(
        ...     completed_workgraph_pk=32945,
        ...     adsorbate_formulas={
        ...         'ooh_on_O_000_O32_top': 'OOH',
        ...         'ooh_on_Ni_000_Ni13_top': 'OOH',
        ...     },
        ...     code=orm.load_code('VASP-6.5.1@cluster'),
        ...     scf_builder_inputs={
        ...         'parameters': {'incar': {'ENCUT': 500, 'EDIFF': 1e-4}},
        ...         'options': {'resources': {'num_machines': 1}},
        ...         'kpoints_spacing': 0.06,
        ...         'potential_family': 'PBE',
        ...         'potential_mapping': {'La': 'La', 'Ni': 'Ni_pv', 'O': 'O', 'H': 'H'},
        ...     },
        ...     max_number_jobs=4,
        ... )
        >>> wg.name = 'adsorption_energy_from_completed'
        >>> wg.submit()
    """
    from aiida_workgraph import get_current_graph

    # Set max_number_jobs on this workgraph to control concurrency
    if max_number_jobs is not None:
        wg = get_current_graph()
        max_jobs_value = extract_max_jobs_value(max_number_jobs)
        wg.max_number_jobs = max_jobs_value

    # ===== PHASE 1: INPUT EXTRACTION =====
    # Get complete system data from either PK or direct inputs
    if complete_structures is not None and complete_energies is not None:
        # Direct inputs provided
        structures = complete_structures
        energies = complete_energies
        retrieved = complete_retrieved if complete_retrieved is not None else {}
    elif completed_workgraph_pk is not None:
        # Extract from completed WorkGraph
        extracted = extract_from_completed_workgraph(workgraph_pk=orm.Int(completed_workgraph_pk))
        structures = extracted.structures
        energies = extracted.energies
        retrieved = extracted.retrieved if hasattr(extracted, 'retrieved') else {}
    else:
        raise ValueError(
            "Must provide either completed_workgraph_pk or "
            "(complete_structures + complete_energies)"
        )

    # Validate adsorbate_formulas
    if adsorbate_formulas is None:
        raise ValueError("adsorbate_formulas is required")

    # Convert structures namespace to dict if needed
    structures_dict = {}
    for key in structures:
        structures_dict[key] = structures[key]

    energies_dict = {}
    for key in energies:
        energies_dict[key] = energies[key]

    retrieved_dict = {}
    for key in retrieved:
        retrieved_dict[key] = retrieved[key]

    # Validate keys match
    missing_formulas = set(structures_dict.keys()) - set(adsorbate_formulas.keys())
    if missing_formulas:
        raise ValueError(
            f"Missing adsorbate formulas for structures: {missing_formulas}. "
            f"Provide formulas for all structures."
        )

    # ===== PHASE 2: STRUCTURE SEPARATION =====
    separated_structures_dict: dict[str, dict] = {}

    for key, structure in structures_dict.items():
        adsorbate_str = adsorbate_formulas[key]

        separated = separate_adsorbate_structure(
            structure=structure,
            adsorbate_formula=orm.Str(adsorbate_str)
        )

        separated_structures_dict[key] = {
            'substrate': separated.substrate,
            'molecule': separated.molecule,
            'complete': separated.complete,
        }

    # ===== PHASE 3: SCF CALCULATIONS (SUBSTRATE + MOLECULE ONLY) =====
    # Get vasp.v2.vasp workflow
    try:
        vasp_wc = WorkflowFactory('vasp.v2.vasp')
    except Exception as e:
        raise RuntimeError(
            f"Failed to load vasp.v2.vasp workflow plugin. "
            f"Ensure aiida-vasp is installed. Error: {e}"
        ) from e

    vasp_task_cls = task(vasp_wc)

    # Build default SCF builder inputs
    if scf_builder_inputs is not None:
        default_scf_builder = copy.deepcopy(scf_builder_inputs)
    elif scf_parameters is not None:
        default_scf_builder = {
            'parameters': {'incar': dict(scf_parameters)},
            'options': dict(options) if options is not None else {},
            'potential_family': potential_family,
            'potential_mapping': dict(potential_mapping) if potential_mapping is not None else {},
            'clean_workdir': clean_workdir,
            'settings': orm.Dict(dict=get_settings()),
        }
        if kpoints_spacing is not None:
            default_scf_builder['kpoints_spacing'] = kpoints_spacing
    else:
        raise ValueError(
            "Must provide either scf_builder_inputs or scf_parameters"
        )

    # Add Bader files to retrieve list if requested
    if retrieve_bader_files:
        bader_files = ['CHGCAR', 'AECCAR0', 'AECCAR1', 'AECCAR2']
        if 'settings' not in default_scf_builder:
            default_scf_builder['settings'] = {}
        if isinstance(default_scf_builder['settings'], orm.Dict):
            settings_dict = default_scf_builder['settings'].get_dict()
        else:
            settings_dict = dict(default_scf_builder['settings'])

        if 'ADDITIONAL_RETRIEVE_LIST' not in settings_dict:
            settings_dict['ADDITIONAL_RETRIEVE_LIST'] = []

        for f in bader_files:
            if f not in settings_dict['ADDITIONAL_RETRIEVE_LIST']:
                settings_dict['ADDITIONAL_RETRIEVE_LIST'].append(f)

        default_scf_builder['settings'] = settings_dict

    # Convert structure_specific keys
    structure_specific_scf = {}
    if structure_specific_scf_builder_inputs is not None:
        for idx_key, override in structure_specific_scf_builder_inputs.items():
            str_key = str(idx_key)
            structure_specific_scf[str_key] = override

    structure_component_specific_scf = {}
    if structure_component_specific_scf_builder_inputs is not None:
        for idx_key, component_overrides in structure_component_specific_scf_builder_inputs.items():
            str_key = str(idx_key)
            structure_component_specific_scf[str_key] = component_overrides

    # Create key_to_index mapping
    key_to_index = {key: idx for idx, key in enumerate(structures_dict.keys())}

    # Output dictionaries
    substrate_energies: dict[str, orm.Float] = {}
    molecule_energies: dict[str, orm.Float] = {}
    substrate_retrieved: dict[str, orm.FolderData] = {}
    molecule_retrieved: dict[str, orm.FolderData] = {}
    adsorption_energies: dict[str, orm.Float] = {}

    # Run SCF for substrate and molecule only (NOT complete - already done)
    for key, separated in separated_structures_dict.items():
        idx = key_to_index[key]
        idx_str = str(idx)

        # Deep merge structure-specific overrides
        if idx_str in structure_specific_scf:
            base_scf_inputs = deep_merge_dicts(
                default_scf_builder,
                structure_specific_scf[idx_str]
            )
        else:
            base_scf_inputs = copy.deepcopy(default_scf_builder)

        # Only run SCF for substrate and molecule (skip complete)
        scf_components = [
            ('substrate', separated['substrate'], substrate_energies, substrate_retrieved),
            ('molecule', separated['molecule'], molecule_energies, molecule_retrieved),
        ]

        for component_name, struct, energy_dict, retrieved_dict_out in scf_components:
            scf_inputs = copy.deepcopy(base_scf_inputs)

            # Apply component-specific overrides
            if idx_str in structure_component_specific_scf:
                component_overrides = structure_component_specific_scf[idx_str]
                if component_name in component_overrides:
                    scf_inputs = deep_merge_dicts(
                        scf_inputs,
                        component_overrides[component_name]
                    )

            # Set structure and code
            scf_inputs['structure'] = struct
            scf_inputs['code'] = code

            # Force single-point calculation
            if 'parameters' not in scf_inputs:
                scf_inputs['parameters'] = {'incar': {}}
            if 'incar' not in scf_inputs['parameters']:
                scf_inputs['parameters']['incar'] = {}

            scf_inputs['parameters']['incar']['NSW'] = 0
            scf_inputs['parameters']['incar']['IBRION'] = -1

            # Run SCF
            scf_calc = vasp_task_cls(**scf_inputs)
            energy_dict[key] = extract_total_energy(energies=scf_calc.misc).result
            retrieved_dict_out[key] = scf_calc.retrieved

    # ===== PHASE 4: ADSORPTION ENERGY CALCULATION =====
    for key in structures_dict.keys():
        E_ads = calculate_adsorption_energy(
            E_complete=energies_dict[key],
            E_substrate=substrate_energies[key],
            E_molecule=molecule_energies[key],
        ).result

        adsorption_energies[key] = E_ads

    # Return all results
    return {
        'separated_structures': separated_structures_dict,
        'substrate_energies': substrate_energies,
        'molecule_energies': molecule_energies,
        'complete_energies': energies_dict,
        'adsorption_energies': adsorption_energies,
        'substrate_retrieved': substrate_retrieved,
        'molecule_retrieved': molecule_retrieved,
        'complete_retrieved': retrieved_dict,
    }


def _extract_from_completed_workgraph_impl(workgraph_pk: int) -> dict:
    """
    Extract structures, energies, and retrieved files from a completed WorkGraph.

    This is a non-task version for use at graph build time.

    Args:
        workgraph_pk: PK of the completed WorkGraph node

    Returns:
        Dictionary with keys: 'structures', 'energies', 'retrieved'
    """
    wg = orm.load_node(workgraph_pk)

    structures = {}
    energies = {}
    retrieved = {}

    # Handle namespaced outputs (multiple structures)
    if hasattr(wg.outputs, 'energies'):
        # Multiple structures with namespaced outputs
        for key in wg.outputs.energies:
            energies[key] = wg.outputs.energies[key]

        # Try to get structures from outputs first
        if hasattr(wg.outputs, 'structures'):
            for key in wg.outputs.structures:
                structures[key] = wg.outputs.structures[key]
        else:
            # Fall back to graph_inputs.structures (common for SCF calculations)
            if hasattr(wg.inputs, 'graph_inputs') and hasattr(wg.inputs.graph_inputs, 'structures'):
                for key in wg.inputs.graph_inputs.structures:
                    structures[key] = wg.inputs.graph_inputs.structures[key]

        if hasattr(wg.outputs, 'retrieved'):
            for key in wg.outputs.retrieved:
                retrieved[key] = wg.outputs.retrieved[key]

    elif hasattr(wg.outputs, 'energy'):
        # Single structure - wrap in dict with default key
        energies['structure'] = wg.outputs.energy

        # Try outputs first, then inputs
        if hasattr(wg.outputs, 'structure'):
            structures['structure'] = wg.outputs.structure
        elif hasattr(wg.inputs, 'graph_inputs') and hasattr(wg.inputs.graph_inputs, 'structure'):
            structures['structure'] = wg.inputs.graph_inputs.structure

        if hasattr(wg.outputs, 'retrieved'):
            retrieved['structure'] = wg.outputs.retrieved

    else:
        raise ValueError(
            f"WorkGraph {workgraph_pk} does not have expected output structure. "
            f"Expected 'energy' or 'energies' namespace in outputs."
        )

    # Validate that we have structures for all energies
    if not structures:
        raise ValueError(
            f"Could not find structures in WorkGraph {workgraph_pk}. "
            f"Checked outputs.structures and inputs.graph_inputs.structures."
        )

    return {
        'structures': structures,
        'energies': energies,
        'retrieved': retrieved,
    }


def build_adsorption_energy_from_completed_workgraph(
    completed_workgraph_pk: int = None,
    complete_structures: dict[str, orm.StructureData] = None,
    complete_energies: dict[str, orm.Float] = None,
    complete_retrieved: dict[str, orm.FolderData] = None,
    adsorbate_formulas: dict[str, str] = None,
    code_label: str = None,
    scf_builder_inputs: dict = None,
    retrieve_bader_files: bool = True,
    max_number_jobs: int = None,
    name: str = 'adsorption_energy_from_completed',
):
    """
    Build a WorkGraph for calculating adsorption energies from completed SCF calculations.

    This function builds a flat WorkGraph that:
    1. Separates complete structures into substrate and molecule components
    2. Runs SCF calculations on substrate and molecule only (complete already done)
    3. Calculates adsorption energies: E_ads = E_complete - E_substrate - E_molecule

    Args:
        completed_workgraph_pk: PK of completed WorkGraph with complete system SCF results
        complete_structures: Direct input of complete structures (alternative to PK)
        complete_energies: Direct input of complete system energies (alternative to PK)
        complete_retrieved: Direct input of retrieved files (alternative to PK)
        adsorbate_formulas: Dictionary mapping structure keys to adsorbate formulas
                           Example: {'ooh_on_O_000': 'OOH', 'ooh_on_Ni_000': 'OOH'}
        code_label: VASP code label (e.g., 'VASP-6.5.1@cluster')
        scf_builder_inputs: VASP builder parameters for SCF calculations
        retrieve_bader_files: If True, add CHGCAR/AECCAR files to retrieve list
        max_number_jobs: Maximum concurrent VASP calculations
        name: WorkGraph name

    Returns:
        WorkGraph ready to submit

    Example:
        >>> wg = build_adsorption_energy_from_completed_workgraph(
        ...     completed_workgraph_pk=32945,
        ...     adsorbate_formulas={'ooh_on_O_000': 'OOH'},
        ...     code_label='VASP-6.5.1@cluster',
        ...     scf_builder_inputs={...},
        ... )
        >>> wg.submit()
    """
    from aiida_workgraph import WorkGraph

    # Load code
    code = orm.load_code(code_label)

    # If completed_workgraph_pk is provided, extract data at build time
    if completed_workgraph_pk is not None and complete_structures is None:
        extracted = _extract_from_completed_workgraph_impl(completed_workgraph_pk)
        complete_structures = extracted['structures']
        complete_energies = extracted['energies']
        complete_retrieved = extracted['retrieved']

    # Validate inputs
    if complete_structures is None or complete_energies is None:
        raise ValueError("Must provide either completed_workgraph_pk or (complete_structures + complete_energies)")

    if adsorbate_formulas is None:
        raise ValueError("adsorbate_formulas is required")

    if scf_builder_inputs is None:
        raise ValueError("scf_builder_inputs is required")

    # Validate keys match
    missing_formulas = set(complete_structures.keys()) - set(adsorbate_formulas.keys())
    if missing_formulas:
        raise ValueError(f"Missing adsorbate formulas for structures: {missing_formulas}")

    # Get vasp.v2.vasp workflow
    try:
        vasp_wc = WorkflowFactory('vasp.v2.vasp')
    except Exception as e:
        raise RuntimeError(f"Failed to load vasp.v2.vasp workflow plugin: {e}") from e

    vasp_task_cls = task(vasp_wc)

    # Build default SCF inputs
    default_scf_builder = copy.deepcopy(scf_builder_inputs)

    # Add Bader files to retrieve list if requested
    if retrieve_bader_files:
        bader_files = ['CHGCAR', 'AECCAR0', 'AECCAR1', 'AECCAR2']
        if 'settings' not in default_scf_builder:
            default_scf_builder['settings'] = {}
        if isinstance(default_scf_builder['settings'], orm.Dict):
            settings_dict = default_scf_builder['settings'].get_dict()
        else:
            settings_dict = dict(default_scf_builder['settings'])

        if 'ADDITIONAL_RETRIEVE_LIST' not in settings_dict:
            settings_dict['ADDITIONAL_RETRIEVE_LIST'] = []

        for f in bader_files:
            if f not in settings_dict['ADDITIONAL_RETRIEVE_LIST']:
                settings_dict['ADDITIONAL_RETRIEVE_LIST'].append(f)

        default_scf_builder['settings'] = settings_dict

    # Force single-point calculation
    if 'parameters' not in default_scf_builder:
        default_scf_builder['parameters'] = {'incar': {}}
    if 'incar' not in default_scf_builder['parameters']:
        default_scf_builder['parameters']['incar'] = {}
    default_scf_builder['parameters']['incar']['NSW'] = 0
    default_scf_builder['parameters']['incar']['IBRION'] = -1

    # Create the WorkGraph
    wg = WorkGraph(name=name)

    # Set max concurrent jobs
    if max_number_jobs is not None:
        wg.max_number_jobs = max_number_jobs

    # Track tasks for each structure
    separation_tasks = {}
    substrate_scf_tasks = {}
    molecule_scf_tasks = {}
    energy_tasks = {}

    # Add tasks for each structure
    for key, structure in complete_structures.items():
        formula = adsorbate_formulas[key]
        complete_energy = complete_energies[key]

        # Task 1: Separate adsorbate structure
        sep_task = wg.add_task(
            separate_adsorbate_structure,
            name=f'separate_{key}',
            structure=structure,
            adsorbate_formula=orm.Str(formula),
        )
        separation_tasks[key] = sep_task

        # Task 2: SCF for substrate
        substrate_inputs = copy.deepcopy(default_scf_builder)
        substrate_inputs['code'] = code
        # Structure will be linked from separation task output

        substrate_task = wg.add_task(
            vasp_task_cls,
            name=f'scf_substrate_{key}',
            **substrate_inputs,
        )
        # Link substrate structure from separation task
        wg.add_link(sep_task.outputs['substrate'], substrate_task.inputs['structure'])
        substrate_scf_tasks[key] = substrate_task

        # Task 3: SCF for molecule
        molecule_inputs = copy.deepcopy(default_scf_builder)
        molecule_inputs['code'] = code
        # Structure will be linked from separation task output

        molecule_task = wg.add_task(
            vasp_task_cls,
            name=f'scf_molecule_{key}',
            **molecule_inputs,
        )
        # Link molecule structure from separation task
        wg.add_link(sep_task.outputs['molecule'], molecule_task.inputs['structure'])
        molecule_scf_tasks[key] = molecule_task

        # Task 4: Extract substrate energy
        substrate_energy_task = wg.add_task(
            extract_total_energy,
            name=f'extract_substrate_energy_{key}',
        )
        wg.add_link(substrate_task.outputs['misc'], substrate_energy_task.inputs['energies'])

        # Task 5: Extract molecule energy
        molecule_energy_task = wg.add_task(
            extract_total_energy,
            name=f'extract_molecule_energy_{key}',
        )
        wg.add_link(molecule_task.outputs['misc'], molecule_energy_task.inputs['energies'])

        # Task 6: Calculate adsorption energy
        energy_task = wg.add_task(
            calculate_adsorption_energy,
            name=f'adsorption_energy_{key}',
            E_complete=complete_energy,  # From input
        )
        wg.add_link(substrate_energy_task.outputs['result'], energy_task.inputs['E_substrate'])
        wg.add_link(molecule_energy_task.outputs['result'], energy_task.inputs['E_molecule'])
        energy_tasks[key] = energy_task

    return wg

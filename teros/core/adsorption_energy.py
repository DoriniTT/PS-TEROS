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

from .slabs import extract_total_energy


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
) -> dict:
    """
    Build VASP WorkChain inputs from either builder_inputs or old-style parameters.

    This helper provides backward compatibility by accepting either:
    1. New-style: builder_inputs dict (full control, recommended)
    2. Old-style: parameters/options/etc dicts (automatic conversion)

    **Priority Behavior:**
    If both builder_inputs and parameters are provided, builder_inputs takes priority
    and parameters is ignored. This allows explicit control over the input construction.

    Example showing priority:
        >>> # If both are provided, builder_inputs is used
        >>> inputs = _build_vasp_inputs(
        ...     structure=struct,
        ...     code=code,
        ...     builder_inputs={'parameters': {'incar': {'ENCUT': 600}}},
        ...     parameters={'ENCUT': 400},  # This is ignored
        ... )
        >>> inputs['parameters']['incar']['ENCUT']  # Returns 600

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

    Returns:
        Complete input dictionary for VASP WorkChain

    Raises:
        ValueError: If neither builder_inputs nor parameters is provided
    """
    if builder_inputs is not None:
        # Use new-style builder inputs
        # Make a deep copy to avoid modifying original
        inputs = copy.deepcopy(builder_inputs)

        # Override structure and code (always set by workflow)
        inputs['structure'] = structure
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

    # Force SCF mode if requested
    if force_scf:
        # Ensure parameters dict exists and has 'incar' key
        if 'parameters' not in inputs:
            inputs['parameters'] = {'incar': {}}
        if 'incar' not in inputs['parameters']:
            inputs['parameters']['incar'] = {}

        # Override NSW and IBRION for single-point calculation
        inputs['parameters']['incar']['NSW'] = 0
        inputs['parameters']['incar']['IBRION'] = -1

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
) -> dict:
    """
    Internal implementation for separating substrate+adsorbate structure.

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
    from pymatgen.io.ase import AseAtomsAdaptor

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

    # Step 2: Build connectivity graph for ALL atoms using CrystalNN
    sg = StructureGraph.from_local_env_strategy(pmg_structure, CrystalNN())

    # Step 3: Build subgraph containing only edges between candidate atoms
    # This allows us to find adsorbate clusters even if bonded to substrate
    import networkx as nx
    candidate_graph = nx.Graph()
    candidate_graph.add_nodes_from(candidate_indices)

    for i in candidate_indices:
        # Get all bonded neighbors from full graph
        neighbors = sg.get_connected_sites(i)
        for neighbor_info in neighbors:
            j = neighbor_info.index
            # Only add edge if both atoms are candidates (adsorbate elements)
            if j in candidate_indices and i < j:  # avoid duplicate edges
                candidate_graph.add_edge(i, j)

    # Step 4: Find connected components in candidate subgraph
    components = list(nx.connected_components(candidate_graph))

    # Step 5: Find cluster matching adsorbate formula
    matching_clusters = []
    for component in components:
        node_indices = sorted(list(component))

        # Get composition of this cluster
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
            f"Note: Only considering bonds between atoms matching {target_elements}."
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

    # NEW: Relaxation control
    relax_before_adsorption: bool = False,
    relax_builder_inputs: dict | None = None,

    # NEW: SCF control
    scf_builder_inputs: dict | None = None,

    # DEPRECATED (kept for backward compatibility)
    parameters: t.Mapping[str, t.Any] | None = None,
    options: t.Mapping[str, t.Any] | None = None,
    kpoints_spacing: float | None = None,
    clean_workdir: bool = True,
) -> t.Annotated[dict, namespace(
    relaxed_complete_structures=dynamic(orm.StructureData),  # NEW
    separated_structures=dynamic(dict),
    substrate_energies=dynamic(orm.Float),
    molecule_energies=dynamic(orm.Float),
    complete_energies=dynamic(orm.Float),
    adsorption_energies=dynamic(orm.Float),
)]:
    """
    Scatter-gather workflow for calculating adsorption energies.

    This workflow supports two modes:
    1. Direct SCF: separate → SCF (3N jobs)
    2. Relax + SCF: relax → separate → SCF (N+3N jobs)

    Workflow phases:
    - Phase 1 (optional): Relax complete structures using vasp.v2.relax
    - Phase 2: Separate relaxed (or original) structures into substrate/molecule/complete
    - Phase 3: SCF calculations using vasp.v2.vasp for all three components

    Args:
        structures: Dynamic namespace of complete structures (substrate + adsorbate)
        adsorbate_formulas: Dictionary mapping structure keys to adsorbate formulas
                           Example: {'system1': 'OOH', 'system2': 'OH'}
        code: AiiDA code for VASP
        potential_family: Pseudopotential family name
        potential_mapping: Element to potential mapping

        relax_before_adsorption: If True, relax complete structures before separation
        relax_builder_inputs: Full builder dict for vasp.v2.relax (NSW, IBRION, ISIF)
        scf_builder_inputs: Full builder dict for vasp.v2.vasp (NSW=0 enforced)

        parameters: [DEPRECATED] Old-style INCAR parameters dict (for backward compat)
        options: [DEPRECATED] Old-style scheduler options dict
        kpoints_spacing: [DEPRECATED] K-points spacing in Angstrom^-1
        clean_workdir: Whether to clean remote working directories

    Returns:
        Dictionary with namespaces:
        - relaxed_complete_structures: Relaxed structures (empty if relax=False)
        - separated_structures: Dict of separated systems for each input
        - substrate_energies: Energies of bare substrates
        - molecule_energies: Energies of isolated molecules
        - complete_energies: Energies of complete systems
        - adsorption_energies: Final adsorption energies (E_ads = E_complete - E_substrate - E_molecule)

    Example (new API):
        >>> relax_inputs = {
        ...     'parameters': {'incar': {'NSW': 100, 'IBRION': 2, 'ENCUT': 520}},
        ...     'options': {'resources': {'num_machines': 1}},
        ...     'potential_family': 'PBE',
        ...     'potential_mapping': {'La': 'La', 'Mn': 'Mn_pv', 'O': 'O', 'H': 'H'},
        ... }
        >>> scf_inputs = {
        ...     'parameters': {'incar': {'ENCUT': 520}},  # NSW=0 added automatically
        ...     'options': {'resources': {'num_machines': 1}},
        ...     'potential_family': 'PBE',
        ...     'potential_mapping': {'La': 'La', 'Mn': 'Mn_pv', 'O': 'O', 'H': 'H'},
        ... }
        >>> results = compute_adsorption_energies_scatter(
        ...     structures={'oh_lamno3': structure},
        ...     adsorbate_formulas={'oh_lamno3': 'OH'},
        ...     code=code,
        ...     potential_family='PBE',
        ...     potential_mapping={...},
        ...     relax_before_adsorption=True,
        ...     relax_builder_inputs=relax_inputs,
        ...     scf_builder_inputs=scf_inputs,
        ... )

    Example (old API, backward compatible):
        >>> results = compute_adsorption_energies_scatter(
        ...     structures={'oh_ag': structure},
        ...     adsorbate_formulas={'oh_ag': 'OH'},
        ...     code=code,
        ...     potential_family='PBE',
        ...     potential_mapping={'Ag': 'Ag', 'O': 'O', 'H': 'H'},
        ...     parameters={'PREC': 'Accurate', 'ENCUT': 520},
        ...     options={'resources': {'num_machines': 1}},
        ...     kpoints_spacing=0.3,
        ... )
    """
    # Validate inputs
    if structures.keys() != adsorbate_formulas.keys():
        raise ValueError(
            f"Mismatch between structures and adsorbate_formulas keys. "
            f"Structures: {list(structures.keys())}, "
            f"Formulas: {list(adsorbate_formulas.keys())}"
        )

    # Get workflow plugins with error handling
    try:
        vasp_relax = WorkflowFactory('vasp.v2.relax')
        vasp_scf = WorkflowFactory('vasp.v2.vasp')
    except Exception as e:
        raise RuntimeError(
            f"Failed to load VASP workflow plugins. "
            f"Ensure aiida-vasp is installed and registered. Error: {e}"
        ) from e

    relax_task_cls = task(vasp_relax)
    scf_task_cls = task(vasp_scf)

    # ===== PHASE 1: OPTIONAL RELAXATION =====
    # Relax complete structures before separation if requested
    if relax_before_adsorption:
        relaxed_structures: dict[str, orm.StructureData] = {}

        for key, structure in structures.items():
            # Build relaxation inputs
            relax_inputs = _build_vasp_inputs(
                structure=structure,
                code=code,
                builder_inputs=relax_builder_inputs,
                parameters=parameters,
                options=options,
                potential_family=potential_family,
                potential_mapping=potential_mapping,
                kpoints_spacing=kpoints_spacing,
                clean_workdir=clean_workdir,
                force_scf=False,  # Allow relaxation (NSW > 0)
            )

            # Run relaxation
            relax_task = relax_task_cls(**relax_inputs)

            # Validate that relaxation produces structure output
            if not hasattr(relax_task.outputs, 'structure'):
                raise RuntimeError(
                    f"Relaxation workflow for '{key}' did not produce structure output. "
                    f"Check that vasp.v2.relax is configured correctly."
                )

            relaxed_structures[key] = relax_task.outputs.structure

        # Validate we have relaxed structures for all inputs
        if set(relaxed_structures.keys()) != set(structures.keys()):
            missing = set(structures.keys()) - set(relaxed_structures.keys())
            raise RuntimeError(
                f"Relaxation failed to produce structures for: {missing}. "
                f"Check workflow logs for details."
            )

        # Use relaxed structures for separation
        structures_to_separate = relaxed_structures
    else:
        # No relaxation: use original structures
        relaxed_structures = {}
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
    for key, separated in separated_structures_dict.items():
        # Helper function to create SCF inputs for each structure
        def create_scf_inputs(struct: orm.StructureData) -> dict:
            return _build_vasp_inputs(
                structure=struct,
                code=code,
                builder_inputs=scf_builder_inputs,
                parameters=parameters,
                options=options,
                potential_family=potential_family,
                potential_mapping=potential_mapping,
                kpoints_spacing=kpoints_spacing,
                clean_workdir=clean_workdir,
                force_scf=True,  # Enforce NSW=0, IBRION=-1
            )

        # Run SCF calculations for all three components
        scf_components = [
            ('substrate', separated['substrate'], substrate_energies),
            ('molecule', separated['molecule'], molecule_energies),
            ('complete', separated['complete'], complete_energies),
        ]

        for component_name, structure, energy_dict in scf_components:
            scf_calc = scf_task_cls(**create_scf_inputs(structure))
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
        'relaxed_complete_structures': relaxed_structures,  # NEW: Empty dict if relax=False
        'separated_structures': separated_structures_dict,
        'substrate_energies': substrate_energies,
        'molecule_energies': molecule_energies,
        'complete_energies': complete_energies,
        'adsorption_energies': adsorption_energies,
    }

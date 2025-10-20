# teros/core/adsorption_energy.py
"""
Adsorption Energy Calculations Module

This module provides functions to calculate adsorption energies from
substrate+adsorbate structures using AiiDA-WorkGraph and VASP.
"""

from __future__ import annotations

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
        import copy
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
    parameters: t.Mapping[str, t.Any],
    options: t.Mapping[str, t.Any],
    kpoints_spacing: float | None = None,
    clean_workdir: bool = True,
) -> t.Annotated[dict, namespace(
    separated_structures=dynamic(dict),
    substrate_energies=dynamic(orm.Float),
    molecule_energies=dynamic(orm.Float),
    complete_energies=dynamic(orm.Float),
    adsorption_energies=dynamic(orm.Float),
)]:
    """
    Scatter-gather workflow for calculating adsorption energies.

    This workflow:
    1. Separates all substrate+adsorbate structures in parallel
    2. Relaxes all systems (3N VASP jobs) in parallel
    3. Calculates adsorption energies for each system

    Args:
        structures: Dynamic namespace of complete structures
        adsorbate_formulas: Dictionary mapping structure keys to adsorbate formulas
                           Example: {'system1': 'OOH', 'system2': 'OH'}
        code: AiiDA code for VASP
        potential_family: Pseudopotential family name
        potential_mapping: Element to potential mapping
        parameters: VASP INCAR parameters
        options: Scheduler options for VASP calculations
        kpoints_spacing: K-points spacing in Angstrom^-1 (optional)
        clean_workdir: Whether to clean remote working directories

    Returns:
        Dictionary with namespaces:
        - separated_structures: Dict of separated systems for each input
        - substrate_energies: Energies of bare substrates
        - molecule_energies: Energies of isolated molecules
        - complete_energies: Energies of complete systems
        - adsorption_energies: Final adsorption energies
    """
    # Validate inputs
    if structures.keys() != adsorbate_formulas.keys():
        raise ValueError(
            f"Mismatch between structures and adsorbate_formulas keys. "
            f"Structures: {list(structures.keys())}, "
            f"Formulas: {list(adsorbate_formulas.keys())}"
        )

    # Get VASP workflow
    vasp_wc = WorkflowFactory('vasp.v2.vasp')
    vasp_task_cls = task(vasp_wc)

    # Output dictionaries
    separated_dict: dict[str, dict] = {}
    substrate_energies: dict[str, orm.Float] = {}
    molecule_energies: dict[str, orm.Float] = {}
    complete_energies: dict[str, orm.Float] = {}
    adsorption_energies: dict[str, orm.Float] = {}

    # Phase 1: Separate structures (parallel)
    for key, structure in structures.items():
        adsorbate_str = adsorbate_formulas[key]

        separated = separate_adsorbate_structure(
            structure=structure,
            adsorbate_formula=orm.Str(adsorbate_str)
        )

        separated_dict[key] = {
            'substrate': separated.substrate,
            'molecule': separated.molecule,
            'complete': separated.complete,
        }

    # Phase 2: VASP relaxations (parallel, 3N jobs)
    for key, separated in separated_dict.items():
        # Helper function to create VASP inputs
        def create_vasp_inputs(struct: orm.StructureData) -> dict:
            inputs: dict[str, t.Any] = {
                'structure': struct,
                'code': code,
                'parameters': {'incar': dict(parameters)},
                'options': dict(options),
                'potential_family': potential_family,
                'potential_mapping': dict(potential_mapping),
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
            return inputs

        # Substrate relaxation
        substrate_calc = vasp_task_cls(**create_vasp_inputs(separated['substrate']))
        substrate_energies[key] = extract_total_energy(energies=substrate_calc.misc).result

        # Molecule relaxation
        molecule_calc = vasp_task_cls(**create_vasp_inputs(separated['molecule']))
        molecule_energies[key] = extract_total_energy(energies=molecule_calc.misc).result

        # Complete system relaxation
        complete_calc = vasp_task_cls(**create_vasp_inputs(separated['complete']))
        complete_energies[key] = extract_total_energy(energies=complete_calc.misc).result

    # Phase 3: Calculate adsorption energies (parallel)
    for key in structures.keys():
        E_ads = calculate_adsorption_energy(
            E_complete=complete_energies[key],
            E_substrate=substrate_energies[key],
            E_molecule=molecule_energies[key],
        ).result

        adsorption_energies[key] = E_ads

    # Return all results
    return {
        'separated_structures': separated_dict,
        'substrate_energies': substrate_energies,
        'molecule_energies': molecule_energies,
        'complete_energies': complete_energies,
        'adsorption_energies': adsorption_energies,
    }

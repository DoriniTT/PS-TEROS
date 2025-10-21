"""Task functions (CalcFunctions) for surface_hydroxylation module."""

import tempfile
from pathlib import Path

from ase.io import read
from aiida.engine import calcfunction
from aiida.orm import StructureData, Dict
from aiida_workgraph import task

from .utils import aiida_to_ase, ase_to_aiida
from .surface_modes import SurfaceModifier


@calcfunction
def extract_manifest(result: dict) -> Dict:
    """
    Extract manifest from generate_structures result.

    Helper function to extract the manifest Dict from the full result
    returned by generate_structures.

    Args:
        result: Full result dict from generate_structures

    Returns:
        Manifest Dict
    """
    return result['manifest']


@calcfunction
def generate_structures(structure: StructureData, params: Dict) -> dict:
    """
    Generate surface variants using surface_modes.py.

    Args:
        structure: Input relaxed slab structure
        params: Dict with surface_modes parameters:
            - mode: str ('vacancies'/'hydrogen'/'combine')
            - species: str (default 'O')
            - z_window: float (default 0.5)
            - which_surface: str ('top'/'bottom'/'both')
            - oh_dist: float (default 0.98)
            - include_empty: bool (default False)
            - supercell: list[int] or None
            - deduplicate_by_coverage: bool
            - coverage_bins: int or None

    Returns:
        dict with:
            - manifest: Dict (parsed manifest from surface_modes)
            - structure_0, structure_1, ..., structure_N: StructureData (generated variants)
    """
    # Convert AiiDA → ASE
    atoms = aiida_to_ase(structure)

    # Extract parameters with defaults
    p = params.get_dict()
    mode = p.get('mode', 'hydrogen')
    species = p.get('species', 'O')
    z_window = p.get('z_window', 0.5)
    which_surface = p.get('which_surface', 'top')
    oh_dist = p.get('oh_dist', 0.98)
    include_empty = p.get('include_empty', False)
    supercell = p.get('supercell', None)
    deduplicate = p.get('deduplicate_by_coverage', False)
    coverage_bins = p.get('coverage_bins', None)

    # Convert supercell to tuple if provided
    if supercell is not None:
        supercell = tuple(supercell)

    # Create temporary directory for outputs
    with tempfile.TemporaryDirectory() as tmpdir:
        outdir = Path(tmpdir)

        # Create SurfaceModifier instance
        sm = SurfaceModifier(
            atoms=atoms,
            species=species,
            z_window=z_window,
            which_surface=which_surface,
            oh_dist=oh_dist,
            include_empty=include_empty,
            outdir=outdir,
            fmt='vasp',
            supercell=supercell,
            deduplicate_by_coverage=deduplicate,
            coverage_bins=coverage_bins
        )

        # Run appropriate mode
        if mode == 'vacancies':
            manifest = sm.run_vacancies()
        elif mode == 'hydrogen':
            manifest = sm.run_hydrogen()
        elif mode == 'combine':
            manifest = sm.run_combine()
        else:
            raise ValueError(f"Unknown mode: {mode}")

        # Read generated structures
        structures = []
        for variant in manifest['variants']:
            filepath = Path(variant['file'])
            variant_atoms = read(filepath.as_posix())
            structures.append(ase_to_aiida(variant_atoms))

    # Return manifest and structures as namespace
    # AiiDA can't store List of unstored StructureData, so return as dict
    result = {
        'manifest': Dict(dict=manifest),
    }

    # Add each structure as a separate output
    for idx, struct in enumerate(structures):
        result[f'structure_{idx}'] = struct

    return result


@calcfunction
def collect_results(
    manifest: Dict,
    structures: dict,
    energies: dict,
    exit_statuses: dict,
    errors: dict,
) -> dict:
    """
    Collect and organize relaxation results from namespace outputs.

    This function accepts the raw namespace dictionaries from relax_slabs_with_semaphore
    and extracts data from AiiDA nodes to create organized result lists.

    Args:
        manifest: Original manifest Dict from generate_structures
        structures: Namespace dict mapping indices to relaxed StructureData {0: StructureData, 1: ...}
        energies: Namespace dict mapping indices to Float nodes {0: Float, 1: ...}
        exit_statuses: Namespace dict mapping indices to Int nodes {0: Int, 1: ...}
        errors: Namespace dict mapping indices to Str nodes {0: Str, 1: ...}

    Returns:
        dict with:
            - successful_relaxations: Dict with list of successful results
                Each result contains: name, structure_pk, energy, coverage, metadata
            - failed_relaxations: Dict with list of failed results
                Each result contains: name, coverage, exit_status, error_message
            - statistics: Dict with total, succeeded, failed counts

    Note:
        Structures are stored by PK in successful_relaxations.
        To retrieve: orm.load_node(result['successful_relaxations']['results'][i]['structure_pk'])
    """
    from aiida.orm import Int, Float, Str, StructureData

    manifest_dict = manifest.get_dict()
    variants = manifest_dict['variants']

    successful = []
    failed = []

    # Iterate through variants by index
    for idx, variant in enumerate(variants):
        # Convert idx to string for dict lookup (AiiDA Dict only supports string keys)
        idx_key = str(idx)

        if idx_key not in exit_statuses:
            # Relaxation not found (shouldn't happen in normal operation)
            failed.append({
                'name': variant['name'],
                'coverage': variant.get('OH_coverage') or variant.get('vacancy_coverage', 0.0),
                'error_message': 'Relaxation output not found'
            })
            continue

        # Extract values from AiiDA nodes
        exit_status_node = exit_statuses[idx_key]
        exit_status_value = exit_status_node.value if isinstance(exit_status_node, Int) else int(exit_status_node)

        if exit_status_value == 0:
            # Success - extract structure and energy
            if idx_key not in structures or idx_key not in energies:
                # Missing data for successful relaxation
                failed.append({
                    'name': variant['name'],
                    'coverage': variant.get('OH_coverage') or variant.get('vacancy_coverage', 0.0),
                    'error_message': 'Missing structure or energy data'
                })
            else:
                # Extract PK and values from nodes
                structure_node = structures[idx_key]
                energy_node = energies[idx_key]

                structure_pk = structure_node.pk if isinstance(structure_node, StructureData) else int(structure_node)
                energy_value = energy_node.value if isinstance(energy_node, Float) else float(energy_node)

                # Store JSON-serializable data including structure PK
                # Users can load structure later with: orm.load_node(structure_pk)
                successful.append({
                    'name': variant['name'],
                    'structure_pk': structure_pk,
                    'energy': energy_value,
                    'coverage': variant.get('OH_coverage') or variant.get('vacancy_coverage', 0.0),
                    'metadata': variant
                })
        else:
            # Failure - record error
            error_node = errors.get(idx_key, None)
            if error_node:
                error_msg = error_node.value if isinstance(error_node, Str) else str(error_node)
            else:
                error_msg = 'Unknown error'

            failed.append({
                'name': variant['name'],
                'coverage': variant.get('OH_coverage') or variant.get('vacancy_coverage', 0.0),
                'exit_status': exit_status_value,
                'error_message': error_msg
            })

    # Calculate statistics
    statistics = Dict(dict={
        'total': len(variants),
        'succeeded': len(successful),
        'failed': len(failed)
    })

    # Store results as JSON-serializable dicts
    successful_dict = Dict(dict={'results': successful})
    failed_dict = Dict(dict={'results': failed})

    # Return only new nodes (CalcFunctions can't return input nodes)
    return {
        'successful_relaxations': successful_dict,
        'failed_relaxations': failed_dict,
        'statistics': statistics
    }

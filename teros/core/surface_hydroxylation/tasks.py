"""Task functions (CalcFunctions) for surface_hydroxylation module."""

import tempfile
from pathlib import Path
import typing as t

from ase.io import read
from aiida.orm import StructureData, Dict
from aiida_workgraph import task, namespace, dynamic

from .utils import aiida_to_ase, ase_to_aiida
from .surface_modes import SurfaceModifier


@task.calcfunction
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


@task.calcfunction
def extract_successful_relaxations(result: dict) -> Dict:
    """Extract successful_relaxations from collect_results result."""
    return result['successful_relaxations']


@task.calcfunction
def extract_failed_relaxations(result: dict) -> Dict:
    """Extract failed_relaxations from collect_results result."""
    return result['failed_relaxations']


@task.calcfunction
def extract_statistics(result: dict) -> Dict:
    """Extract statistics from collect_results result."""
    return result['statistics']


@task.calcfunction
def generate_structures(
    structure: StructureData,
    params: Dict
) -> t.Annotated[dict, namespace(
    manifest=Dict,
    structures=dynamic(StructureData),
)]:
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
            - structures: dict of StructureData (keyed by index: '0', '1', etc.)
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
        structures_dict = {}
        for idx, variant in enumerate(manifest['variants']):
            filepath = Path(variant['file'])
            variant_atoms = read(filepath.as_posix())
            structures_dict[str(idx)] = ase_to_aiida(variant_atoms)

    # Return manifest and structures as namespace with dynamic structures
    return {
        'manifest': Dict(dict=manifest),
        'structures': structures_dict,
    }


@task.calcfunction
def collect_results(
    manifest,
    structures,
    energies,
) -> t.Annotated[dict, namespace(
    successful_relaxations=Dict,
    failed_relaxations=Dict,
    statistics=Dict,
)]:
    """
    Collect and organize relaxation results from namespace outputs.

    This function accepts the namespace dictionaries from relax_slabs_with_semaphore
    and organizes them into successful and failed results. Failures are detected
    by checking for missing structure/energy data.

    Args:
        manifest: Original manifest Dict from generate_structures
        structures: Namespace dict mapping indices to relaxed StructureData {0: StructureData, 1: ...}
                   Only successful relaxations will have entries
        energies: Namespace dict mapping indices to Float nodes {0: Float, 1: ...}
                 Only successful relaxations will have entries

    Returns:
        dict with:
            - successful_relaxations: Dict with list of successful results
                Each result contains: name, structure_pk, energy, coverage, metadata
            - failed_relaxations: Dict with list of failed results
                Each result contains: name, coverage, error_message
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

        # Check if this relaxation succeeded (has structure and energy outputs)
        if idx_key in structures and idx_key in energies:
            # Success - extract structure and energy
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
            # Failure - missing structure or energy output
            failed.append({
                'name': variant['name'],
                'coverage': variant.get('OH_coverage') or variant.get('vacancy_coverage', 0.0),
                'error_message': 'Relaxation failed - no structure/energy output'
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

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
            - structures: dict of StructureData (keyed by: '0_<variant_name>', '1_<variant_name>', etc.)
                         e.g., '0_oh_000_3_7572', '1_oh_001_7_5145'
                         Note: Dots are replaced with underscores for AiiDA link label compatibility
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
    max_vacancies = p.get('max_vacancies', None)
    max_OH = p.get('max_OH', None)

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
            coverage_bins=coverage_bins,
            max_vacancies=max_vacancies,
            max_OH=max_OH,
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
            # Use descriptive key: index_variantname (e.g., "0_oh_000_3_7572")
            # Replace dots with underscores for AiiDA link label compatibility
            variant_name_safe = variant['name'].replace('.', '_')
            key = f"{idx}_{variant_name_safe}"
            structures_dict[key] = ase_to_aiida(variant_atoms)

    # Return manifest and structures as namespace with dynamic structures
    return {
        'manifest': Dict(dict=manifest),
        'structures': structures_dict,
    }


@task.calcfunction
def _process_results(
    manifest: Dict,
    structure_pks: Dict,
    energy_values: Dict,
) -> t.Annotated[dict, namespace(
    successful_relaxations=Dict,
    failed_relaxations=Dict,
    statistics=Dict,
)]:
    """
    Process results with PKs and values (CalcFunction - provenance tracked).

    This is a helper function called by collect_results @task.graph wrapper.
    """
    manifest_dict = manifest.get_dict()
    variants = manifest_dict['variants']

    pks_dict = structure_pks.get_dict()
    vals_dict = energy_values.get_dict()

    successful = []
    failed = []

    for idx, variant in enumerate(variants):
        # Construct the same descriptive key format used in outputs
        # Format: "idx_variantname" (e.g., "0_oh_000_3_7572")
        # Replace dots with underscores for AiiDA link label compatibility
        variant_name_safe = variant['name'].replace('.', '_')
        key = f"{idx}_{variant_name_safe}"

        if key in pks_dict and key in vals_dict:
            successful.append({
                'name': variant['name'],
                'structure_pk': pks_dict[key],
                'energy': vals_dict[key],
                'coverage': variant.get('OH_coverage') or variant.get('vacancy_coverage', 0.0),
                'metadata': variant
            })
        else:
            failed.append({
                'name': variant['name'],
                'coverage': variant.get('OH_coverage') or variant.get('vacancy_coverage', 0.0),
                'error_message': 'Relaxation failed - no structure/energy output'
            })

    return {
        'successful_relaxations': Dict(dict={'results': successful}),
        'failed_relaxations': Dict(dict={'results': failed}),
        'statistics': Dict(dict={
            'total': len(variants),
            'succeeded': len(successful),
            'failed': len(failed)
        })
    }


# Note: collect_results is now handled directly in workgraph.py
# using helper functions that extract PKs/values from namespace outputs
# and pass them to _process_results CalcFunction.
# This avoids namespace serialization issues.

"""Task functions (CalcFunctions) for surface_hydroxylation module."""

import tempfile
from pathlib import Path

from ase.io import read
from aiida.engine import calcfunction
from aiida.orm import StructureData, Dict

from .utils import aiida_to_ase, ase_to_aiida
from .surface_modes import SurfaceModifier


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

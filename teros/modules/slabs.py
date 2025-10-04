"""
Slab Generation Module

This module provides functions to generate surface slabs from bulk structures
using Pymatgen's SlabGenerator.
"""

from aiida import orm
from aiida_workgraph import task, spec
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from typing import Annotated
from ase import Atoms


@task
def get_slabs(
    relaxed_structure: Atoms,
    miller_indices: list,
    min_slab_thickness: float,
    min_vacuum_thickness: float,
    lll_reduce: bool = True,
    center_slab: bool = True,
    symmetrize: bool = True,
    primitive: bool = True,
    in_unit_planes: bool = False,
    max_normal_search: int = None,
) -> Annotated[dict, spec.namespace(slabs=spec.dynamic(orm.StructureData))]:
    """
    Generate slab structures from a bulk crystal structure using Pymatgen's SlabGenerator.

    This function wraps Pymatgen's slab generation capabilities to produce various
    surface terminations for a given bulk material and Miller index. All slabs are
    generated with symmetric, reduced, c-axis orthogonal cells.

    Args:
        relaxed_structure: ASE Atoms object of the bulk crystal (automatically unwrapped from StructureData)
        miller_indices: List of Miller indices for slab generation (e.g., [1, 0, 0])
        min_slab_thickness: Minimum slab thickness in Angstroms
        min_vacuum_thickness: Minimum vacuum thickness in Angstroms
        lll_reduce: Reduce cell using LLL algorithm before slab generation. Default: False
        center_slab: Center the slab in the c direction of the cell. Default: True
        symmetrize: Generate symmetrically distinct terminations. Default: False
        primitive: Find primitive cell before generating slabs. Default: True
        in_unit_planes: Restrict Miller indices to unit planes. Default: False
        max_normal_search: Max normal search for finding Miller indices. Default: None

    Returns:
        Dictionary with key 'slabs' containing a dict of slab structures.
        Each slab is keyed by termination identifier (e.g., "term_0", "term_1")
        and contains AiiDA StructureData nodes.
    """
    # --- Helper functions for structure conversion ---
    adaptor = AseAtomsAdaptor()

    def get_pymatgen_structure(atoms):
        """Convert ASE Atoms to pymatgen Structure."""
        return adaptor.get_structure(atoms)

    def get_aiida_structure(structure):
        """Convert pymatgen Structure to AiiDA StructureData."""
        ase_atoms = adaptor.get_atoms(structure)
        return orm.StructureData(ase=ase_atoms)

    # --- Convert input structure to pymatgen format ---
    bulk_structure = get_pymatgen_structure(relaxed_structure)

    # --- Slab generation parameters ---
    py_miller_indices = tuple(miller_indices)

    # --- Optionally reduce to primitive cell for cleaner slabs ---
    if primitive:
        analyzer = SpacegroupAnalyzer(bulk_structure)
        bulk_structure = analyzer.get_primitive_standard_structure()

    # --- Create the slab generator object ---
    slab_gen = SlabGenerator(
        bulk_structure,
        py_miller_indices,
        min_slab_thickness,
        min_vacuum_thickness,
        lll_reduce=lll_reduce,
        center_slab=center_slab,
        max_normal_search=max_normal_search,
        in_unit_planes=in_unit_planes,
    )

    # --- Generate all possible slabs for the given orientation ---
    slabs = slab_gen.get_slabs(symmetrize=symmetrize)

    # --- Convert slabs to orthogonal cells and then to AiiDA StructureData ---
    slab_structures = {}

    for i, slab in enumerate(slabs):
        ortho_slab = (
            slab.get_orthogonal_c_slab()
        )  # Convert to orthogonal cell along c-axis
        super_slab = ortho_slab.make_supercell(
            (1, 1, 1)
        )  # No expansion, but could be changed
        slab_structures[f"term_{i}"] = get_aiida_structure(
            super_slab
        )  # Convert to AiiDA StructureData

    # --- Return dictionary with slabs nested in namespace ---
    return {"slabs": slab_structures}

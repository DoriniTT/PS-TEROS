import pytest

pytest.importorskip('aiida')

from pymatgen.analysis.wulff import WulffShape
from pymatgen.core import Lattice

from teros.core.surface_energy.wulff_geometry import extract_wulff_geometry


@pytest.mark.tier1
def test_extract_wulff_geometry_returns_vertices_and_faces():
    lattice = Lattice.cubic(4.0)
    miller_list = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]
    e_surf_list = [1.0, 1.2, 0.8]
    wulff = WulffShape(lattice, miller_list, e_surf_list)

    geometry = extract_wulff_geometry(wulff)

    assert geometry['vertices'], "Expected at least one vertex"
    assert geometry['faces'], "Expected at least one face"
    vertex_count = len(geometry['vertices'])
    assert all(
        all(0 <= idx < vertex_count for idx in face)
        for face in geometry['faces']
    ), "Face indices should reference existing vertices"
    assert geometry['bounding_radius'] > 0

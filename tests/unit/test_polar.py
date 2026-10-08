"""Tier-1 tests for polar surfaces with a pseudo-hydrogen passivated bottom.

The pseudo-hydrogen charges and the formulas follow Zhang et al.,
Sci. Rep. 6, 20055 (2016) and arXiv:1510.08961.
"""

from __future__ import annotations

import pytest

pytest.importorskip("pymatgen")

from pymatgen.core import Lattice, Structure  # noqa: E402

from psteros import polar  # noqa: E402


def _sg(group, lattice, species, coords):
    return Structure.from_spacegroup(group, lattice, species, coords)


def gaas():
    return _sg("F-43m", Lattice.cubic(5.653), ["Ga", "As"], [[0, 0, 0], [0.25, 0.25, 0.25]])


def zno():
    return _sg("P6_3mc", Lattice.hexagonal(3.25, 5.21), ["Zn", "O"], [[1 / 3, 2 / 3, 0], [1 / 3, 2 / 3, 0.382]])


def gan():
    return _sg("P6_3mc", Lattice.hexagonal(3.189, 5.185), ["Ga", "N"], [[1 / 3, 2 / 3, 0], [1 / 3, 2 / 3, 0.377]])


def si():
    return _sg("Fd-3m", Lattice.cubic(5.431), ["Si"], [[0, 0, 0]])


# ---------------------------------------------------------------------------
# Step 4: pseudo-hydrogen model
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("element, charge", [
    ("O", 0.5), ("S", 0.5), ("N", 0.75), ("P", 0.75), ("As", 0.75),
    ("Si", 1.0), ("Ga", 1.25), ("Al", 1.25), ("Zn", 1.5), ("Cd", 1.5),
])
def test_pseudo_hydrogen_charge_is_two_minus_z_over_four(element, charge):
    """The charges used in the papers: 0.5 (S, O), 0.75 (P, As, N), 1.25 (Ga), 1.5 (Zn)."""
    assert polar.pseudo_hydrogen_charge(element) == pytest.approx(charge)


def test_transition_metals_have_no_pseudo_hydrogen():
    with pytest.raises(ValueError, match="no well-defined"):
        polar.pseudo_hydrogen_charge("Fe")


def test_names_for_aiida_and_vasp():
    hydrogens = polar.pseudo_hydrogens(gaas())
    on_as, on_ga = hydrogens["As"], hydrogens["Ga"]
    assert (on_as.charge, on_as.kind_name, on_as.vasp_potential) == (0.75, "H0p75", "H.75")
    assert (on_ga.charge, on_ga.kind_name, on_ga.vasp_potential) == (1.25, "H1p25", "H1.25")
    assert on_as.label == "H0.75(As)"
    assert on_as.site_properties() == {"kind_name": "H0p75", "pseudo_hydrogen": "As"}
    zinc = polar.pseudo_hydrogens(zno())
    assert (zinc["O"].vasp_potential, zinc["Zn"].vasp_potential) == ("H.5", "H1.5")
    silicon = polar.pseudo_hydrogens(si())["Si"]
    assert (silicon.charge, silicon.kind_name, silicon.vasp_potential) == (1.0, "H", "H")


def test_formal_charge_completes_the_electron_count():
    """A pseudo-H stands for a quarter of the missing neighbour: +3/4 on As, -3/4 on Ga."""
    hydrogens = polar.pseudo_hydrogens(gaas())
    assert hydrogens["As"].formal_charge == pytest.approx(0.75)
    assert hydrogens["Ga"].formal_charge == pytest.approx(-0.75)
    assert polar.pseudo_hydrogens(gan())["N"].formal_charge == pytest.approx(0.75)
    assert polar.pseudo_hydrogens(si())["Si"].formal_charge == 0.0


def test_only_tetrahedral_compounds_are_accepted():
    assert polar.tetrahedral_coordination(zno()) == {"Zn": 4, "O": 4}
    rock_salt = _sg("Fm-3m", Lattice.cubic(4.21), ["Mg", "O"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    with pytest.raises(ValueError, match="four-fold"):
        polar.pseudo_hydrogens(rock_salt)


def test_kind_name_reaches_aiida_structures():
    aiida = pytest.importorskip("aiida")  # noqa: F841
    from aiida import orm

    hydrogen = polar.pseudo_hydrogens(gaas())["As"]
    structure = Structure(Lattice.cubic(10.0), ["As", "H"], [[0, 0, 0], [0.15, 0, 0]])
    structure.add_site_property("kind_name", ["As", hydrogen.kind_name])
    node = orm.StructureData(pymatgen=structure)
    assert sorted(node.get_kind_names()) == ["As", "H0p75"]
    assert {kind.name: kind.symbol for kind in node.kinds}["H0p75"] == "H"


# ---------------------------------------------------------------------------
# Step 5: polar slabs on one passivated bottom
# ---------------------------------------------------------------------------

def _check_slab(termination, states):
    """Geometry and bookkeeping every slab must satisfy."""
    import numpy as np

    structure = termination.structure
    normal = structure.lattice.matrix[2] / structure.lattice.c
    heights = structure.cart_coords @ normal
    pseudo = structure.site_properties["pseudo_hydrogen"]
    hydrogen = [i for i, p in enumerate(pseudo) if p]
    atoms = [i for i, p in enumerate(pseudo) if not p]
    assert hydrogen, "bottom is not passivated"
    assert max(heights[hydrogen]) < min(heights[atoms]), "pseudo-H must sit below the slab"
    distances = structure.distance_matrix
    for i in hydrogen:
        nearest = sorted(distances[i][atoms])
        assert nearest[0] < 1.7 and nearest[1] > 2.0, "each pseudo-H bonds to one atom"
        assert structure[int(np.array(atoms)[np.argmin(distances[i][atoms])])].specie.symbol == pseudo[i]
    assert min(distances[i][j] for i in atoms for j in atoms if i < j) > 1.8
    assert set(hydrogen) <= set(termination.bottom_indices)
    assert 0.1 < min(structure.frac_coords[:, 2]) and max(structure.frac_coords[:, 2]) < 0.9


@pytest.mark.parametrize("bulk, miller, top, bottom, hydrogen", [
    (gaas, (1, 1, 1), "Ga", "As", "H0p75"),
    (gaas, (-1, -1, -1), "As", "Ga", "H1p25"),
    (zno, (0, 0, 1), "Zn", "O", "H0p5"),
    (zno, (0, 0, -1), "O", "Zn", "H1p5"),
    (gan, (0, 0, 1), "Ga", "N", "H0p75"),
])
def test_polar_faces_and_their_passivated_bottoms(bulk, miller, top, bottom, hydrogen):
    """(hkl) is the top face; (-h-k-l) puts the other species on top."""
    structure = bulk()
    slabs = polar.find_polar_terminations(structure, miller)
    assert slabs.polar
    assert [t.label for t in slabs] == ["term_0", "term_1"]
    ideal, counted = slabs
    assert (ideal.top_element, ideal.bottom_element) == (top, bottom)
    assert set(ideal.structure.site_properties["kind_name"]) == {top, bottom, hydrogen}
    assert slabs.supercell == (2, 2)  # one vacancy per 2x2 cell satisfies electron counting
    assert ideal.composition == {top: 36, bottom: 36}  # 9 bilayers x 4 cells
    assert ideal.pseudo_hydrogen_counts == {bottom: 4}
    assert not ideal.electron_counting and counted.electron_counting
    assert counted.removed_per_cell == (top,)
    assert counted.composition[top] == 35
    assert ideal.bottom_fingerprint == counted.bottom_fingerprint == slabs.bottom_fingerprint
    states = slabs.oxidation_states
    for termination in slabs:
        _check_slab(termination, states)
        charge = sum(states[e] * n for e, n in termination.composition.items()) + sum(
            slabs.hydrogens[e].formal_charge * n for e, n in termination.pseudo_hydrogen_counts.items()
        )
        assert (abs(charge) < 1e-9) == termination.electron_counting


def test_the_bottom_is_identical_in_every_slab():
    import numpy as np

    slabs = polar.find_polar_terminations(gaas(), (1, 1, 1))
    first, second = slabs
    assert first.structure.lattice == second.structure.lattice
    bottom_first = [(first.structure[i].specie.symbol, tuple(first.structure[i].coords)) for i in first.bottom_indices]
    bottom_second = [(second.structure[i].specie.symbol, tuple(second.structure[i].coords)) for i in second.bottom_indices]
    assert len(bottom_first) == len(bottom_second) > 4
    for (a, x), (b, y) in zip(sorted(bottom_first), sorted(bottom_second)):
        assert a == b and np.allclose(x, y)


def test_nonpolar_validation_slab():
    """GaAs(110) passivated on one side: both species on the bottom, no repair needed."""
    slabs = polar.find_polar_terminations(gaas(), (1, 1, 0))
    assert not slabs.polar and slabs.supercell == (1, 1) and len(slabs) == 1
    termination = slabs[0]
    assert termination.electron_counting and termination.is_stoichiometric
    assert termination.pseudo_hydrogen_counts == {"Ga": 1, "As": 1}
    assert termination.composition == {"Ga": 12, "As": 12}
    _check_slab(termination, slabs.oxidation_states)


def test_thickness_and_options():
    thin = polar.find_polar_terminations(gaas(), (1, 1, 1), bilayers=4, electron_counting=False)
    assert len(thin) == 1 and thin.supercell == (1, 1)
    assert thin[0].composition == {"Ga": 4, "As": 4}
    with pytest.raises(ValueError, match="whole bilayers"):
        polar.find_polar_terminations(gaas(), (1, 1, 1), layers=7, electron_counting=False)
    fixed = polar.find_polar_terminations(gaas(), (1, 1, 1), bilayers=3, supercell=(2, 2), include_ideal=False)
    assert [t.removed_per_cell for t in fixed] == [("Ga",)]


def test_summary_and_files(tmp_path):
    import json

    slabs = polar.find_polar_terminations(gaas(), (1, 1, 1), bilayers=3)
    text = slabs.summary()
    assert text.splitlines()[0] == (
        "GaAs(111) | polar | 2x2 surface cell | top Ga, bottom As + H0.75(As) (shared)"
    )
    assert "Ga-terminated minus Ga per 2x2 cell" in text and "4 H0.75(As)" in text
    assert "<table>" in slabs._repr_html_()
    paths = slabs.write(str(tmp_path))
    summary = json.loads((tmp_path / "terminations.json").read_text())
    assert summary["potcar_order"]["term_0"] == ["Ga", "As", "H.75"]
    poscar = (tmp_path / "term_0_Ga12As12.vasp").read_text().splitlines()
    assert poscar[5].split() == ["Ga", "As", "H"] and poscar[6].split() == ["12", "12", "4"]
    assert len(paths) == 3
    figure = slabs.plot(str(tmp_path / "slabs.png"))
    assert (tmp_path / "slabs.png").stat().st_size > 5000 and figure is not None


def test_polar_slabs_become_aiida_structures_with_pseudo_hydrogen_kinds():
    pytest.importorskip("aiida")
    from aiida import orm

    slab = polar.find_polar_terminations(zno(), (0, 0, 1), bilayers=2, electron_counting=False)[0]
    node = orm.StructureData(pymatgen=slab.structure)
    assert sorted(node.get_kind_names()) == ["H0p5", "O", "Zn"]

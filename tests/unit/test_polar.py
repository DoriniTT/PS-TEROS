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

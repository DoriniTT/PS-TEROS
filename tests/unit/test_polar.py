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


# ---------------------------------------------------------------------------
# Step 6: pseudo chemical potentials
# ---------------------------------------------------------------------------

def test_pseudo_molecule_is_tetrahedral_with_eight_electrons():
    import numpy as np

    molecule = polar.pseudo_molecule(gaas(), "As")
    assert molecule.composition.formula == "As1 H4"
    assert molecule.site_properties["kind_name"] == ["As"] + ["H0p75"] * 4
    distances = molecule.distance_matrix[0][1:]
    assert np.allclose(distances, polar.default_hydrogen_bond_length("As"))
    centre = molecule[0].coords
    vectors = [site.coords - centre for site in molecule[1:]]
    for i in range(4):
        for j in range(i + 1, 4):
            cosine = vectors[i] @ vectors[j] / (np.linalg.norm(vectors[i]) * np.linalg.norm(vectors[j]))
            assert cosine == pytest.approx(-1 / 3)
    assert molecule.lattice.a - 2 * distances[0] >= 15.0 - 1e-9
    # 5 valence electrons of As + 4 x 0.75 from the pseudo-H = 8
    assert 5 + 4 * polar.pseudo_hydrogen_charge("As") == 8
    assert 3 + 4 * polar.pseudo_hydrogen_charge("Ga") == 8


def test_pseudo_molecule_references_follow_eq_8():
    references = polar.PseudoHydrogenReferences.from_pseudo_molecules({"As": -20.0, "Ga": -14.0})
    assert references["As"].value_ev(-4.7) == pytest.approx((-20.0 + 4.7) / 4)
    assert references["Ga"].value_ev(-3.1) == pytest.approx((-14.0 + 3.1) / 4)
    reservoir = references.reservoir_energy_ev({"As": 4}, {"As": -4.7, "Ga": -3.8})
    assert reservoir == pytest.approx(4 * (-20.0 + 4.7) / 4)
    with pytest.raises(ValueError, match="no pseudo chemical potential"):
        references.reservoir_energy_ev({"N": 1}, {"N": -8.0})


@pytest.mark.parametrize("size", [2, 3, 4, 6])
@pytest.mark.parametrize("outer", ["Ga", "As"])
def test_tetrahedral_clusters_match_eq_9_counts(size, outer):
    structure = polar.tetrahedral_cluster(gaas(), outer, size)
    counts = polar.cluster_counts(size)
    inner = "As" if outer == "Ga" else "Ga"
    composition = structure.composition.get_el_amt_dict()
    assert composition[outer] == counts["outer"] and composition.get(inner, 0) == counts["inner"]
    assert composition["H"] == counts["face"] + counts["edge"] + counts["corner"]
    distances = structure.distance_matrix
    for i, site in enumerate(structure):
        if site.specie.symbol != "H":
            assert sum(1 for j in range(len(structure)) if j != i and distances[i][j] < 2.6) == 4
    hydrogen = [i for i, site in enumerate(structure) if site.specie.symbol == "H"]
    if len(hydrogen) > 1:
        assert min(distances[i][j] for i in hydrogen for j in hydrogen if i < j) > 1.5


def test_wurtzite_clusters_use_the_zinc_blende_analogue():
    import numpy as np

    analogue = polar.zinc_blende_analogue(zno())
    bond = min(n.nn_distance for found in zno().get_all_neighbors(4.0) for n in found)
    assert analogue.lattice.a == pytest.approx(4 * bond / np.sqrt(3))
    cluster = polar.tetrahedral_cluster(zno(), "O", 4)
    assert cluster.composition.get_el_amt_dict() == {"O": 20.0, "Zn": 10.0, "H": 40.0}
    assert set(cluster.site_properties["kind_name"]) == {"O", "Zn", "H0p5"}


def test_cluster_fit_recovers_the_parameters_and_is_independent_of_mu():
    face, edge, corner, bulk = -1.2, -1.0, -0.8, -8.5

    def energy(size, mu):
        c = polar.cluster_counts(size)
        return (c["outer"] * mu + c["inner"] * (bulk - mu)
                + c["face"] * face + c["edge"] * edge + c["corner"] * corner)

    for mu in (-3.0, -3.6):
        # muhat shifts by -1/4 per unit of mu: energies built at the shifted values
        shift = -(mu + 3.0) / 4
        energies = {n: energy(n, mu) + shift * sum(polar.cluster_counts(n)[k] for k in ("face", "edge", "corner"))
                    for n in (2, 3, 8, 9)}
        fit = polar.fit_cluster_pseudo_chemical_potentials("Ga", energies, mu)
        assert fit.bulk_energy_ev == pytest.approx(bulk)
        assert (fit.face_ev, fit.edge_ev, fit.corner_ev) == pytest.approx((face + shift, edge + shift, corner + shift))
        assert fit.residual_ev == pytest.approx(0.0, abs=1e-9)
        assert fit.reference.method == "cluster"
        assert fit.reference.value_ev(-3.0) == pytest.approx(face)
    with pytest.raises(ValueError, match="at least four"):
        polar.fit_cluster_pseudo_chemical_potentials("Ga", {2: -1.0, 3: -2.0, 4: -3.0}, -3.0)


def test_doubly_passivated_slab_for_the_eq_7_check():
    import numpy as np

    structure = polar.doubly_passivated_slab(gaas(), (1, 1, 1), bilayers=4)
    kinds = structure.site_properties["kind_name"]
    assert kinds.count("H0p75") == 1 and kinds.count("H1p25") == 1
    normal = structure.lattice.matrix[2] / structure.lattice.c
    heights = structure.cart_coords @ normal
    atoms = [i for i, p in enumerate(structure.site_properties["pseudo_hydrogen"]) if p is None]
    top_h = [i for i, k in enumerate(kinds) if k == "H1p25"]
    assert heights[top_h[0]] > max(heights[atoms])
    assert np.isclose(structure.distance_matrix[top_h[0]][atoms].min(), polar.default_hydrogen_bond_length("Ga"))


# ---------------------------------------------------------------------------
# Step 7: absolute gamma of polar faces in the phase diagram
# ---------------------------------------------------------------------------

import psteros  # noqa: E402

E_GAAS, E_GA, E_AS = -8.5, -3.0, -4.7          # per formula unit / atom
E_MOLECULE = {"As": -16.0, "Ga": -12.0}         # pseudo-molecules As(H.75)4, Ga(H1.25)4


def gaas_references():
    return psteros.BinaryReferences(
        bulk_energy_ev=2 * E_GAAS, bulk_composition={"Ga": 2, "As": 2},
        reference_energies_per_atom_ev={"Ga": E_GA, "As": E_AS},
    )


def test_polar_gamma_follows_sci_rep_eq_5():
    slabs = polar.find_polar_terminations(gaas(), (1, 1, 1), bilayers=3)
    ideal, vacancy = slabs
    energy = -60.0
    termination = psteros.SlabTermination.from_polar(ideal, energy)
    assert termination.surfaces == 1 and termination.face == "GaAs(111)"
    assert termination.pseudo_hydrogen == {"As": 4} and termination.composition == {"Ga": 12, "As": 12}
    hydrogen = psteros.PseudoHydrogenReferences.from_pseudo_molecules(E_MOLECULE)
    references = gaas_references()
    diagram = psteros.surface_phase_diagram([termination], references, pseudo_hydrogen=hydrogen, points=5)
    for point in diagram.curves[termination.label]:
        mu = references.chemical_potentials_ev(point.delta_mu_ev)
        muhat = (E_MOLECULE["As"] - mu["As"]) / 4
        expected = (energy - 12 * mu["Ga"] - 12 * mu["As"] - 4 * muhat) / ideal.area
        assert point.gamma_ev_per_angstrom2 == pytest.approx(expected)
    # A stoichiometric Ga-terminated face still depends on Delta mu through the pseudo-H: slope n_H / (4 A).
    curve = diagram.curves[termination.label]
    slope = (curve[-1].gamma_ev_per_angstrom2 - curve[0].gamma_ev_per_angstrom2) / (
        curve[-1].delta_mu_ev - curve[0].delta_mu_ev)
    assert slope == pytest.approx(4 / (4 * ideal.area))


def test_polar_and_symmetric_slabs_share_one_absolute_scale():
    ideal, vacancy = polar.find_polar_terminations(gaas(), (1, 1, 1), bilayers=3)
    hydrogen = psteros.PseudoHydrogenReferences.from_pseudo_molecules(E_MOLECULE)
    symmetric = psteros.SlabTermination("GaAs(110)", 12 * E_GAAS + 2 * 30.0 * 0.05, {"Ga": 12, "As": 12}, 30.0)
    terminations = [
        psteros.SlabTermination.from_polar(ideal, -60.0),
        psteros.SlabTermination.from_polar(vacancy, -56.5, label="V_Ga 2x2"),
        symmetric,
    ]
    diagram = psteros.surface_phase_diagram(terminations, gaas_references(), pseudo_hydrogen=hydrogen, points=9)
    assert set(diagram.curves) == {"term_0", "V_Ga 2x2", "GaAs(110)"}
    flat = [p.gamma_ev_per_angstrom2 for p in diagram.curves["GaAs(110)"]]
    assert flat == pytest.approx([0.05] * 9)


def test_slabs_of_one_face_must_share_one_bottom():
    hydrogen = psteros.PseudoHydrogenReferences.from_pseudo_molecules(E_MOLECULE)
    thick = polar.find_polar_terminations(gaas(), (1, 1, 1), bilayers=3)[0]
    thin = polar.find_polar_terminations(gaas(), (1, 1, 1), bilayers=2)[0]
    # The same bottom under a thinner slab is the same bottom.
    assert thick.bottom_fingerprint == thin.bottom_fingerprint
    moved = polar.find_polar_terminations(gaas(), (1, 1, 1), bilayers=3, hydrogen_bond_lengths={"As": 1.3})[0]
    small = polar.find_polar_terminations(gaas(), (1, 1, 1), bilayers=3, electron_counting=False)[0]
    assert len({thick.bottom_fingerprint, moved.bottom_fingerprint, small.bottom_fingerprint}) == 3
    mixed = [psteros.SlabTermination.from_polar(thick, -60.0, label="a"),
             psteros.SlabTermination.from_polar(moved, -60.1, label="b")]
    with pytest.raises(ValueError, match="must share one bottom"):
        psteros.surface_phase_diagram(mixed, gaas_references(), pseudo_hydrogen=hydrogen)
    mixed_cells = [psteros.SlabTermination.from_polar(thick, -60.0, label="2x2"),
                   psteros.SlabTermination.from_polar(small, -15.0, label="1x1")]
    with pytest.raises(ValueError, match="must share one bottom"):
        psteros.surface_phase_diagram(mixed_cells, gaas_references(), pseudo_hydrogen=hydrogen)
    other_face = polar.find_polar_terminations(gaas(), (-1, -1, -1), bilayers=3)[0]
    both_faces = [psteros.SlabTermination.from_polar(thick, -60.0, label="A"),
                  psteros.SlabTermination.from_polar(other_face, -61.0, label="B")]
    diagram = psteros.surface_phase_diagram(both_faces, gaas_references(), pseudo_hydrogen=hydrogen)
    assert set(diagram.curves) == {"A", "B"}


def test_passivated_slabs_need_references_and_a_fingerprint():
    ideal = polar.find_polar_terminations(gaas(), (1, 1, 1), bilayers=3)[0]
    termination = psteros.SlabTermination.from_polar(ideal, -60.0)
    with pytest.raises(ValueError, match="need pseudo_hydrogen"):
        psteros.surface_phase_diagram([termination], gaas_references())
    hydrogen = psteros.PseudoHydrogenReferences.from_pseudo_molecules(E_MOLECULE)
    manual = psteros.SlabTermination("x", -60.0, {"Ga": 12, "As": 12}, 20.0, surfaces=1, pseudo_hydrogen={"As": 4})
    with pytest.raises(ValueError, match="no bottom_fingerprint"):
        psteros.surface_phase_diagram([manual], gaas_references(), pseudo_hydrogen=hydrogen)
    with pytest.raises(ValueError, match="surfaces=1"):
        psteros.SlabTermination("y", -60.0, {"Ga": 12, "As": 12}, 20.0, pseudo_hydrogen={"As": 4})
    with pytest.raises(ValueError, match="binary compounds only"):
        refs = psteros.TernaryReferences(-64.0, {"Sr": 2, "Ti": 2, "O": 6}, {"Sr": -1.0, "Ti": -2.0, "O": -4.0})
        psteros.ternary_surface_phase_diagram([manual], refs)


def test_polar_oxide_face_with_oxygen_reference():
    slabs = polar.find_polar_terminations(zno(), (0, 0, 1), bilayers=3)
    hydrogen = psteros.PseudoHydrogenReferences.from_pseudo_molecules({"O": -12.0})
    oxide = psteros.BinaryOxideReferences(-18.0, {"Zn": 2, "O": 2}, -9.8, -1.3)
    diagram = psteros.surface_phase_diagram(
        [psteros.SlabTermination.from_polar(t, -50.0 - i) for i, t in enumerate(slabs)], oxide,
        pseudo_hydrogen=hydrogen, points=5,
    )
    assert diagram.references.variable == "O" and len(diagram.curves) == 2
    mu = oxide.chemical_potentials_ev(-1.0)
    assert mu["Zn"] + mu["O"] == pytest.approx(-9.0)


# ---------------------------------------------------------------------------
# Step 8: bottom check after relaxation
# ---------------------------------------------------------------------------

def _relax(termination, *, bottom_jitter=0.0, shift=(0.0, 0.0, 0.0), top_jitter=0.1, seed=0, reorder=True):
    """A fake relaxation: move the top freely, the bottom by ``bottom_jitter``, shift and reorder sites."""
    import numpy as np

    rng = np.random.default_rng(seed)
    structure = termination.structure.copy()
    bottom = set(termination.bottom_indices)
    for index in range(len(structure)):
        amount = bottom_jitter if index in bottom else top_jitter
        structure.translate_sites([index], rng.normal(0, amount / np.sqrt(3), 3) + np.array(shift),
                                  frac_coords=False, to_unit_cell=False)
    if reorder:
        structure = structure.get_sorted_structure(key=lambda site: site.specie.Z)
    return structure


def test_identical_bottoms_pass_and_a_distorted_bottom_fails():
    slabs = polar.find_polar_terminations(gaas(), (1, 1, 1), bilayers=3)
    ideal, vacancy = slabs
    common = _relax(ideal, bottom_jitter=0.0, shift=(0.1, 0.0, 0.3))
    relaxed = {"term_0": common, "term_1": _relax(vacancy, bottom_jitter=0.0, shift=(0.0, 0.05, -0.2), seed=1)}
    report = polar.check_bottoms(slabs, relaxed)
    assert report.all_passed and report.reference == "term_0"
    assert report[1].rmsd == pytest.approx(0.0, abs=1e-9)
    assert "All bottoms agree." in report.summary()

    relaxed["term_1"] = _relax(vacancy, bottom_jitter=0.15, seed=2)
    report = polar.check_bottoms(slabs, relaxed)
    assert report.failed == ["term_1"]
    assert "FAILED" in report.summary() and "Left out of comparisons: term_1." in report.summary()


def test_changed_cell_fails_and_failed_slabs_are_left_out():
    slabs = polar.find_polar_terminations(gaas(), (1, 1, 1), bilayers=3)
    ideal, vacancy = slabs
    strained = vacancy.structure.copy()
    strained.lattice = strained.lattice.__class__(strained.lattice.matrix * [[1.01], [1.0], [1.0]])
    relaxed = {"term_0": ideal.structure, "term_1": strained}
    report = polar.check_bottoms(slabs, relaxed)
    assert report[1].cell_unchanged is False and report.failed == ["term_1"]
    with pytest.warns(UserWarning, match="left out"):
        kept, _ = polar.polar_slab_terminations(slabs, {"term_0": -60.0, "term_1": -56.0}, relaxed)
    assert [t.label for t in kept] == ["term_0"]


def test_bottom_check_accepts_aiida_structures_and_rejects_mixed_bottoms():
    pytest.importorskip("aiida")
    from aiida import orm

    slabs = polar.find_polar_terminations(gaas(), (1, 1, 1), bilayers=3)
    relaxed = {t.label: orm.StructureData(pymatgen=_relax(t, seed=i)) for i, t in enumerate(slabs)}
    assert polar.check_bottoms(slabs, relaxed).all_passed
    other = polar.find_polar_terminations(gaas(), (-1, -1, -1), bilayers=3)
    with pytest.raises(ValueError, match="one bottom"):
        polar.check_bottoms([slabs[0], other[0]], {"term_0": slabs[0].structure})
    with pytest.raises(ValueError, match="no relaxed structure"):
        polar.check_bottoms(slabs, {"term_0": slabs[0].structure})


# ---------------------------------------------------------------------------
# Step 10: consistency checks
# ---------------------------------------------------------------------------

def test_eq7_check_is_zero_for_consistent_energies():
    references = gaas_references()
    hydrogen = psteros.PseudoHydrogenReferences.from_pseudo_molecules(E_MOLECULE)
    slab = polar.doubly_passivated_slab(gaas(), (1, 1, 1), bilayers=3)
    mu = references.chemical_potentials_ev(0.0)
    consistent = 3 * mu["Ga"] + 3 * mu["As"] + hydrogen.reservoir_energy_ev({"As": 1, "Ga": 1}, mu)
    check = polar.eq7_check(slab, consistent, references, hydrogen)
    assert check.difference_mev_per_angstrom2 == pytest.approx(0.0, abs=1e-9)
    shifted = polar.eq7_check(slab, consistent - 0.1, references, hydrogen)
    assert shifted.difference_mev_per_angstrom2 == pytest.approx(1000 * 0.1 / check.area)
    assert "Eq. 7" in shifted.summary() and "meV/Å²" in shifted.summary()
    # The sum does not depend on the chemical potential.
    assert polar.eq7_check(slab, consistent, references, hydrogen, delta_mu_ev=-0.5).difference_mev_per_angstrom2 == \
        pytest.approx(0.0, abs=1e-9)


def test_nonpolar_check_compares_both_slab_types():
    references = gaas_references()
    hydrogen = psteros.PseudoHydrogenReferences.from_pseudo_molecules(E_MOLECULE)
    one_sided = polar.find_polar_terminations(gaas(), (1, 1, 0))[0]
    area = one_sided.area
    gamma = 0.05
    mu = references.chemical_potentials_ev(0.0)
    symmetric = psteros.SlabTermination("sym", 12 * E_GAAS + 2 * area * gamma, {"Ga": 12, "As": 12}, area)
    energy = 12 * E_GAAS + area * gamma + hydrogen.reservoir_energy_ev(one_sided.pseudo_hydrogen_counts, mu)
    passivated = psteros.SlabTermination.from_polar(one_sided, energy)
    check = polar.nonpolar_check(symmetric, passivated, references, hydrogen)
    assert check.difference_mev_per_angstrom2 == pytest.approx(0.0, abs=1e-9)
    worse = psteros.SlabTermination.from_polar(one_sided, energy + 0.02 * area)
    assert polar.nonpolar_check(symmetric, worse, references, hydrogen).difference_mev_per_angstrom2 == \
        pytest.approx(20.0)

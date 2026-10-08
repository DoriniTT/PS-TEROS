"""Tier-1/2 tests for the charge-neutral calculation set of the v2 workflow."""

from __future__ import annotations

import pytest

pytest.importorskip("pymatgen")

from pymatgen.core import Lattice, Structure  # noqa: E402

import psteros  # noqa: E402
from psteros.surface_study import ChargeNeutralSurfaceStudy, read_qe_results  # noqa: E402


def _sg(group, lattice, species, coords):
    return Structure.from_spacegroup(group, lattice, species, coords)


def mgo():
    return _sg("Fm-3m", Lattice.cubic(4.21), ["Mg", "O"], [[0, 0, 0], [0.5, 0.5, 0.5]])


def mg_metal():
    return _sg("P6_3/mmc", Lattice.hexagonal(3.21, 5.21), ["Mg"], [[1 / 3, 2 / 3, 0.25]])


def o2_molecule():
    return Structure(Lattice.cubic(12.0), ["O", "O"], [[6, 6, 5.4], [6, 6, 6.6]], coords_are_cartesian=True)


def srtio3():
    return _sg("Pm-3m", Lattice.cubic(3.905), ["Sr", "Ti", "O"], [[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.5, 0]])


def srtio3_study():
    return ChargeNeutralSurfaceStudy(
        srtio3(), [(1, 0, 0)],
        references={
            "Sr": _sg("Fm-3m", Lattice.cubic(6.08), ["Sr"], [[0, 0, 0]]),
            "Ti": _sg("P6_3/mmc", Lattice.hexagonal(2.95, 4.69), ["Ti"], [[1 / 3, 2 / 3, 0.25]]),
            "O": o2_molecule(),
        },
        competing_phases={
            "SrO": _sg("Fm-3m", Lattice.cubic(5.16), ["Sr", "O"], [[0, 0, 0], [0.5, 0.5, 0.5]]),
            "TiO2": _sg("P4_2/mnm", Lattice.tetragonal(4.59, 2.96), ["Ti", "O"], [[0, 0, 0], [0.305, 0.305, 0]]),
        },
        min_slab_thickness=8,
    )


def mgo_study(**options):
    values = dict(references={"Mg": mg_metal(), "O": o2_molecule()}, min_slab_thickness=8)
    values.update(options)
    return ChargeNeutralSurfaceStudy(mgo(), [(1, 0, 0), (1, 1, 0)], **values)


PER_ATOM = {"Mg": -1.5, "O": -4.9, "Sr": -1.6, "Ti": -7.8, "Cu": -3.7}


def fake_energies(study, stabilisation=None):
    """Energies from a per-atom model; ``stabilisation`` (eV/atom) makes compounds stable."""

    stabilisation = stabilisation or {"bulk": -1.0, "phase_SrO": -0.9, "phase_TiO2": -0.9}
    energies = {}
    for label, structure in study.structures.items():
        energy = sum(PER_ATOM[site.specie.symbol] for site in structure)
        energies[label] = energy + stabilisation.get(label, 0.0) * len(structure)
    return energies


def test_binary_study_collects_every_calculation():
    study = mgo_study()
    assert study.roles == {
        "bulk": "bulk", "ref_Mg": "reference", "ref_O": "reference",
        "MgO_100_term_0": "slab", "MgO_110_term_0": "slab",
    }
    assert study.slabs == ["MgO_100_term_0", "MgO_110_term_0"]
    assert study.structures["MgO_100_term_0"] is study.termination_sets["MgO_100"][0].structure
    text = str(study)
    assert "MgO: 5 calculations (2 slabs, 2 references)" in text and "MgO(110)" in text


def test_qe_overrides_by_stage():
    study = mgo_study()
    relax = study.qe_overrides("relax")
    assert relax["bulk"].parameters == {} and relax["MgO_100_term_0"].parameters == {}
    assert relax["ref_Mg"].parameters["CONTROL"] == {"calculation": "vc-relax"}
    o2 = relax["ref_O"]
    assert o2.parameters["SYSTEM"] == {"nspin": 2, "tot_magnetization": 2, "starting_magnetization": {"O": 0.5}}
    assert o2.kpoints_distance == 10.0
    static = study.qe_overrides("static")
    assert static["ref_Mg"].parameters == {} and static["ref_O"].parameters["SYSTEM"]["nspin"] == 2
    tuned = study.qe_overrides(extra={"ref_O": psteros.CalculationOverride(parameters={"SYSTEM": {"ecutwfc": 60}})})
    assert tuned["ref_O"].parameters["SYSTEM"]["ecutwfc"] == 60 and tuned["ref_O"].parameters["SYSTEM"]["nspin"] == 2
    with pytest.raises(ValueError, match="stage"):
        study.qe_overrides("scf")
    with pytest.raises(ValueError, match="unknown calculation label"):
        study.qe_overrides(extra={"nothing": psteros.CalculationOverride()})


def test_vasp_overrides_and_potential_mapping():
    study = srtio3_study()
    overrides = study.vasp_overrides()
    assert dict(overrides["bulk"].parameters["INCAR"]) == {"ISIF": 2}
    assert dict(overrides["SrTiO3_100_term_0"].parameters["INCAR"]) == {"ISIF": 2}
    assert dict(overrides["ref_Ti"].parameters["INCAR"]) == {"ISIF": 3}
    assert dict(overrides["phase_SrO"].parameters["INCAR"]) == {"ISIF": 3}
    assert dict(overrides["ref_O"].parameters["INCAR"]) == {"ISIF": 2, "ISPIN": 2, "MAGMOM": [1.0, 1.0]}
    assert study.potential_mapping({"Sr": "Sr_sv", "Ti": "Ti_pv"}) == {"Sr": "Sr_sv", "Ti": "Ti_pv", "O": "O"}


def test_binary_analysis_gives_the_phase_diagram():
    study = mgo_study()
    energies = fake_energies(study)
    result = study.analyse(energies, dict(study.structures), points=5)
    assert set(result.diagram.curves) == set(study.slabs)
    assert result.references.variable == "O"
    assert result.references.reservoir_labels["O"].startswith("$\\frac{1}{2}$")
    text = result.summary()
    assert text.startswith("MgO: 2 terminations, Delta mu_O from") and "eV: MgO_1" in text
    # Without relaxed structures the built slabs are used.
    assert set(study.analyse(energies, points=5).diagram.curves) == set(study.slabs)
    with pytest.raises(ValueError, match="no energy"):
        study.analyse({"bulk": -1.0})


def test_ternary_analysis_uses_the_competing_phases():
    study = srtio3_study()
    assert {label for label, role in study.roles.items() if role == "competing_phase"} == {"phase_SrO", "phase_TiO2"}
    result = study.analyse(fake_energies(study), points=5)
    assert {phase.label for phase in result.references.competing_phases} == {"SrO", "TiO2"}
    assert set(result.diagram.planes) == set(study.slabs)
    assert "stable somewhere" in result.summary()


def test_element_is_its_own_reference():
    copper = _sg("Fm-3m", Lattice.cubic(3.61), ["Cu"], [[0, 0, 0]])
    study = ChargeNeutralSurfaceStudy(copper, [(1, 1, 1)], {}, min_slab_thickness=8)
    assert study.roles == {"bulk": "bulk", "Cu_111_term_0": "slab"}
    energies = {"bulk": -3.7 * len(copper), "Cu_111_term_0": -3.7 * len(study.structures["Cu_111_term_0"]) + 1.0}
    result = study.analyse(energies)
    area = result.terminations[0].surface_area_angstrom2
    assert result.diagram is None
    assert result.surface_energies_j_per_m2["Cu_111_term_0"] == pytest.approx(
        1.0 / (2 * area) * psteros.EV_PER_ANGSTROM2_TO_J_PER_M2)
    assert "Cu_111_term_0" in result.summary()


def test_study_rejects_bad_input():
    with pytest.raises(ValueError, match="lacks a reference"):
        ChargeNeutralSurfaceStudy(mgo(), [(1, 0, 0)], {"Mg": mg_metal()})
    with pytest.raises(ValueError, match="ternary compounds only"):
        mgo_study(competing_phases={"MgO2": mgo()})
    with pytest.raises(ValueError, match="at least one Miller index"):
        ChargeNeutralSurfaceStudy(mgo(), [], {"Mg": mg_metal(), "O": o2_molecule()})


def test_stoichiometric_only_keeps_the_stoichiometric_slabs():
    assert mgo_study(stoichiometric_only=True).slabs == mgo_study().slabs
    study = srtio3_study()
    # Both symmetric SrTiO3(100) slabs (SrO and TiO2 faces) are off-stoichiometric.
    assert not any(t.is_stoichiometric for t in study.termination_sets["SrTiO3_100"])
    with pytest.raises(ValueError, match="no slab to calculate"):
        ChargeNeutralSurfaceStudy(
            srtio3(), [(1, 0, 0)], study.references, competing_phases=study.competing_phases,
            min_slab_thickness=8, stoichiometric_only=True,
        )


def test_read_qe_results_from_both_graph_layouts():
    class Graph:
        outputs = {
            "a_static_parameters": {"energy": -12.5},
            "a_relaxed_structure": "relaxed-a",
            "b_parameters": {"energy": -3.0},
            "b_structure": "relaxed-b",
        }

    energies, structures = read_qe_results(Graph(), ["a", "b"])
    assert energies == {"a": -12.5, "b": -3.0} and structures == {"a": "relaxed-a", "b": "relaxed-b"}
    with pytest.raises(ValueError, match="no output parameters"):
        read_qe_results(Graph(), ["c"])


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_qe_relax_static_workgraph_from_the_study(tmp_path):
    pytest.importorskip("aiida_quantumespresso")
    from aiida import orm

    computer = orm.Computer(
        label=f"local-{tmp_path.name}", hostname="localhost", transport_type="core.local",
        scheduler_type="core.direct", workdir=str(tmp_path),
    ).store()
    orm.InstalledCode(computer=computer, filepath_executable="/bin/true", label="pw",
                      default_calc_job_plugin="quantumespresso.pw").store()
    family = _pseudo_family(tmp_path, ["Mg", "O"])
    study = mgo_study()

    def recipe(control, overrides, name):
        return psteros.SurfaceWorkflowConfig(
            backend="qe",
            calculation=psteros.QeCalculationConfig(
                code_label=f"pw@{computer.label}", pseudo_family=family,
                parameters={"CONTROL": control, "SYSTEM": {"ecutwfc": 40.0}, "ELECTRONS": {"conv_thr": 1e-6}},
            ),
            execution=psteros.ExecutionPolicy(computer=computer.label, queue="debug", max_concurrent_jobs=1),
            name=name, role_overrides=overrides,
        )

    graph = psteros.build_qe_relax_static_workgraph(
        study.structures,
        recipe({"calculation": "relax"}, study.qe_overrides("relax"), "mgo"),
        recipe({"calculation": "scf"}, study.qe_overrides("static"), "mgo_static"),
    )
    relax = graph.tasks["ref_Mg_relax_qe"].inputs.pw.parameters.value.get_dict()
    assert relax["CONTROL"]["calculation"] == "vc-relax"
    slab = graph.tasks["MgO_100_term_0_relax_qe"].inputs.pw.parameters.value.get_dict()
    assert slab["CONTROL"]["calculation"] == "relax"
    o2 = graph.tasks["ref_O_static_qe"].inputs.pw.parameters.value.get_dict()
    assert o2["SYSTEM"]["nspin"] == 2 and o2["CONTROL"]["calculation"] == "scf"


def _pseudo_family(tmp_path, elements):
    """A minimal aiida-pseudo UPF family for the test profile."""

    from aiida_pseudo.groups.family import PseudoPotentialFamily

    directory = tmp_path / "pseudos"
    directory.mkdir()
    for element in elements:
        (directory / f"{element}.upf").write_text(
            f'<UPF version="2.0.1">\n<PP_HEADER element="{element}" z_valence="2.0"/>\n</UPF>\n')
    label = f"test_family_{tmp_path.name}"
    from aiida_pseudo.data.pseudo import UpfData

    PseudoPotentialFamily.create_from_folder(directory, label, pseudo_type=UpfData)
    return label

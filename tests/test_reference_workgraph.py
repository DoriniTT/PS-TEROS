"""Tier-2 tests: VASP reference WorkGraph construction and its calcfunctions.

Graphs are built on the throwaway profile of ``conftest.py`` and never submitted.
"""

from __future__ import annotations

import pytest

pytest.importorskip("aiida_vasp")

import psteros  # noqa: E402

pytestmark = pytest.mark.requires_aiida

INCAR = {"ENCUT": 520, "EDIFF": 1e-7, "ISMEAR": 0, "SIGMA": 0.05, "IBRION": 2, "NSW": 100, "EDIFFG": -0.005}


@pytest.fixture(scope="module")
def code_label() -> str:
    from aiida import orm

    label = "vasp-reference-test"
    try:
        computer = orm.load_computer("localhost-test")
    except Exception:
        computer = orm.Computer(
            label="localhost-test", hostname="localhost", transport_type="core.local", scheduler_type="core.direct"
        ).store()
    try:
        orm.load_code(f"{label}@localhost-test")
    except Exception:
        orm.InstalledCode(
            computer=computer, filepath_executable="/bin/true", default_calc_job_plugin="vasp.vasp", label=label
        ).store()
    return f"{label}@localhost-test"


def recipe(code_label: str, incar=None) -> psteros.SurfaceWorkflowConfig:
    return psteros.SurfaceWorkflowConfig(
        backend="vasp",
        calculation=psteros.VaspCalculationConfig(
            code_label=code_label,
            incar=INCAR if incar is None else incar,
            potential_mapping={"Sn": "Sn_d", "O": "O"},
            kpoints_spacing=0.25,
        ),
        execution=psteros.ExecutionPolicy(computer="localhost-test", queue="debug"),
        name="sno2",
    )


def references() -> dict:
    cell_relax = psteros.CalculationOverride(parameters={"INCAR": {"ISIF": 3}})
    return {
        "o2": psteros.ReferenceSystem(
            psteros.triplet_o2_cell(),
            "gas",
            symmetry_number=2,
            spin=1.0,
            override=psteros.CalculationOverride(
                parameters={"INCAR": {"ISPIN": 2, "NUPDOWN": 2}},
                kpoints_distance=5.0,
                metadata={"max_wallclock_seconds": 600},
            ),
        ),
        "sno2": psteros.ReferenceSystem(
            psteros.rutile_sno2_bulk(), "solid", supercell=(2, 2, 3), block_overrides={"relax": cell_relax}
        ),
    }


def links(workgraph) -> set[tuple[str, str, str, str]]:
    return {
        (link["from_task"], link["from_socket"], link["to_task"], link["to_socket"])
        for link in (link.to_dict() for link in workgraph.links)
    }


def incar(workgraph, task: str) -> dict:
    return workgraph.tasks[task].inputs.parameters.value.get_dict()["incar"]


def test_default_blocks_give_named_tasks_and_outputs(code_label) -> None:
    workgraph = psteros.build_vasp_reference_workgraph(references(), recipe(code_label))
    names = {task.name for task in workgraph.tasks}
    for label in ("o2", "sno2"):
        assert {
            f"{label}_relax_vasp", f"{label}_relax_energy", f"{label}_static_vasp", f"{label}_static_energy",
            f"{label}_vibrations_vasp", f"{label}_vibrations_frequencies",
        } <= names
    assert "sno2_vibrations_supercell" in names
    assert "o2_vibrations_supercell" not in names
    outputs = set(workgraph.outputs._sockets)
    assert {"o2_static_energy", "o2_relax_structure", "o2_vibrations_frequencies", "sno2_vibrations_misc"} <= outputs
    assert workgraph.max_number_jobs == 1


def test_relaxed_structure_feeds_static_and_vibrations(code_label) -> None:
    graph_links = links(psteros.build_vasp_reference_workgraph(references(), recipe(code_label)))
    assert ("o2_relax_vasp", "structure", "o2_static_vasp", "structure") in graph_links
    assert ("o2_relax_vasp", "structure", "o2_vibrations_vasp", "structure") in graph_links
    assert ("sno2_relax_vasp", "structure", "sno2_vibrations_supercell", "structure") in graph_links
    assert ("sno2_vibrations_supercell", "result", "sno2_vibrations_vasp", "structure") in graph_links
    assert ("o2_vibrations_vasp", "retrieved", "o2_vibrations_frequencies", "retrieved") in graph_links


def test_incar_layers_and_reference_overrides(code_label) -> None:
    workgraph = psteros.build_vasp_reference_workgraph(references(), recipe(code_label))
    assert incar(workgraph, "sno2_relax_vasp")["isif"] == 3
    assert incar(workgraph, "sno2_relax_vasp")["nsw"] == 100
    static = incar(workgraph, "sno2_static_vasp")
    assert (static["nsw"], static["ibrion"]) == (0, -1)
    assert "isif" not in static
    vibrations = incar(workgraph, "sno2_vibrations_vasp")
    assert (vibrations["ibrion"], vibrations["nsw"], vibrations["potim"], vibrations["nfree"]) == (6, 1, 0.015, 2)
    assert vibrations["isif"] == 2
    o2 = incar(workgraph, "o2_vibrations_vasp")
    assert (o2["ibrion"], o2["ispin"], o2["nupdown"]) == (5, 2, 2)
    o2_task = workgraph.tasks["o2_vibrations_vasp"]
    assert o2_task.inputs.kpoints_spacing.value.value == 5.0
    assert o2_task.inputs.options.value.get_dict()["max_wallclock_seconds"] == 600
    assert workgraph.tasks["sno2_static_vasp"].inputs.kpoints_spacing.value.value == 0.25
    settings = o2_task.inputs.settings.value.get_dict()
    assert settings["CHECK_IONIC_CONVERGENCE"] is False
    assert settings["parser_settings"] == {"check_ionic_convergence": False}
    assert "OUTCAR" in settings["ADDITIONAL_RETRIEVE_LIST"]
    assert "CHECK_IONIC_CONVERGENCE" not in workgraph.tasks["o2_relax_vasp"].inputs.settings.value.get_dict()


def test_recipe_incar_may_already_use_the_aiida_vasp_namespace(code_label) -> None:
    namespaced = recipe(code_label, {"incar": INCAR})
    workgraph = psteros.build_vasp_reference_workgraph(references(), namespaced)
    assert incar(workgraph, "o2_relax_vasp")["encut"] == 520


def test_custom_blocks_and_structure_from(code_label) -> None:
    blocks = (
        psteros.Relax(name="coarse"),
        psteros.Relax(name="fine", incar={"ediffg": -0.002}),
        psteros.Vibrations(name="phonons", structure_from="fine", nfree=4),
    )
    o2 = {"o2": references()["o2"]}  # sno2 overrides a block named "relax"
    workgraph = psteros.build_vasp_reference_workgraph(o2, recipe(code_label), blocks=blocks)
    graph_links = links(workgraph)
    assert ("o2_coarse_vasp", "structure", "o2_fine_vasp", "structure") in graph_links
    assert ("o2_fine_vasp", "structure", "o2_phonons_vasp", "structure") in graph_links
    assert incar(workgraph, "o2_phonons_vasp")["nfree"] == 4
    assert incar(workgraph, "o2_fine_vasp")["ediffg"] == -0.002


def test_builder_reports_bad_plans_by_name(code_label) -> None:
    no_nsw = {key: value for key, value in INCAR.items() if key != "NSW"}
    with pytest.raises(ValueError, match="needs NSW > 0"):
        psteros.build_vasp_reference_workgraph(references(), recipe(code_label, no_nsw))
    bad = dict(references())
    bad["sn"] = psteros.ReferenceSystem(
        psteros.alpha_sn_bulk(), "solid", block_overrides={"scf": psteros.CalculationOverride()}
    )
    with pytest.raises(ValueError, match="unknown blocks"):
        psteros.build_vasp_reference_workgraph(bad, recipe(code_label))
    with pytest.raises(ValueError, match="invalid reference labels"):
        psteros.build_vasp_reference_workgraph({"o-2": references()["o2"]}, recipe(code_label))


def surface_task(code_label, incar=None, override=None):
    config = recipe(code_label, incar)
    if override is not None:
        config = psteros.SurfaceWorkflowConfig(
            backend="vasp", calculation=config.calculation, execution=config.execution,
            name=config.name, role_overrides={"bulk": override},
        )
    workgraph = psteros.build_surface_workgraph({"bulk": psteros.rutile_sno2_bulk()}, config)
    assert {task.name for task in workgraph.tasks} - {"graph_ctx", "graph_inputs", "graph_outputs"} == {"bulk_vasp"}
    return workgraph.tasks["bulk_vasp"]


def test_surface_builder_puts_a_flat_incar_in_the_aiida_vasp_namespace(code_label) -> None:
    from aiida_vasp.assistant.parameters import ParametersMassage

    parameters = surface_task(code_label).inputs.parameters.value.get_dict()
    assert parameters == {"incar": {key.lower(): value for key, value in INCAR.items()}}
    # aiida-vasp 5 rejects a flat INCAR ("namespace encut is not supported").
    assert ParametersMassage(parameters).parameters.incar["encut"] == 520


def test_surface_builder_keeps_a_namespaced_recipe(code_label) -> None:
    namespaced = {"incar": {"encut": 520, "nsw": 0}, "dynamics": {"positions_dof": [[True] * 3] * 6}}
    parameters = surface_task(code_label, namespaced).inputs.parameters.value.get_dict()
    assert parameters == namespaced


def test_surface_builder_applies_vasp_overrides(code_label) -> None:
    override = psteros.CalculationOverride(
        parameters={"INCAR": {"ISIF": 3}},
        kpoints_distance=0.5,
        settings={"parser_settings": {"add_dos": True}},
        metadata={"max_wallclock_seconds": 600},
    )
    task = surface_task(code_label, {"incar": {"encut": 520, "nsw": 50}}, override)
    assert task.inputs.parameters.value.get_dict() == {"incar": {"encut": 520, "nsw": 50, "isif": 3}}
    assert task.inputs.kpoints_spacing.value.value == 0.5
    settings = task.inputs.settings.value.get_dict()
    assert settings["parser_settings"] == {"add_dos": True}
    assert {"OUTCAR", "vasprun.xml"} <= set(settings["ADDITIONAL_RETRIEVE_LIST"])
    assert task.inputs.options.value.get_dict()["max_wallclock_seconds"] == 600


def test_surface_builder_accepts_psteros_slabs(code_label) -> None:
    slab, _ = psteros.sno2_110_slab(termination="sno", triple_layers=3, vacuum_angstrom=15.0)
    workgraph = psteros.build_surface_workgraph({"slab_sno": slab}, recipe(code_label))
    structure = workgraph.tasks["slab_sno_vasp"].inputs.structure.value
    assert structure.get_pymatgen_structure().composition == slab.composition


def test_surface_builder_defaults_without_override(code_label) -> None:
    task = surface_task(code_label)
    assert task.inputs.kpoints_spacing.value.value == 0.25
    assert set(task.inputs.settings.value.get_dict()) == {"ADDITIONAL_RETRIEVE_LIST"}


OUTCAR = """\
 Eigenvectors and eigenvalues of the dynamical matrix
 ----------------------------------------------------

   1 f  =   47.316210 THz   297.297474 2PiTHz 1578.301557 cm-1   195.683925 meV
   2 f/i=    0.212345 THz     1.334200 2PiTHz    7.083031 cm-1     0.878184 meV

 Eigenvectors after division by SQRT(mass)
"""


def test_calcfunctions_parse_energy_frequencies_and_build_supercells(tmp_path) -> None:
    from aiida import orm

    from psteros.backends.vasp_tasks import make_supercell, vasp_energy, vasp_frequencies

    (tmp_path / "OUTCAR").write_text(OUTCAR + "  free  energy   TOTEN  =       -9.12345600 eV\n")
    retrieved = orm.FolderData(tree=tmp_path)
    assert tuple(vasp_frequencies(retrieved).get_list()) == pytest.approx((1578.301557, -7.083031))

    misc = orm.Dict({"total_energies": {"energy_extrapolated": -9.87, "energy_no_entropy": -9.86}})
    assert vasp_energy(misc, retrieved).value == pytest.approx(-9.87)
    assert vasp_energy(orm.Dict({"total_energies": {}}), retrieved).value == pytest.approx(-9.123456)

    supercell = make_supercell(orm.StructureData(pymatgen=psteros.rutile_sno2_bulk()), orm.List([2, 2, 3]))
    assert len(supercell.sites) == 6 * 12
    assert supercell.get_pymatgen_structure().composition.reduced_formula == "SnO2"

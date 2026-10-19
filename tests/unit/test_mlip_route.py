"""The MLIP route for every study: vibrations from the model, electron-counted slabs, the metal check opt-out."""

from __future__ import annotations

import json
import subprocess
import sys
import warnings

import numpy as np
import pytest

import psteros

pytest.importorskip("pymatgen")

from test_mlip_workflow import (  # noqa: E402
    INSULATOR_DOS,
    METAL_DOS,
    _computer,
    _node_inputs,
    _outputs,
    _recipes,
    _ruo2_study,
    _structures,
)
from test_polar_workflow import fake_energies as polar_energies  # noqa: E402
from test_polar_workflow import gaas_study  # noqa: E402
from test_surface_study import fake_energies, mgo_study, srtio3_study  # noqa: E402

#: The training settings of a graph built without the new options (the keys of 849830f).
DEFAULT_TRAINING_KEYS = {"backend", "name", "device", "metallic_bulk", "frames", "fine_tune", "replay_selection"}


def modes(frequencies, composition, frozen=None, supercell_size=1):
    return psteros.HarmonicVibrations(frequencies, composition=composition, frozen_composition=frozen or {},
                                      supercell_size=supercell_size)


# ------------------------------------------------------------------ free energies of a whole diagram


def test_free_energies_with_vibrations_are_the_manual_ones() -> None:
    energies = {"bulk": -60.0, "slab": -178.0, "ref_Ru": -18.0, "ref_O": -10.0}
    vibrations = {
        "bulk": modes([300.0, 420.0, 510.0], {"Ru": 2, "O": 4}),
        "slab": modes([120.0, 250.0, 640.0, 700.0], {"Ru": 6, "O": 12}, frozen={"Ru": 2, "O": 4}),
        "ref_Ru": modes([180.0, 210.0], {"Ru": 2}),
        "ref_O": modes([1580.0], {"O": 2}),
    }
    result = psteros.free_energies_with_vibrations(energies, vibrations, 600.0, bulk_label="bulk",
                                                   molecule_labels=("ref_O",))
    bulk = vibrations["bulk"]
    assert result == {
        "bulk": psteros.solid_free_energy_ev(-60.0, bulk, 600.0),
        "slab": psteros.solid_free_energy_ev(-178.0, vibrations["slab"], 600.0, bulk=bulk),
        "ref_Ru": psteros.solid_free_energy_ev(-18.0, vibrations["ref_Ru"], 600.0),
        "ref_O": psteros.molecule_reference_energy_ev(-10.0, vibrations["ref_O"]),
    }
    with pytest.raises(ValueError, match=r"no vibrations for \['ref_Ru'\]"):
        psteros.free_energies_with_vibrations(energies, {k: v for k, v in vibrations.items() if k != "ref_Ru"}, 600.0,
                                              bulk_label="bulk")
    with pytest.raises(ValueError, match=r"molecule_labels are not labels of energies_ev: \['o2'\]"):
        psteros.free_energies_with_vibrations(energies, vibrations, 600.0, bulk_label="bulk", molecule_labels=("o2",))
    with pytest.raises(ValueError, match="bulk_label"):
        psteros.free_energies_with_vibrations(energies, vibrations, 600.0, bulk_label="ref_O", molecule_labels=("ref_O",))
    odd = dict(vibrations, slab=modes([120.0], {"Ru": 6, "O": 12}, frozen={"Ru": 2, "O": 3}))
    with pytest.raises(ValueError, match="slab: the frozen region"):
        psteros.free_energies_with_vibrations(energies, odd, 600.0, bulk_label="bulk")


# ------------------------------------------------------------------ the GPU job of the vibrations


EMT_STUB = '''from ase.calculators.emt import EMT


def MACECalculator(model_paths, device, default_dtype, head):
    assert model_paths.startswith("finetune/models/") and head == "Default"
    return EMT()
'''


def _stub_mace(tmp_path):
    """A ``mace.calculators.MACECalculator`` that is ASE's EMT, so the job script runs without a GPU or mace-torch."""

    package = tmp_path / "stub" / "mace"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text("")
    (package / "calculators.py").write_text(EMT_STUB)
    return str(tmp_path / "stub")


def _run_job(folder, stub):
    import os

    environment = dict(os.environ, PYTHONPATH=stub)
    completed = subprocess.run([sys.executable, "-u", "vibrations.py"], cwd=folder, env=environment,
                               capture_output=True, text=True, timeout=300)
    assert completed.returncode == 0, completed.stderr
    return json.loads((folder / "results.json").read_text())


def _cu_slab():
    """Three layers of fcc Cu(100), two atoms per layer (sites 0-1 bottom, 2-3 middle, 4-5 top)."""

    from ase.build import fcc100
    from pymatgen.io.ase import AseAtomsAdaptor

    return AseAtomsAdaptor.get_structure(fcc100("Cu", size=(2, 1, 3), vacuum=6.0, periodic=True))


def _emt_forces(structure, sites, step):
    """The loop of the job script, written out with ASE's EMT."""

    from ase.calculators.emt import EMT
    from pymatgen.io.ase import AseAtomsAdaptor

    atoms = AseAtomsAdaptor.get_atoms(structure)
    atoms.calc = EMT()
    reference = atoms.get_positions()
    plus, minus = [], []
    for site in sites:
        for axis in range(3):
            for sign, forces in ((1.0, plus), (-1.0, minus)):
                positions = reference.copy()
                positions[site, axis] += sign * step
                atoms.set_positions(positions)
                forces.append(atoms.get_forces())
    return np.array(plus), np.array(minus)


def test_vibrations_job_gives_the_forces_of_the_finite_differences(tmp_path) -> None:
    pytest.importorskip("ase")
    slab = _cu_slab()
    sites = [0, 1, 4, 5]  # the middle layer (2, 3) stays fixed
    folder = psteros.prepare_vibrations_job(tmp_path / "job", {"cu": slab}, ("finetune/models/x_seed1.model", "Default"),
                                            {"cu": sites}, displacement_angstrom=0.01, device="cpu")
    plan = json.loads((folder / "job.json").read_text())
    assert plan["structures"] == [{"label": "cu", "file": "cu.extxyz", "sites": sites}]
    assert plan["model"] == {"path": "finetune/models/x_seed1.model", "head": "Default"}
    results = _run_job(folder, _stub_mace(tmp_path))
    assert results["status"] == "completed"
    assert results["structures"]["cu"]["force_evaluations"] == 6 * len(sites) + 1
    data = np.load(folder / "vibrations" / "cu.npz")
    plus, minus = _emt_forces(slab, sites, 0.01)
    assert data["forces_plus"].shape == (3 * len(sites), len(slab), 3)
    assert np.allclose(data["forces_plus"], plus, atol=1e-8) and np.allclose(data["forces_minus"], minus, atol=1e-8)

    with pytest.raises(ValueError, match="every structure"):
        psteros.prepare_vibrations_job(tmp_path / "x", {"cu": slab}, ("m", "Default"), {})
    with pytest.raises(ValueError, match="cu: displaced_sites must be distinct sites"):
        psteros.prepare_vibrations_job(tmp_path / "y", {"cu": slab}, ("m", "Default"), {"cu": [0, 99]})


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_vibrations_calcfunctions_read_the_job_like_the_dft_modes(tmp_path) -> None:
    pytest.importorskip("ase")
    from aiida import orm

    from psteros.backends.mlip_tasks import mlip_harmonic_modes, mlip_vibrations_folder
    from psteros.backends.qe import as_aiida_structure

    slab = _cu_slab()
    node = as_aiida_structure(slab).store()
    sites = [0, 1, 4, 5]
    settings = {"model": ["finetune/models/x_seed1.model", "Default"], "displaced_sites": {"cu": sites},
                "displacement_angstrom": 0.01, "default_dtype": "float64", "device": "cpu"}
    folder = mlip_vibrations_folder._callable(settings=orm.Dict(settings), structure_cu=node)
    assert {"job.json", "cu.extxyz", "vibrations.py"} <= set(folder.list_object_names())
    work = tmp_path / "work"
    work.mkdir()
    for name in folder.list_object_names():
        (work / name).write_text(folder.get_object_content(name))
    _run_job(work, _stub_mace(tmp_path))
    forces = orm.FolderData(tree=str(work / "vibrations")).store()
    options = {"label": "cu", "model": "seed1", "displaced_sites": sites, "displacement_angstrom": 0.01,
               "supercell_size": 1, "molecule": False}
    result = mlip_harmonic_modes._callable(structure=node, settings=orm.Dict(options), forces=forces).get_dict()

    plus, minus = _emt_forces(slab, sites, 0.01)
    expected = psteros.harmonic_vibrations_from_forces(
        masses_amu=[kind.mass for kind in node.kinds] * len(slab), displaced_sites=sites, displacement_angstrom=0.01,
        forces_plus=plus, forces_minus=minus, imaginary_modes="drop")
    assert np.allclose(result["frequencies_cm1"], expected.frequencies_cm1, atol=1e-4)
    assert result["frozen_composition"] == {"Cu": 2} and result["composition"] == {"Cu": 6}
    assert result["zero_modes"] == 0 and result["model"] == "seed1"
    vibrations = psteros.HarmonicVibrations.from_dict(result, imaginary_modes="drop")
    assert vibrations.frozen_composition == {"Cu": 2}

    with pytest.raises(ValueError, match="wrote no ag.npz"):
        mlip_harmonic_modes._callable(structure=node, settings=orm.Dict(dict(options, label="ag")), forces=forces)


# ------------------------------------------------------------------ the builder: old path, vibrations, metal check


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_default_mlip_graph_is_unchanged(tmp_path) -> None:
    """Leaving the new options out builds the graph of before: no vibrations, the metal check on."""

    pytest.importorskip("aiida_cp2k")
    pytest.importorskip("aiida_shell")
    computer = _computer(tmp_path)
    static, fine_tune, job = _recipes(tmp_path, computer)
    structures = _structures()
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        graph = psteros.build_mlip_surface_workgraph(structures, static, fine_tune, job=job, metallic_bulk="bulk")
    names = {task.name for task in graph.tasks}
    assert not any("vibrations" in name or "supercell" in name for name in names)
    assert not any("vibrations" in name for name in _outputs(graph))
    for round_ in (0, 1):
        settings = graph.tasks[f"mlip_round{round_}_training"].inputs.settings.value.get_dict()
        assert set(settings) == DEFAULT_TRAINING_KEYS


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_mlip_graph_adds_the_model_vibrations(tmp_path) -> None:
    pytest.importorskip("aiida_cp2k")
    pytest.importorskip("aiida_shell")
    computer = _computer(tmp_path)
    static, fine_tune, job = _recipes(tmp_path, computer)
    structures = _structures()
    fixed = {"slab": psteros.central_sites(structures["slab"])}
    vibrations = psteros.VibrationsConfig(displacement_angstrom=0.015, supercells={"bulk": (2, 1, 1)},
                                          molecules=("o2",))
    graph = psteros.build_mlip_surface_workgraph(structures, static, fine_tune, job=job, metallic_bulk="bulk",
                                                 fixed_sites=fixed, vibrations=vibrations)
    names = {task.name for task in graph.tasks}
    assert {"mlip_vibrations_job", "mlip_vibrations", "bulk_mlip_supercell"} <= names
    assert {f"{label}_mlip_vibrations" for label in structures} <= names
    assert "slab_mlip_supercell" not in names
    links = {(link.from_task.name, link.to_task.name) for link in graph.links}
    assert ("slab_mlip_round1_relaxed", "mlip_vibrations_job") in links
    assert ("bulk_mlip_round1_relaxed", "bulk_mlip_supercell") in links
    assert ("bulk_mlip_supercell", "mlip_vibrations_job") in links
    assert ("mlip_round1_finetune", "mlip_vibrations") in links  # the model of the last round
    assert ("mlip_vibrations_job", "mlip_vibrations") in links
    assert ("mlip_vibrations", "slab_mlip_vibrations") in links

    settings = graph.tasks["mlip_vibrations_job"].inputs.settings.value.get_dict()
    assert settings["model"] == ["finetune/models/ruo2_round1_seed11.model", "Default"]
    assert settings["displacement_angstrom"] == 0.015
    free = [site for site in range(len(structures["slab"])) if site not in fixed["slab"]]
    assert settings["displaced_sites"]["slab"] == free
    assert settings["displaced_sites"]["bulk"] == list(range(2 * len(structures["bulk"])))
    modes = graph.tasks["o2_mlip_vibrations"].inputs.settings.value.get_dict()
    assert modes["molecule"] is True and modes["model"] == "seed11" and modes["supercell_size"] == 1
    assert graph.tasks["bulk_mlip_vibrations"].inputs.settings.value.get_dict()["supercell_size"] == 2
    job_inputs = graph.tasks["mlip_vibrations"].inputs
    assert dict(job_inputs.filenames.value) == {"finetune": "finetune"}
    assert job_inputs.metadata.options.use_symlinks.value is True
    outputs = _outputs(graph)
    assert {f"{label}_mlip_vibrations" for label in structures} | {"mlip_vibrations_results"} <= outputs

    build = psteros.build_mlip_surface_workgraph
    with pytest.raises(ValueError, match=r"not in structures: \['metal'\]"):
        build(structures, static, fine_tune, job=job, metallic_bulk="bulk",
              vibrations=psteros.VibrationsConfig(molecules=("metal",)))
    with pytest.raises(ValueError, match="supercell displaces every site"):
        build(structures, static, fine_tune, job=job, metallic_bulk="bulk", fixed_sites=fixed,
              vibrations=psteros.VibrationsConfig(supercells={"slab": (1, 1, 1)}))
    with pytest.raises(TypeError, match="VibrationsConfig"):
        build(structures, static, fine_tune, job=job, metallic_bulk="bulk", vibrations={"displacement": 0.01})


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_metal_check_can_be_switched_off_with_a_warning(tmp_path) -> None:
    pytest.importorskip("aiida_cp2k")
    pytest.importorskip("aiida_shell")
    from aiida import orm

    from psteros.backends.mlip_tasks import mlip_training_folder

    computer = _computer(tmp_path)
    static, fine_tune, job = _recipes(tmp_path, computer)
    with pytest.warns(UserWarning, match="electron count"):
        graph = psteros.build_mlip_surface_workgraph(_structures(), static, fine_tune, job=job, metallic_bulk="bulk",
                                                     allow_non_metallic=True)
    assert graph.tasks["mlip_round0_training"].inputs.settings.value.get_dict()["metal_check"] is False
    with pytest.raises(TypeError, match="allow_non_metallic"):
        psteros.build_mlip_surface_workgraph(_structures(), static, fine_tune, job=job, metallic_bulk="bulk",
                                             allow_non_metallic="yes")

    # the calcfunction: an insulating bulk passes only without the check
    structure = psteros.rutile_bulk("Ru", a=4.49, c=3.11)
    node, output, retrieved = _node_inputs(structure, -100.0, INSULATOR_DOS)
    plain = psteros.MaceFineTuneConfig(foundation_model="/models/mh1.model", replay=False, seeds=(5,))
    from dataclasses import asdict

    settings = {"backend": "cp2k", "name": "t_round0", "device": "cuda", "metallic_bulk": "bulk", "frames": ["bulk"],
                "fine_tune": asdict(plain), "replay_selection": {}}
    inputs = {"structure_bulk": node, "output_bulk": output, "retrieved_bulk": retrieved}
    with pytest.raises(ValueError, match="no partly filled states"):
        mlip_training_folder._callable(settings=orm.Dict(settings), **inputs)
    folder = mlip_training_folder._callable(settings=orm.Dict(dict(settings, metal_check=False)), **inputs)
    assert "train.extxyz" in folder.list_object_names()
    metal = _node_inputs(structure, -100.0, METAL_DOS)
    assert mlip_training_folder._callable(
        settings=orm.Dict(settings), structure_bulk=metal[0], output_bulk=metal[1], retrieved_bulk=metal[2])


# ------------------------------------------------------------------ MlipSurfaceStudy with vibrations


def _study_vibrations(study):
    """Made-up modes of every calculation of an MlipSurfaceStudy or an MlipRoute."""

    bulk = getattr(study, "study", study).bulk.composition.get_el_amt_dict()
    result = {"bulk": modes([300.0, 420.0, 510.0], bulk)}
    for label, slab in study.slabs.items():
        frozen = {}
        for index in study.fixed_sites.get(label, ()):
            symbol = slab[index].specie.symbol
            frozen[symbol] = frozen.get(symbol, 0) + 1
        result[label] = modes([90.0 + len(slab), 250.0, 610.0], slab.composition.get_el_amt_dict(), frozen=frozen)
    for label in study.structures:
        if label.startswith("ref_"):
            element = label[len("ref_"):]
            result[label] = modes([1580.0] if element == "O" else [150.0, 190.0], {element: 2})
    return result


def test_mlip_study_analysis_with_vibrations_uses_the_free_energies() -> None:
    study = _ruo2_study()
    energies = {"bulk": -60.0, "ref_O": -10.0, "ref_Ru": -18.0, "o_bridge": -178.0, "o_cus": -190.5, "metal": -164.0}
    vibrations = _study_vibrations(study)
    assert study.molecule_labels == ("ref_O",)
    result = study.analyse(energies, temperature_k=700.0, vibrations=vibrations)
    free = psteros.free_energies_with_vibrations(energies, vibrations, 700.0, bulk_label="bulk",
                                                 molecule_labels=("ref_O",))
    direct = study.analyse(free)
    assert result.temperature_k == 700.0 and direct.temperature_k is None
    assert [(round(mu, 9), a, b) for mu, a, b in result.diagram.transitions] == \
        [(round(mu, 9), a, b) for mu, a, b in direct.diagram.transitions]
    assert result.references.formation_enthalpy_ev == pytest.approx(direct.references.formation_enthalpy_ev)
    assert study.analyse(energies).temperature_k is None  # the old call, unchanged
    with pytest.raises(ValueError, match="only with temperature_k"):
        study.analyse(energies, vibrations=vibrations)
    with pytest.raises(ValueError, match="typed by hand"):
        study.analyse(energies, temperature_k=700.0)
    assert study.calculation_count(vibrations=True) == {"dft_statics": 18, "gpu_jobs": 5}


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_mlip_study_builds_the_vibrations_and_the_metal_opt_out(tmp_path) -> None:
    pytest.importorskip("aiida_cp2k")
    pytest.importorskip("aiida_shell")
    computer = _computer(tmp_path)
    static, fine_tune, job = _recipes(tmp_path, computer)
    static = psteros.SurfaceWorkflowConfig(backend="cp2k", name="ruo2", execution=static.execution,
                                           calculation=static.calculation)
    study = _ruo2_study()
    graph = study.build(static, fine_tune, job=job, vibrations=psteros.VibrationsConfig(supercells={"bulk": (1, 1, 2)}))
    modes_of_o2 = graph.tasks["ref_O_mlip_vibrations"].inputs.settings.value.get_dict()
    assert modes_of_o2["molecule"] is True  # the gas reference is a molecule without being named
    assert graph.tasks["ref_Ru_mlip_vibrations"].inputs.settings.value.get_dict()["molecule"] is False
    with pytest.warns(UserWarning, match="electron count"):
        study.build(static, fine_tune, job=job, allow_non_metallic=True)
    off_centre = _ruo2_study(fixed_half_width_angstrom=0.4)  # freezes Ru2O2, not whole RuO2 units
    with pytest.raises(ValueError, match=r"frozen centre of o_bridge \(O2Ru2\).* bulk stoichiometry RuO2"):
        off_centre.build(static, fine_tune, job=job, vibrations=psteros.VibrationsConfig())
    off_centre.build(static, fine_tune, job=job)  # without vibrations the centre is free to choose


# ------------------------------------------------------------------ MlipRoute: the studies with electron-counted slabs


def test_mlip_route_runs_a_charge_neutral_study() -> None:
    study = mgo_study()
    route = psteros.MlipRoute(study)
    assert route.structures == study.structures
    assert list(route.slabs) == ["MgO_100_term_0", "MgO_110_term_0"]
    assert route.fixed_sites == {label: psteros.central_sites(slab) for label, slab in route.slabs.items()}
    assert route.molecule_labels == ("ref_O",)
    assert route.calculation_count() == {"dft_statics": 15, "gpu_jobs": 4}
    assert route.calculation_count(vibrations=True)["gpu_jobs"] == 5
    assert "MgO: 2 slabs (charge-neutral)" in route.summary()
    rattled = psteros.MlipRoute(study, rattled=2, fixed_half_width_angstrom=None)
    assert list(rattled.training_structures) == ["MgO_100_term_0_r0", "MgO_100_term_0_r1",
                                                 "MgO_110_term_0_r0", "MgO_110_term_0_r1"]
    assert rattled.fixed_sites == {}
    with pytest.raises(TypeError, match="MlipSurfaceStudy"):
        psteros.MlipRoute(_ruo2_study())
    with pytest.raises(ValueError, match="rattled"):
        psteros.MlipRoute(study, rattled=-1)

    energies = fake_energies(study)
    result = route.analyse(energies)
    direct = study.analyse(energies)
    assert isinstance(result, psteros.MlipStudyResult) and result.energy_source == "mlip_dft"
    assert [(round(mu, 9), a, b) for mu, a, b in result.diagram.transitions] == \
        [(round(mu, 9), a, b) for mu, a, b in direct.diagram.transitions]
    with pytest.raises(ValueError, match="no energy for"):
        route.analyse({label: value for label, value in energies.items() if label != "bulk"})

    vibrations = _study_vibrations(route)
    hot = route.analyse(energies, temperature_k=500.0, vibrations=vibrations)
    free = psteros.free_energies_with_vibrations(energies, vibrations, 500.0, bulk_label="bulk",
                                                 molecule_labels=("ref_O",))
    assert hot.temperature_k == 500.0
    assert [(round(mu, 9), a, b) for mu, a, b in hot.diagram.transitions] == \
        [(round(mu, 9), a, b) for mu, a, b in study.analyse(free).diagram.transitions]


def test_mlip_route_runs_ternary_and_polar_studies() -> None:
    ternary = srtio3_study()
    route = psteros.MlipRoute(ternary)
    energies = fake_energies(ternary)
    assert route.analyse(energies).diagram.regions.keys() == ternary.analyse(energies).diagram.regions.keys()

    polar = gaas_study(faces=[(1, 1, 1)], bilayers=2, eq7_check=False)
    route = psteros.MlipRoute(polar)
    assert route.fixed_sites == {}  # every atom relaxes, as in the polar method
    assert set(route.molecule_labels) == {label for label, role in polar.roles.items() if role == "pseudo_molecule"}
    assert "(polar, passivated bottom)" in route.summary()
    energies = polar_energies(polar)
    result = route.analyse(energies)
    direct = polar.analyse(energies, dict(polar.structures))
    assert isinstance(result, psteros.PolarStudyResult)
    assert set(result.diagram.curves) == set(direct.diagram.curves)
    assert all(report.all_passed for report in result.bottom_checks.values())
    with pytest.raises(ValueError, match="no vibrations yet"):
        route.analyse(energies, temperature_k=300.0, vibrations={})


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_mlip_route_builds_without_the_metal_check(tmp_path) -> None:
    pytest.importorskip("aiida_cp2k")
    pytest.importorskip("aiida_shell")
    from test_cp2k_backend import cp2k_config

    computer = _computer(tmp_path)
    _, fine_tune, job = _recipes(tmp_path, computer)
    fine_tune = psteros.MaceFineTuneConfig(foundation_model="/models/mace-mh-1.model", replay=False, seeds=(3,))
    study = mgo_study()
    route = psteros.MlipRoute(study, rattled=1)
    slab = next(iter(route.slabs))
    calculation = cp2k_config(f"cp2k@{computer}", basis_sets={"Mg": "DZVP-MOLOPT-SR-GTH-q10", "O": "DZVP-MOLOPT-PBE-GTH-q6"},
                              potentials={"Mg": "GTH-PBE-q10", "O": "GTH-PBE-q6"})
    static = psteros.SurfaceWorkflowConfig(
        backend="cp2k", name="mgo", calculation=calculation,
        execution=psteros.ExecutionPolicy(computer=computer, queue="standard", max_concurrent_jobs=2),
        role_overrides={"ref_O": psteros.CalculationOverride(kpoints_distance=100.0),
                        slab: psteros.CalculationOverride(kpoints_distance=0.5)},
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)  # electron-counted slabs: no warning
        graph = route.build(static, fine_tune, job=job, vibrations=psteros.VibrationsConfig(supercells={"bulk": (2, 2, 2)}))
    settings = graph.tasks["mlip_round0_training"].inputs.settings.value.get_dict()
    assert settings["metal_check"] is False and settings["metallic_bulk"] == "bulk"
    assert f"{slab}_r0" in settings["frames"]
    # the rattled copy of a slab has the k-points of the slab's override
    from psteros.backends.cp2k import cp2k_kpoints_mesh

    expected = list(cp2k_kpoints_mesh(route.slabs[slab].lattice.matrix, 0.5))
    for name in (slab, f"{slab}_r0"):
        assert graph.tasks[f"{name}_static_cp2k"].inputs.cp2k.kpoints.value.get_kpoints_mesh()[0] == expected
    assert expected != list(cp2k_kpoints_mesh(route.slabs[slab].lattice.matrix, calculation.kpoints_spacing))
    plan = json.loads(graph.tasks["mlip_round0_relax"].inputs.nodes.job.value.get_object_content("job.json"))
    fixed = {item["label"]: item["fixed"] for item in plan["structures"]}
    assert fixed[slab] == list(route.fixed_sites[slab]) and fixed["bulk"] == []
    assert graph.tasks["ref_O_mlip_vibrations"].inputs.settings.value.get_dict()["molecule"] is True

    polar = psteros.MlipRoute(gaas_study(faces=[(1, 1, 1)], bilayers=2, eq7_check=False))
    with pytest.raises(ValueError, match="runs with VASP"):
        polar.build(static, fine_tune, job=job)


# ------------------------------------------------------------------ DFT vibrations with CP2K, to check the model's


def _cp2k_forces(forces) -> str:
    rows = [f"      {i + 1}      1      O   {f[0]:.12f}   {f[1]:.12f}   {f[2]:.12f}" for i, f in enumerate(forces)]
    return "\n".join([" ATOMIC FORCES in [eV/Angstrom]", "", " # Atom   Kind   Element   X   Y   Z", *rows,
                      " SUM OF ATOMIC FORCES  0.0 0.0 0.0 0.0"]) + "\n"


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_cp2k_vibrations_graph_and_modes(tmp_path) -> None:
    pytest.importorskip("aiida_cp2k")
    import math

    from aiida import orm
    from test_vibrations import CM1_PER_SQRT_EV_PER_A2_AMU, _folder, diatomic_constants, displaced_forces

    from psteros.backends.vibration_tasks import harmonic_modes
    from psteros.vibrations_workflow import _vibrations_workgraph, displacement_key

    computer = _computer(tmp_path)
    static, _, _ = _recipes(tmp_path, computer, run_type="ENERGY_FORCE")
    o2 = psteros.triplet_o2_cell(cell_length=10.0)
    with pytest.raises(ValueError, match="CP2K is not supported yet"):  # the public builder, unchanged
        psteros.build_vibrations_workgraph({"o2": o2}, static)
    graph = _vibrations_workgraph({"o2": o2}, static, psteros.VibrationsConfig(molecules=("o2",)), submit=False,
                                  allow_cp2k=True)
    names = {task.name for task in graph.tasks}
    assert {f"o2_vib_{displacement_key(site, axis, sign)}_cp2k" for site in (0, 1) for axis in range(3)
            for sign in ("plus", "minus")} <= names
    assert "o2_vibrations" in names
    relaxing, _, _ = _recipes(tmp_path, computer, run_type="GEO_OPT")
    with pytest.raises(ValueError, match="o2: the displaced CP2K calculations need RUN_TYPE ENERGY_FORCE"):
        _vibrations_workgraph({"o2": o2}, relaxing, None, submit=False, allow_cp2k=True)

    spring, mass, step = 70.0, 15.999, 0.01
    molecule = orm.StructureData(cell=[[10.0, 0, 0], [0, 10.0, 0], [0, 0, 10.0]])
    molecule.append_atom(position=(5.0, 5.0, 5.0), symbols="O")
    molecule.append_atom(position=(6.2, 5.0, 5.0), symbols="O")
    plus, minus = displaced_forces(diatomic_constants(spring), [0, 1], step)
    retrieved = {}
    for index, (site, axis) in enumerate((s, a) for s in (0, 1) for a in range(3)):
        retrieved[displacement_key(site, axis, "plus")] = _folder("psteros_forces.dat", _cp2k_forces(plus[index]))
        retrieved[displacement_key(site, axis, "minus")] = _folder("psteros_forces.dat", _cp2k_forces(minus[index]))
    settings = orm.Dict({"backend": "cp2k", "displaced_sites": [0, 1], "displacement_angstrom": step,
                         "supercell_size": 1, "molecule": True})
    vibrations = psteros.HarmonicVibrations.from_dict(
        harmonic_modes._callable(structure=molecule, settings=settings, **retrieved))
    assert vibrations.frequencies_cm1 == pytest.approx(
        (math.sqrt(spring / (mass / 2)) * CM1_PER_SQRT_EV_PER_A2_AMU,), rel=1e-3)


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_dft_check_of_the_vibrations_keeps_the_model_fixed_sites(tmp_path) -> None:
    pytest.importorskip("aiida_cp2k")
    from psteros.mlip_workflow import _dft_vibrations

    computer = _computer(tmp_path)
    static, _, _ = _recipes(tmp_path, computer, run_type="ENERGY_FORCE")
    study = _ruo2_study()
    static = psteros.SurfaceWorkflowConfig(backend="cp2k", name="ruo2", execution=static.execution,
                                           calculation=static.calculation)
    structures = {label: study.structures[label] for label in ("o_bridge", "ref_O")}
    graph = _dft_vibrations(structures, static, None, {"o_bridge": study.fixed_sites["o_bridge"]},
                            study.molecule_labels, submit=False)
    settings = graph.tasks["o_bridge_vibrations"].inputs.settings.value.get_dict()
    assert not set(settings["displaced_sites"]) & set(study.fixed_sites["o_bridge"])
    assert len(settings["displaced_sites"]) == len(study.slabs["o_bridge"]) - len(study.fixed_sites["o_bridge"])
    assert graph.tasks["ref_O_vibrations"].inputs.settings.value.get_dict()["molecule"] is True


def test_half_width_of_the_fixed_centre_per_slab() -> None:
    study = mgo_study()
    route = psteros.MlipRoute(study, fixed_half_width_angstrom={"MgO_100_term_0": 2.5, "MgO_110_term_0": None})
    assert route.fixed_sites == {"MgO_100_term_0": psteros.central_sites(route.slabs["MgO_100_term_0"], 2.5)}
    assert psteros.MlipRoute(study, fixed_half_width_angstrom=1.5).fixed_sites == psteros.MlipRoute(study).fixed_sites
    with pytest.raises(ValueError, match=r"not slabs: \['bulk'\]"):
        psteros.MlipRoute(study, fixed_half_width_angstrom={"bulk": 1.0})
    with pytest.raises(ValueError, match="positive"):
        psteros.MlipRoute(study, fixed_half_width_angstrom={"MgO_100_term_0": 0.0})
    ruo2 = _ruo2_study(fixed_half_width_angstrom={"o_cus": 0.4})
    assert set(ruo2.fixed_sites) == {"o_cus"}
    assert _ruo2_study(fixed_half_width_angstrom=1.5).fixed_sites == _ruo2_study().fixed_sites  # the old call
    with pytest.raises(ValueError, match="not slabs"):
        _ruo2_study(fixed_half_width_angstrom={"ref_O": 1.0})


# ------------------------------------------------------------------ the examples


def _example(name):
    import importlib.util
    from pathlib import Path

    path = Path(__file__).resolve().parents[2] / "examples" / "mlip_routes" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"example_{name}", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_mlip_route_examples_show_their_routes() -> None:
    sno2 = _example("sno2_charge_neutral_mlip").route()
    assert len(sno2.slabs) == 3 and sno2.molecule_labels == ("ref_O",)
    from psteros.mlip_workflow import _check_frozen_regions

    _check_frozen_regions(sno2.slabs, sno2.fixed_sites, sno2.study.bulk)  # whole SnO2 units, as the vibrations need
    gaas = _example("gaas_polar_mlip").route()
    assert gaas.polar and gaas.fixed_sites == {} and len(gaas.slabs) == 4


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_mlip_route_examples_build_their_graphs(tmp_path) -> None:
    pytest.importorskip("aiida_vasp")
    pytest.importorskip("aiida_shell")
    import argparse

    from aiida import orm
    from ase import Atoms
    from ase.io import write

    computer = orm.Computer(label=f"local-{tmp_path.name}", hostname="localhost", transport_type="core.local",
                            scheduler_type="core.direct", workdir=str(tmp_path)).store()
    orm.InstalledCode(computer=computer, filepath_executable="/bin/true", label="vasp",
                      default_calc_job_plugin="vasp.vasp").store()
    orm.InstalledCode(computer=computer, filepath_executable="/opt/mace/bin/python", label="mace",
                      default_calc_job_plugin="core.shell").store()
    frames = [Atoms(symbols, positions=[[0, 0, 0], [2.0 + 0.01 * i, 0, 0]], cell=[8, 8, 8], pbc=True)
              for i in range(6) for symbols in ("SnO", "GaAs")]
    write(tmp_path / "replay.xyz", frames, format="extxyz")
    args = argparse.Namespace(code=f"vasp@{computer.label}", mace_code=f"mace@{computer.label}",
                              model="/models/mace-mh-1.model", replay=str(tmp_path / "replay.xyz"), queue="gpu",
                              vibrations=True, submit=False)

    example = _example("sno2_charge_neutral_mlip")
    route = example.route()
    static, fine_tune, job = example.recipes(args)
    graph = route.build(static, fine_tune, job=job, vibrations=example.vibrations())
    names = {task.name for task in graph.tasks}
    assert {"mlip_vibrations", "bulk_mlip_supercell", "ref_Sn_mlip_supercell"} <= names
    assert graph.tasks["mlip_round0_training"].inputs.settings.value.get_dict()["metal_check"] is False

    example = _example("gaas_polar_mlip")
    route = example.route()
    static, fine_tune, job = example.recipes(args, route.study)
    graph = route.build(static, fine_tune, job=job)
    slab = graph.tasks["GaAs_111_term_0_static_vasp"].inputs
    assert "H0p75" in slab.structure.value.get_kind_names() or "H1p25" in slab.structure.value.get_kind_names()
    assert slab.parameters.value.get_dict()["incar"]["ldipol"] is True
    plan = json.loads(graph.tasks["mlip_round0_relax"].inputs.nodes.job.value.get_object_content("job.json"))
    assert all(item["fixed"] == [] for item in plan["structures"])  # every atom of a polar slab relaxes


# ------------------------------------------------------------------ reading a finished graph


class _FinishedGraph:
    """What the readers use of a finished MLIP graph: called processes by task name, and the graph outputs."""

    def __init__(self, energies, structures, vibrations):
        from types import SimpleNamespace

        def node(pk, **outputs):
            return SimpleNamespace(pk=pk, is_finished_ok=True, outputs=SimpleNamespace(**outputs))

        children = {"mlip_round0_relax": node(1), "mlip_round1_relax": node(2)}
        for index, label in enumerate(energies):
            children[f"{label}_mlip_dft_energy"] = node(10 + index, result=SimpleNamespace(value=energies[label]))
            children[f"{label}_mlip_round1_relaxed"] = node(100 + index, structure=structures[label])
        links = [SimpleNamespace(link_label=name, node=child) for name, child in children.items()]
        self.base = SimpleNamespace(links=SimpleNamespace(get_outgoing=lambda link_type=None: SimpleNamespace(
            all=lambda: links)))
        self.outputs = {f"{label}_mlip_vibrations": modes.to_dict() for label, modes in vibrations.items()}


def test_analysis_reads_energies_structures_and_modes_of_the_graph() -> None:
    study = mgo_study()
    route = psteros.MlipRoute(study)
    energies = fake_energies(study)
    vibrations = _study_vibrations(route)
    graph = _FinishedGraph(energies, study.structures, vibrations)

    read_energies, relaxed = route.read(graph)
    assert read_energies == energies and set(relaxed) == set(study.structures)
    cold = route.analyse(graph)
    assert [t[:3] for t in cold.diagram.transitions] == [t[:3] for t in route.analyse(energies).diagram.transitions]
    hot = route.analyse(graph, temperature_k=600.0)  # the modes of the graph
    typed = route.analyse(energies, temperature_k=600.0, vibrations=vibrations)
    assert hot.temperature_k == 600.0
    assert [(round(mu, 9), a, b) for mu, a, b in hot.diagram.transitions] == \
        [(round(mu, 9), a, b) for mu, a, b in typed.diagram.transitions]

    ruo2 = _ruo2_study()
    ruo2_energies = {"bulk": -60.0, "ref_O": -10.0, "ref_Ru": -18.0, "o_bridge": -178.0, "o_cus": -190.5,
                     "metal": -164.0}
    ruo2_vibrations = _study_vibrations(ruo2)
    graph = _FinishedGraph(ruo2_energies, ruo2.structures, ruo2_vibrations)
    result = ruo2.analyse(graph, temperature_k=600.0)
    expected = ruo2.analyse(ruo2_energies, temperature_k=600.0, vibrations=ruo2_vibrations)
    assert result.references.formation_enthalpy_ev == pytest.approx(expected.references.formation_enthalpy_ev)
    del graph.outputs["metal_mlip_vibrations"]
    with pytest.raises(ValueError, match="no vibrations output for metal"):
        ruo2.analyse(graph, temperature_k=600.0)


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_vibrations_job_merges_stderr_into_stdout(tmp_path) -> None:
    """The vibrations job also writes MACE's stderr into stdout (aiida-shell exit 410 otherwise)."""

    pytest.importorskip("aiida_cp2k")
    pytest.importorskip("aiida_shell")
    computer = _computer(tmp_path)
    static, fine_tune, job = _recipes(tmp_path, computer)
    graph = psteros.build_mlip_surface_workgraph(_structures(), static, fine_tune, job=job, metallic_bulk="bulk",
                                                 vibrations=psteros.VibrationsConfig())
    assert graph.tasks["mlip_vibrations"].inputs.metadata.options.redirect_stderr.value is True

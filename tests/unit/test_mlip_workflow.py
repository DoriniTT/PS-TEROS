"""The MLIP route as a WorkGraph: DFT labels of each backend, the calcfunctions and the builder."""

from __future__ import annotations

import json

import numpy as np
import pytest

import psteros
from psteros import mlip
from psteros.backends.cp2k import BOHR_TO_ANGSTROM, HARTREE_TO_EV

pymatgen = pytest.importorskip("pymatgen")

FORCES_AU = """
 ATOMIC FORCES in [a.u.]

 # Atom   Kind   Element          X                      Y                      Z
      1      1      Ru          0.0100000000000000     0.0000000000000000    -0.0200000000000000
      2      2      O          -0.0100000000000000     0.0050000000000000     0.0200000000000000
 SUM OF ATOMIC FORCES           0.0000000000000000     0.0050000000000000     0.0000000000000000     0.0050000
"""


def dos(rows) -> str:
    lines = [" # DOS at iteration step i = 12", " #  Energy[a.u.]       Density     Occupation"]
    return "\n".join(lines + ["".join(f"{value:15.8f}" for value in row) for row in rows]) + "\n"


METAL_DOS = dos([(-0.30, 0.25, 0.50), (-0.20, 0.25, 0.25), (-0.10, 0.50, 0.0)])
INSULATOR_DOS = dos([(-0.30, 0.50, 1.0), (-0.20, 0.0, 0.0), (0.10, 0.50, 0.0)])


def pair():
    from pymatgen.core import Lattice, Structure

    return Structure(Lattice.cubic(6.0), ["Ru", "O"], [[0, 0, 0], [0.3, 0, 0]])


def vasprun(symbols, forces, *, free_energy=-12.5, occupations=(1.0, 0.5, 0.0)) -> str:
    atoms = "".join(f"<rc><c>{symbol:<2}</c><c>   1</c></rc>" for symbol in symbols)
    rows = "".join(f"<v> {f[0]:.8f} {f[1]:.8f} {f[2]:.8f} </v>" for f in forces)
    bands = "".join(f"<r> {-1.0 + i:.4f} {value:.4f} </r>" for i, value in enumerate(occupations))
    step = (
        f"<calculation><scstep><energy><i name='e_fr_energy'> 99.0 </i></energy></scstep>"
        f"<energy><i name='e_fr_energy'> {free_energy:.6f} </i><i name='e_0_energy'> -12.4 </i></energy>"
        f"<varray name='forces'>{rows}</varray>"
        f"<eigenvalues><array><set><set comment='spin 1'><set comment='kpoint 1'>{bands}</set></set></set></array>"
        "</eigenvalues></calculation>"
    )
    return (
        "<?xml version='1.0' encoding='ISO-8859-1'?><modeling>"
        f"<atominfo><atoms>{len(symbols)}</atoms><array name='atoms'><set>{atoms}</set></array></atominfo>"
        "<calculation><energy><i name='e_fr_energy'> 0.0 </i></energy></calculation>"
        f"{step}</modeling>"
    )


def pw_output(forces) -> str:
    ry_per_bohr = 13.605693122994 / 0.529177210903
    lines = ["     Forces acting on atoms (cartesian axes, Ry/au):", ""]
    lines += [
        f"     atom {i + 1:4d} type  1   force = {f[0] / ry_per_bohr:14.8f}{f[1] / ry_per_bohr:14.8f}{f[2] / ry_per_bohr:14.8f}"
        for i, f in enumerate(forces)
    ]
    return "\n".join(lines) + "\n     Total force =     0.001\n"


QE_XML = """<?xml version="1.0"?>
<qes:espresso xmlns:qes="http://www.quantum-espresso.org/ns/qes/qes-1.0"><output><band_structure>
<ks_energies><eigenvalues size="3">-0.5 -0.2 0.1</eigenvalues><occupations size="3">1.0 0.4 0.0</occupations></ks_energies>
<ks_energies><eigenvalues size="3">-0.5 -0.2 0.1</eigenvalues><occupations size="3">1.0 1.0 0.0</occupations></ks_energies>
</band_structure></output></qes:espresso>"""


# ------------------------------------------------------------------ AiiDA-free layer


def test_rattled_copies_reproduce_the_validation_rattling() -> None:
    structures = {"a": pair(), "b": pair()}
    copies = psteros.rattled_copies(structures, count=2, sigma_angstrom=0.05, seed=7)
    assert list(copies) == ["a_r0", "a_r1", "b_r0", "b_r1"]
    rng = np.random.default_rng(7)  # the loop of validation/ruo2_mace_lovelace/c2_cp2k.py
    for label in ("a_r0", "a_r1", "b_r0", "b_r1"):
        expected = pair()
        for site in range(len(expected)):
            expected.translate_sites([site], rng.normal(0.0, 0.05, 3), frac_coords=False, to_unit_cell=True)
        assert np.allclose(copies[label].cart_coords, expected.cart_coords)
    assert np.allclose(structures["a"].cart_coords, pair().cart_coords)  # the originals are not moved
    with pytest.raises(ValueError, match="count"):
        psteros.rattled_copies(structures, count=0)
    with pytest.raises(ValueError, match="sigma"):
        psteros.rattled_copies(structures, sigma_angstrom=0.0)


def test_dft_labels_of_each_backend() -> None:
    structure = pair()
    cp2k = psteros.dft_label("cp2k", label="x", structure=structure,
                             parameters={"energy": -10.0, "energy_units": "a.u."},
                             files={"psteros_forces.dat": FORCES_AU})
    assert cp2k.energy_ev == pytest.approx(-10.0 * HARTREE_TO_EV)
    assert cp2k.forces_ev_per_angstrom[0][0] == pytest.approx(0.01 * HARTREE_TO_EV / BOHR_TO_ANGSTROM)

    forces = [[0.1, -0.2, 0.3], [-0.1, 0.2, -0.3]]
    qe = psteros.dft_label("qe", label="x", structure=structure, parameters={"energy": -50.0, "energy_units": "eV"},
                           files={"aiida.out": pw_output(forces)})
    assert qe.energy_ev == -50.0 and np.allclose(qe.forces_ev_per_angstrom, forces, atol=1e-6)

    vasp = psteros.dft_label("vasp", label="x", structure=structure, parameters={},
                             files={"vasprun.xml": vasprun(["Ru", "O"], forces)})
    assert vasp.energy_ev == pytest.approx(-12.5)  # the free energy of the last step, not sigma -> 0
    assert np.allclose(vasp.forces_ev_per_angstrom, forces)
    with pytest.raises(ValueError, match="order"):
        psteros.dft_label("vasp", label="x", structure=structure, parameters={},
                          files={"vasprun.xml": vasprun(["O", "Ru"], forces)})
    with pytest.raises(ValueError, match="no psteros_forces.dat"):
        psteros.dft_label("cp2k", label="x", structure=structure, parameters={"energy": -1.0}, files={})
    with pytest.raises(ValueError, match="backend"):
        psteros.dft_label("siesta", label="x", structure=structure, parameters={}, files={})


def test_diagram_energy_follows_the_psteros_readers() -> None:
    assert psteros.dft_diagram_energy("cp2k", label="x", parameters={"energy": -1.0, "energy_units": "a.u."}) == \
        pytest.approx(-HARTREE_TO_EV)
    assert psteros.dft_diagram_energy("qe", label="x", parameters={"energy": -3.0, "energy_units": "eV"}) == -3.0
    misc = {"total_energies": {"energy_extrapolated": -4.0, "energy_no_entropy": -4.1}}
    assert psteros.dft_diagram_energy("vasp", label="x", parameters=misc) == -4.0
    with pytest.raises(ValueError, match="Hartree"):
        psteros.dft_diagram_energy("cp2k", label="x", parameters={"energy": -1.0, "energy_units": "eV"})
    with pytest.raises(ValueError, match="no total energy"):
        psteros.dft_diagram_energy("vasp", label="x", parameters={})


def test_partly_filled_states_of_each_backend() -> None:
    assert psteros.dft_partially_occupied("cp2k", label="bulk", files={"psteros.dos": METAL_DOS}) == 1
    assert psteros.dft_partially_occupied("cp2k", label="bulk", files={"psteros.dos": INSULATOR_DOS}) == 0
    forces = [[0, 0, 0], [0, 0, 0]]
    metal = vasprun(["Ru", "O"], forces, occupations=(1.0, 0.5, 0.0))
    insulator = vasprun(["Ru", "O"], forces, occupations=(1.0, 1.0, 0.0))
    assert psteros.dft_partially_occupied("vasp", label="bulk", files={"vasprun.xml": metal}) == 1
    assert psteros.dft_partially_occupied("vasp", label="bulk", files={"vasprun.xml": insulator}) == 0
    assert mlip.qe_xml_occupations(QE_XML) == [1.0, 0.4, 0.0, 1.0, 1.0, 0.0]
    assert psteros.dft_partially_occupied("qe", label="bulk", files={"data-file-schema.xml": QE_XML}) == 1
    with pytest.raises(ValueError, match="no data-file-schema.xml"):
        psteros.dft_partially_occupied("qe", label="bulk", files={})


def test_job_config_runs_one_python_process() -> None:
    execution = psteros.ExecutionPolicy(queue="testegpu", with_mpi=False)
    assert psteros.MlipJobConfig("mace@gpu", execution).device == "cuda"
    with pytest.raises(ValueError, match="with_mpi=False"):
        psteros.MlipJobConfig("mace@gpu", psteros.ExecutionPolicy())
    with pytest.raises(ValueError, match="code_label"):
        psteros.MlipJobConfig(" ", execution)
    with pytest.raises(ValueError, match="device"):
        psteros.MlipJobConfig("mace@gpu", execution, device="tpu")


# ------------------------------------------------------------------ AiiDA layer


def _computer(tmp_path):
    from aiida import orm

    computer = orm.Computer(
        label=f"local-{tmp_path.name}", hostname="localhost", transport_type="core.local",
        scheduler_type="core.direct", workdir=str(tmp_path),
    ).store()
    orm.InstalledCode(computer=computer, filepath_executable="/bin/true", label="cp2k",
                      default_calc_job_plugin="cp2k").store()
    orm.InstalledCode(computer=computer, filepath_executable="/opt/mace/bin/python", label="mace",
                      default_calc_job_plugin="core.shell").store()
    return computer.label


def _replay_file(tmp_path):
    from ase import Atoms
    from ase.io import write

    frames = [Atoms("RuO2", positions=[[0, 0, 0], [1.9 + 0.01 * i, 0, 0], [-1.9, 0, 0]], cell=[8, 8, 8], pbc=True)
              for i in range(6)]
    frames.append(Atoms("Ag2", positions=[[0, 0, 0], [2.8, 0, 0]], cell=[8, 8, 8], pbc=True))
    path = tmp_path / "replay.xyz"
    write(path, frames, format="extxyz")
    return path


def _recipes(tmp_path, computer, run_type=None):
    from test_cp2k_backend import cp2k_config

    execution = psteros.ExecutionPolicy(computer=computer, queue="standard", max_concurrent_jobs=2)
    static = psteros.SurfaceWorkflowConfig(
        backend="cp2k", name="ruo2", execution=execution, calculation=cp2k_config(f"cp2k@{computer}", run_type),
        role_overrides={"o2": psteros.CalculationOverride(kpoints_distance=100.0)},
    )
    fine_tune = psteros.MaceFineTuneConfig(foundation_model="/models/mace-mh-1.model",
                                           replay_file=str(_replay_file(tmp_path)), epochs=2, seeds=(11, 12))
    job = psteros.MlipJobConfig(f"mace@{computer}", psteros.ExecutionPolicy(
        queue="testegpu", with_mpi=False, max_wallclock_seconds=1800,
        resources={"num_machines": 1, "num_mpiprocs_per_machine": 8},
        extra_options={"environment_variables": {"OMP_NUM_THREADS": "8"}}))
    return static, fine_tune, job


def _structures():
    bulk = psteros.rutile_bulk("Ru", a=4.49, c=3.11)
    slab, _ = psteros.rutile_110_slab(metal="Ru", termination="o_bridge", a=4.49, c=3.11, vacuum_angstrom=12.0)
    from pymatgen.core import Structure

    return {"bulk": bulk, "slab": Structure(slab.lattice, slab.species, slab.frac_coords),
            "o2": psteros.triplet_o2_cell(cell_length=10.0)}


def _outputs(graph) -> set:
    return set(graph.outputs._get_keys()) if hasattr(graph.outputs, "_get_keys") else set(graph.outputs)


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_mlip_workgraph_wires_the_validated_route(tmp_path) -> None:
    pytest.importorskip("aiida_cp2k")
    pytest.importorskip("aiida_shell")
    pytest.importorskip("ase")
    computer = _computer(tmp_path)
    static, fine_tune, job = _recipes(tmp_path, computer)
    structures = _structures()
    fixed = {"slab": psteros.central_sites(structures["slab"])}
    rattled = psteros.rattled_copies({"slab": structures["slab"]}, count=2)
    graph = psteros.build_mlip_surface_workgraph(
        structures, static, fine_tune, job=job, metallic_bulk="bulk", fixed_sites=fixed,
        training_structures=rattled)

    names = {task.name for task in graph.tasks}
    labels = list(structures)
    assert {f"{label}_static_cp2k" for label in [*labels, *rattled]} <= names
    assert {"mlip_round0_training", "mlip_round0_finetune", "mlip_round0_relax",
            "mlip_round1_training", "mlip_round1_finetune", "mlip_round1_relax"} <= names
    assert {f"{label}_al1_static_cp2k" for label in labels} <= names
    assert {f"{label}_mlip_dft_cp2k" for label in labels} <= names
    assert {f"{label}_mlip_dft_energy" for label in labels} <= names
    assert "mlip_round2_training" not in names and not any(name.endswith("_al2_static_cp2k") for name in names)
    assert not any(name.startswith(tuple(rattled)) and "mlip" in name for name in names)  # extra labels: statics only
    assert graph.max_number_jobs == 2

    links = {(link.from_task.name, link.to_task.name) for link in graph.links}
    assert ("slab_static_cp2k", "mlip_round0_training") in links
    assert ("slab_r1_static_cp2k", "mlip_round0_training") in links
    assert ("mlip_round0_training", "mlip_round0_finetune") in links
    assert ("mlip_round0_finetune", "mlip_round0_relax") in links
    assert ("mlip_round0_relax", "slab_mlip_round0_relaxed") in links
    assert ("slab_mlip_round0_relaxed", "slab_al1_static_cp2k") in links
    assert ("slab_al1_static_cp2k", "mlip_round1_training") in links
    assert ("slab_static_cp2k", "mlip_round1_training") in links  # the labels accumulate
    assert ("mlip_round1_relax", "slab_mlip_round1_relaxed") in links
    assert ("slab_mlip_round1_relaxed", "slab_mlip_dft_cp2k") in links
    assert ("slab_mlip_dft_cp2k", "slab_mlip_dft_energy") in links

    settings = graph.tasks["mlip_round1_training"].inputs.settings.value.get_dict()
    assert settings["frames"] == [*labels, *rattled, "bulk_al1", "slab_al1", "o2_al1"]
    assert settings["metallic_bulk"] == "bulk" and settings["backend"] == "cp2k"
    assert settings["fine_tune"]["seeds"] == [11, 12] and settings["name"] == "ruo2_round1"
    assert settings["replay_selection"]["kept"] == 6  # Ru-O structures only

    relax = graph.tasks["mlip_round0_relax"].inputs
    plan = json.loads(relax.nodes.job.value.get_object_content("job.json"))
    assert plan["models"] == {"seed11": {"path": "finetune/models/ruo2_round0_seed11.model", "head": "Default"},
                              "seed12": {"path": "finetune/models/ruo2_round0_seed12.model", "head": "Default"}}
    slab_item = next(item for item in plan["structures"] if item["label"] == "slab")
    assert slab_item["fixed"] == list(fixed["slab"])
    assert sorted(item["label"] for item in plan["structures"]) == sorted(labels)
    assert dict(relax.filenames.value) == {"finetune": "finetune"}
    assert relax.metadata.options.use_symlinks.value is True
    finetune = graph.tasks["mlip_round0_finetune"].inputs
    assert finetune.metadata.options.queue_name.value == "testegpu"
    assert finetune.metadata.options.withmpi.value is False
    assert finetune.metadata.options.environment_variables.value == {"OMP_NUM_THREADS": "8"}
    assert graph.tasks["slab_mlip_round1_relaxed"].inputs.settings.value.get_dict() == {"label": "slab", "model": "seed11"}

    outputs = _outputs(graph)
    for label in labels:
        assert {f"{label}_static_parameters", f"{label}_static_retrieved", f"{label}_mlip_relaxed_structure",
                f"{label}_mlip_energy", f"{label}_mlip_dft_energy", f"{label}_mlip_dft_parameters",
                f"{label}_mlip_dft_retrieved"} <= outputs
    assert {f"mlip_round{r}_{kind}" for r in (0, 1) for kind in ("finetune_summary", "model_folder", "relax_results")} \
        <= outputs

    plain = psteros.build_mlip_surface_workgraph(structures, static, fine_tune, job=job, metallic_bulk="bulk",
                                                 active_learning_rounds=0, dft_check=False)
    plain_names = {task.name for task in plain.tasks}
    assert "mlip_round1_training" not in plain_names and not any("mlip_dft" in name for name in plain_names)
    assert "slab_mlip_energy" in _outputs(plain) and "slab_mlip_dft_energy" not in _outputs(plain)


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_mlip_workgraph_refuses_inconsistent_inputs(tmp_path) -> None:
    pytest.importorskip("aiida_cp2k")
    pytest.importorskip("aiida_shell")
    computer = _computer(tmp_path)
    static, fine_tune, job = _recipes(tmp_path, computer)
    structures = _structures()
    build = psteros.build_mlip_surface_workgraph
    with pytest.raises(ValueError, match="metallic_bulk 'ruo2'"):
        build(structures, static, fine_tune, job=job, metallic_bulk="ruo2")
    with pytest.raises(ValueError, match=r"fixed_sites names labels that are not in structures: \['slab2'\]"):
        build(structures, static, fine_tune, job=job, metallic_bulk="bulk", fixed_sites={"slab2": (0,)})
    with pytest.raises(ValueError, match="slab: fixed_sites must be within"):
        build(structures, static, fine_tune, job=job, metallic_bulk="bulk", fixed_sites={"slab": (999,)})
    with pytest.raises(ValueError, match=r"repeat labels of structures: \['slab'\]"):
        build(structures, static, fine_tune, job=job, metallic_bulk="bulk", training_structures={"slab": structures["slab"]})
    with pytest.raises(ValueError, match=r"already used: \['slab_al1'\]"):
        build(structures, static, fine_tune, job=job, metallic_bulk="bulk",
              training_structures={"slab_al1": structures["slab"]})
    with pytest.raises(ValueError, match="non-negative"):
        build(structures, static, fine_tune, job=job, metallic_bulk="bulk", active_learning_rounds=-1)
    relax_static, _, _ = _recipes(tmp_path, computer, run_type="GEO_OPT")
    with pytest.raises(ValueError, match="bulk: the CP2K statics of the MLIP route need RUN_TYPE ENERGY_FORCE"):
        build(structures, relax_static, fine_tune, job=job, metallic_bulk="bulk")
    with pytest.raises(TypeError, match="MlipJobConfig"):
        build(structures, static, fine_tune, job="mace@gpu", metallic_bulk="bulk")


def _node_inputs(structure, energy_hartree, dos_text):
    from aiida import orm

    from psteros.backends.qe import as_aiida_structure

    folder = orm.FolderData()
    forces = "\n".join(
        [" ATOMIC FORCES in [eV/Angstrom]", "", " # Atom   Kind   Element          X      Y      Z"]
        + [f"      {i + 1}      1      {site.specie.symbol}      0.0100000000000000     0.0     -0.0100000000000000"
           for i, site in enumerate(structure)]
        + [" SUM OF ATOMIC FORCES  0.0 0.0 0.0 0.0"]) + "\n"
    import io

    folder.base.repository.put_object_from_filelike(io.BytesIO(forces.encode()), "psteros_forces.dat")
    folder.base.repository.put_object_from_filelike(io.BytesIO(dos_text.encode()), "psteros.dos")
    return as_aiida_structure(structure).store(), orm.Dict({"energy": energy_hartree, "energy_units": "a.u."}).store(), \
        folder.store()


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_training_folder_and_relaxed_structure_calcfunctions(tmp_path) -> None:
    pytest.importorskip("ase")
    from aiida import orm
    from ase.io import read

    from psteros.backends.mlip_tasks import mlip_relaxed, mlip_training_folder

    structures = {"bulk": psteros.rutile_bulk("Ru", a=4.49, c=3.11), "pair": pair()}
    inputs = {}
    for index, (label, structure) in enumerate(structures.items()):
        node, output, retrieved = _node_inputs(structure, -100.0 - index, METAL_DOS)
        inputs.update({f"structure_{label}": node, f"output_{label}": output, f"retrieved_{label}": retrieved})
    replay = read(_replay_file(tmp_path), index=":")
    for key, atoms in (("replay_train", replay[:4]), ("replay_valid", replay[4:6])):
        from ase.io import write

        write(tmp_path / f"{key}.extxyz", atoms, format="extxyz")
        inputs[key] = orm.SinglefileData(file=str(tmp_path / f"{key}.extxyz")).store()
    fine_tune = psteros.MaceFineTuneConfig(foundation_model="/models/mh1.model", replay_file="replay.xyz", seeds=(5,))
    from dataclasses import asdict

    settings = {"backend": "cp2k", "name": "t_round0", "device": "cuda", "metallic_bulk": "bulk",
                "frames": ["bulk", "pair"], "fine_tune": asdict(fine_tune), "replay_selection": {}}
    folder = mlip_training_folder._callable(settings=orm.Dict(settings), **inputs)
    names = set(folder.list_object_names())
    assert {"job.json", "train.extxyz", "finetune.py", "replay_train_geometries.extxyz",
            "replay_valid_geometries.extxyz"} <= names
    job = json.loads(folder.get_object_content("job.json"))
    assert [run["name"] for run in job["runs"]] == ["t_round0_seed5"]
    train = folder.get_object_content("train.extxyz")
    assert train.count("REF_energy=") == 2 and f"REF_energy={-100.0 * HARTREE_TO_EV:.10f}" in train

    insulating = dict(inputs)
    insulating["retrieved_bulk"] = _node_inputs(structures["bulk"], -100.0, INSULATOR_DOS)[2]
    with pytest.raises(ValueError, match="bulk: the bulk has no partly filled states"):
        mlip_training_folder._callable(settings=orm.Dict(settings), **insulating)

    # the relaxation job's results for one structure, with the model 'seed5'
    from ase import Atoms
    from ase.io import write

    start = inputs["structure_pair"]
    atoms = Atoms("RuO", positions=[[0.1, 0.0, 0.0], [1.9, 0.0, 0.0]], cell=np.array(start.cell), pbc=True)
    (tmp_path / "relaxed" / "seed5").mkdir(parents=True)
    write(tmp_path / "relaxed" / "seed5" / "pair.extxyz", atoms, format="extxyz")
    row = {"energy_ev": -7.5, "converged": True, "max_force_ev_per_angstrom": 0.004, "steps": 12}
    (tmp_path / "results.json").write_text(json.dumps({"models": {"seed5": {"structures": {"pair": row}}}}))
    results = orm.SinglefileData(file=str(tmp_path / "results.json")).store()
    relaxed = orm.FolderData(tree=str(tmp_path / "relaxed")).store()
    out = mlip_relaxed._callable(structure=start, results=results, relaxed=relaxed,
                       settings=orm.Dict({"label": "pair", "model": "seed5"}))
    assert out["energy"].value == -7.5
    assert [kind.name for kind in out["structure"].kinds] == [kind.name for kind in start.kinds]
    assert np.allclose([site.position for site in out["structure"].sites], atoms.positions)
    row["converged"] = False
    (tmp_path / "results.json").write_text(json.dumps({"models": {"seed5": {"structures": {"pair": row}}}}))
    unconverged = orm.SinglefileData(file=str(tmp_path / "results.json")).store()
    with pytest.raises(ValueError, match="pair: the relaxation with seed5 did not converge"):
        mlip_relaxed._callable(structure=start, results=unconverged, relaxed=relaxed,
                     settings=orm.Dict({"label": "pair", "model": "seed5"}))


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_existing_cp2k_builder_is_unchanged(tmp_path) -> None:
    """The MLIP route only adds: the CP2K static builder still makes one task and three outputs per label."""

    pytest.importorskip("aiida_cp2k")
    computer = _computer(tmp_path)
    static, _, _ = _recipes(tmp_path, computer)
    graph = psteros.build_cp2k_surface_workgraph(_structures(), static)
    assert {task.name for task in graph.tasks if not task.name.startswith("graph_")} == {"bulk_cp2k", "slab_cp2k", "o2_cp2k"}
    assert _outputs(graph) == {f"{label}_{kind}" for label in ("bulk", "slab", "o2")
                               for kind in ("parameters", "structure", "retrieved")}

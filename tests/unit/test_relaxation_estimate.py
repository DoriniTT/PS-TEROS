"""The second-order estimate of the DFT energy still to lose at a model geometry, and the error bars of gamma."""

from __future__ import annotations

import io
import json
from types import SimpleNamespace

import numpy as np
import pytest

import psteros
from psteros import relaxation_estimate


def spd(size, seed=1, shift=3.0):
    matrix = np.random.default_rng(seed).normal(size=(size, size))
    return matrix @ matrix.T + shift * np.eye(size)


def test_estimate_is_exact_for_a_quadratic_surface() -> None:
    hessian = spd(9)
    offset = np.random.default_rng(2).normal(scale=0.02, size=9)  # model geometry - DFT minimum
    forces = (-hessian @ offset).reshape(3, 3)
    estimate = relaxation_estimate(forces, hessian, sites=[0, 1, 2])
    assert estimate.energy_ev == pytest.approx(-0.5 * offset @ hessian @ offset, rel=1e-12)
    assert estimate.max_step_angstrom == pytest.approx(np.linalg.norm(offset.reshape(3, 3), axis=1).max())
    assert estimate.modes == 9 and estimate.negative_modes == estimate.soft_modes == 0 and estimate.reliable
    assert estimate.max_force_ev_per_angstrom == pytest.approx(np.linalg.norm(forces, axis=1).max())
    # the atoms that do not relax are left out: the forces of site 3 do not enter
    padded = np.vstack([forces, [[5.0, 5.0, 5.0]]])
    assert relaxation_estimate(padded, hessian, sites=[0, 1, 2]).energy_ev == pytest.approx(estimate.energy_ev)


def test_zero_soft_and_negative_modes_are_left_out() -> None:
    # two atoms on a spring along x: one zero mode (rigid translation along x) and stiff, soft and negative ones
    spring, soft, negative = 10.0, 0.01, -1.0
    hessian = np.zeros((6, 6))
    hessian[np.ix_([0, 3], [0, 3])] = spring * np.array([[1, -1], [-1, 1]])
    hessian[1, 1] = hessian[4, 4] = 5.0
    hessian[2, 2], hessian[5, 5] = soft, negative
    forces = np.array([[0.3, 0.0, 0.02], [-0.1, 0.0, 0.04]])  # net force 0.2 along x: on the translation
    estimate = relaxation_estimate(forces, hessian, sites=[0, 1], zero_modes=1)
    stretch = (forces[0, 0] - forces[1, 0]) / np.sqrt(2)
    assert estimate.energy_ev == pytest.approx(-0.5 * stretch ** 2 / (2 * spring))
    assert (estimate.modes, estimate.negative_modes, estimate.soft_modes) == (3, 1, 1)
    assert estimate.left_out_force_ev_per_angstrom == pytest.approx(np.hypot(0.02, 0.04))
    assert not estimate.reliable  # a negative mode: the model geometry is a saddle point
    assert "1 negative and 1 soft modes left out" in estimate.summary() and "NOT RELIABLE" in estimate.summary()
    calm = relaxation_estimate(forces * [1, 1, 0], np.diag([20, 5, 5, 20, 5, 5.0]), sites=[0, 1])
    assert calm.reliable and calm.summary().startswith("dE")
    far = relaxation_estimate(forces * [1, 1, 0], np.diag([0.1, 5, 5, 0.1, 5, 5.0]), sites=[0, 1])
    assert far.max_step_angstrom == pytest.approx(3.0) and not far.reliable  # a 3 A step: outside the quadratic model


def test_estimate_checks_its_inputs() -> None:
    with pytest.raises(ValueError, match="shape"):
        relaxation_estimate(np.zeros((2, 3)), np.eye(3), sites=[0, 1])
    with pytest.raises(ValueError, match="symmetric"):
        relaxation_estimate(np.zeros((1, 3)), np.triu(np.ones((3, 3))), sites=[0])
    with pytest.raises(ValueError, match="sites"):
        relaxation_estimate(np.zeros((1, 3)), np.eye(3), sites=[1])
    with pytest.raises(ValueError, match="zero_modes"):
        relaxation_estimate(np.zeros((1, 3)), np.eye(3), sites=[0], zero_modes=3)


def test_hessian_from_forces_is_the_matrix_of_the_modes() -> None:
    hessian = spd(6, seed=4)
    step = 0.01
    plus = np.array([(-hessian @ (np.eye(6)[i] * step)).reshape(2, 3) for i in range(6)])
    built = psteros.hessian_from_forces(displaced_sites=[0, 1], displacement_angstrom=step, forces_plus=plus,
                                        forces_minus=-plus)
    assert np.allclose(built, hessian)
    modes = psteros.harmonic_vibrations_from_forces(masses_amu=[1.0, 1.0], displaced_sites=[0, 1],
                                                    displacement_angstrom=step, forces_plus=plus, forces_minus=-plus)
    expected = np.sqrt(np.linalg.eigvalsh(hessian)) * 521.47  # cm-1 for unit masses (eV/A2/amu)
    assert np.allclose(sorted(modes.frequencies_cm1), sorted(expected), rtol=1e-3)
    with pytest.raises(ValueError, match="shape"):
        psteros.hessian_from_forces(displaced_sites=[0], displacement_angstrom=step, forces_plus=plus, forces_minus=plus)


def _result(curves):
    points = {label: tuple(SimpleNamespace(gamma_ev_per_angstrom2=value) for value in values)
              for label, values in curves.items()}
    return SimpleNamespace(diagram=SimpleNamespace(curves=points))


def test_gamma_error_bars_from_two_results() -> None:
    good = relaxation_estimate(np.zeros((1, 3)), np.eye(3) * 5, sites=[0])
    bad = psteros.RelaxationEstimate(-0.01, 0.2, 0.3, 2, 1, 0, 0.1, False)
    estimates = {"bulk": good, "slab_a": good, "slab_b": bad}
    base = _result({"slab_a": [0.05, 0.06], "slab_b": [0.07, 0.07]})
    lowered = _result({"slab_a": [0.049, 0.0595], "slab_b": [0.068, 0.068]})
    bars = psteros.gamma_error_bars(base, lowered, estimates, ["slab_a", "slab_b"])
    assert bars["slab_a"].delta_gamma_ev_per_angstrom2 == pytest.approx(-0.001)  # the largest change on the curve
    assert bars["slab_a"].reliable and not bars["slab_b"].reliable
    assert bars["slab_b"].summary().startswith("slab_b: d gamma -2.000 meV/A2") and "[not reliable]" in bars["slab_b"].summary()
    unreliable_reference = psteros.gamma_error_bars(base, lowered, {**estimates, "bulk": bad}, ["slab_a"])
    assert not unreliable_reference["slab_a"].reliable  # a reference that is far from its minimum spoils every gamma


def test_read_relaxation_estimates_from_a_finished_graph(monkeypatch) -> None:
    """The reader with stand-ins for the graph's processes: DFT forces of the final statics, the Hessian job."""

    from psteros import mlip_workflow

    pymatgen = pytest.importorskip("pymatgen.core")
    slab = pymatgen.Structure(pymatgen.Lattice.tetragonal(4, 20), ["Mg", "O", "Mg"], [[0, 0, 0], [0.5, 0, 0.1], [0, 0, 0.2]])
    o2 = pymatgen.Structure(pymatgen.Lattice.cubic(10), ["O", "O"], [[0, 0, 0], [0.12, 0, 0]])
    structures = {"slab": slab, "ref_O": o2}
    sites = {"slab": [0, 2], "ref_O": [0, 1]}  # site 1 of the slab is fixed
    bond = np.zeros((6, 6))
    bond[np.ix_([0, 3], [0, 3])] = 30.0 * np.array([[1, -1], [-1, 1]])  # O2 along x: one stretch, five zero modes
    stiffness = {"slab": spd(6, seed=5), "ref_O": bond}
    offsets = {"slab": np.array([0.01, 0, 0, 0, 0, -0.01]), "ref_O": np.array([0.005, 0, 0, -0.005, 0, 0])}

    def npz(label):
        rows = []
        for index in range(len(sites[label]) * 3):
            full = np.zeros((len(structures[label]), 3))
            column = -stiffness[label][:, index] * 0.01
            full[sites[label]] = column.reshape(len(sites[label]), 3)
            rows.append(full)
        buffer = io.BytesIO()
        np.savez(buffer, forces_plus=np.array(rows), forces_minus=-np.array(rows))
        return buffer.getvalue()

    def dft_forces(label):
        full = np.zeros((len(structures[label]), 3))
        full[sites[label]] = (-stiffness[label] @ offsets[label]).reshape(len(sites[label]), 3)
        return full

    folder = SimpleNamespace(base=SimpleNamespace(repository=SimpleNamespace(
        get_object_content=lambda name, mode="r": npz(name[:-4]))))
    nodes = {
        "mlip_round0_training": SimpleNamespace(inputs=SimpleNamespace(settings=SimpleNamespace(
            get_dict=lambda: {"backend": "vasp"}))),
        "mlip_hessian_job": SimpleNamespace(inputs=SimpleNamespace(settings=SimpleNamespace(
            get_dict=lambda: {"displaced_sites": sites, "displacement_angstrom": 0.01}))),
        "mlip_hessian": SimpleNamespace(is_finished_ok=True, outputs=SimpleNamespace(vibrations=folder)),
        **{f"{label}_mlip_dft_vasp": SimpleNamespace(is_finished_ok=True, outputs=SimpleNamespace(
            misc=SimpleNamespace(get_dict=lambda label=label: {"label": label}), retrieved=label))
           for label in structures},
    }
    relaxed = {label: SimpleNamespace(get_pymatgen_structure=lambda s=s: s, sites=list(s)) for label, s in structures.items()}
    monkeypatch.setattr(mlip_workflow, "_task_node", lambda graph, name: nodes.get(name))
    monkeypatch.setattr(mlip_workflow, "_relaxed_structures", lambda graph, labels: relaxed)
    monkeypatch.setattr("psteros.backends.mlip_tasks.retrieved_files", lambda folder, backend: {})
    monkeypatch.setattr("psteros.backends.vibration_tasks._is_linear", lambda structure: True)
    monkeypatch.setattr("psteros.mlip.dft_label", lambda backend, *, label, structure, parameters, files:
                        SimpleNamespace(forces_ev_per_angstrom=dft_forces(parameters["label"])))

    estimates = psteros.read_relaxation_estimates(object(), ["slab", "ref_O"], molecules=["ref_O"])
    for label in structures:
        expected = -0.5 * offsets[label] @ stiffness[label] @ offsets[label]
        assert estimates[label].energy_ev == pytest.approx(expected, rel=1e-6)
    assert estimates["slab"].modes == 6  # a fixed centre: no zero mode
    assert estimates["ref_O"].modes == 1  # a linear molecule: 5 zero modes out of 6
    nodes.pop("mlip_hessian_job")
    with pytest.raises(ValueError, match="build it with relaxation_estimate=True"):
        psteros.read_relaxation_estimates(object(), ["slab"])


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_mlip_graph_with_the_hessian_job(tmp_path) -> None:
    pytest.importorskip("aiida_cp2k")
    pytest.importorskip("aiida_shell")
    pytest.importorskip("ase")
    from test_mlip_workflow import _computer, _recipes, _structures

    computer = _computer(tmp_path)
    static, fine_tune, job = _recipes(tmp_path, computer)
    structures = _structures()
    fixed = {"slab": psteros.central_sites(structures["slab"])}
    graph = psteros.build_mlip_surface_workgraph(structures, static, fine_tune, job=job, metallic_bulk="bulk",
                                                 fixed_sites=fixed, relaxation_estimate=True)
    names = {task.name for task in graph.tasks}
    assert {"mlip_hessian_job", "mlip_hessian"} <= names
    settings = graph.tasks["mlip_hessian_job"].inputs.settings.value.get_dict()
    free = [index for index in range(len(structures["slab"])) if index not in set(fixed["slab"])]
    assert settings["displaced_sites"]["slab"] == free and settings["displaced_sites"]["o2"] == [0, 1]
    assert settings["model"] == ["finetune/models/ruo2_round1_seed11.model", "Default"]
    assert settings["displacement_angstrom"] == psteros.mlip_workflow.HESSIAN_DISPLACEMENT_ANGSTROM
    links = {(link.from_task.name, link.to_task.name) for link in graph.links}
    assert ("slab_mlip_round1_relaxed", "mlip_hessian_job") in links and ("mlip_round1_finetune", "mlip_hessian") in links

    plain = psteros.build_mlip_surface_workgraph(structures, static, fine_tune, job=job, metallic_bulk="bulk")
    assert not {"mlip_hessian_job", "mlip_hessian"} & {task.name for task in plain.tasks}  # left out: as before
    with pytest.raises(ValueError, match="relaxation_estimate needs dft_check"):
        psteros.build_mlip_surface_workgraph(structures, static, fine_tune, job=job, metallic_bulk="bulk",
                                             dft_check=False, relaxation_estimate=True)
    assert json.loads(graph.tasks["mlip_round0_relax"].inputs.nodes.job.value.get_object_content("job.json"))

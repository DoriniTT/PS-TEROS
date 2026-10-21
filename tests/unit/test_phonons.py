"""Phonons on a q-point mesh with phonopy (issue #3): the pure layer, the DFT graph and the MLIP route.

ASE's EMT potential stands in for the forces of DFT or of a model.
"""

from __future__ import annotations

import json
import warnings

import numpy as np
import pytest

import psteros
from psteros.phonons import THZ_TO_CM1, vacuum_along_c_angstrom

pytest.importorskip("pymatgen.core")
pytest.importorskip("ase")


def _pymatgen(atoms):
    from pymatgen.io.ase import AseAtomsAdaptor

    return AseAtomsAdaptor.get_structure(atoms)


def cu_bulk():
    from ase.build import bulk

    return _pymatgen(bulk("Cu", "fcc", a=3.59))


def cu3au_bulk():
    """L1_2 Cu3Au, four atoms in the cubic cell."""

    from ase.build import bulk

    atoms = bulk("Cu", "fcc", a=3.75, cubic=True)
    atoms.symbols[0] = "Au"
    return _pymatgen(atoms)


def cu_slab(layers=3):
    from ase.build import fcc100
    from ase.calculators.emt import EMT
    from ase.optimize import BFGS

    atoms = fcc100("Cu", size=(1, 1, layers), a=3.59, vacuum=7.0, periodic=True)
    atoms.calc = EMT()
    BFGS(atoms, logfile=None).run(fmax=1e-4)
    return _pymatgen(atoms)


def emt_forces(structures):
    from ase.calculators.emt import EMT
    from pymatgen.io.ase import AseAtomsAdaptor

    result = []
    for structure in structures:
        atoms = AseAtomsAdaptor.get_atoms(structure)
        atoms.calc = EMT()
        result.append(atoms.get_forces())
    return np.array(result)


def emt_modes(structure, supercell, mesh, **options):
    pytest.importorskip("phonopy")
    plan = psteros.phonopy_displacements(structure, supercell)
    perfect, displaced = plan.structures(structure)
    return psteros.phonopy_modes(structure, plan, emt_forces(displaced),
                                 reference_forces_ev_per_angstrom=emt_forces([perfect])[0], mesh=mesh, **options)


# ------------------------------------------------------------------ the default is unchanged


def test_default_vibrations_config_keeps_the_gamma_point_modes() -> None:
    """Leaving phonopy out is the config of before: no phonopy, the same fields and equality."""

    config = psteros.VibrationsConfig(displacement_angstrom=0.015, supercells={"bulk": (2, 2, 2)}, molecules=("o2",))
    assert config.phonopy is None
    assert config == psteros.VibrationsConfig(displacement_angstrom=0.015, supercells={"bulk": (2, 2, 2)},
                                              molecules=("o2",), phonopy=None)
    assert list(psteros.VibrationsConfig.__dataclass_fields__)[:4] == [
        "displacement_angstrom", "displaced_sites", "supercells", "molecules"]
    assert config.labels == {"bulk", "o2"}


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_gamma_point_graphs_are_those_of_before(tmp_path) -> None:
    """Without phonopy both builders give the Gamma-point graphs of before: the same tasks, job script and settings."""

    pytest.importorskip("aiida_cp2k")
    pytest.importorskip("aiida_shell")
    pytest.importorskip("aiida_vasp")
    from test_mlip_workflow import _computer, _recipes, _structures
    from test_vibrations import _computer as _vasp_computer

    static, fine_tune, job = _recipes(tmp_path, _computer(tmp_path))
    structures = _structures()
    graph = psteros.build_mlip_surface_workgraph(
        structures, static, fine_tune, job=job, metallic_bulk="bulk",
        vibrations=psteros.VibrationsConfig(supercells={"bulk": (2, 1, 1)}, molecules=("o2",)))
    assert list(graph.tasks["mlip_vibrations"].inputs.arguments.value) == ["-u", "vibrations.py"]
    assert set(graph.tasks["mlip_vibrations_job"].inputs.settings.value.get_dict()) == {
        "model", "displaced_sites", "device", "displacement_angstrom", "default_dtype"}
    assert set(graph.tasks["slab_mlip_vibrations"].inputs.settings.value.get_dict()) == {
        "label", "model", "displaced_sites", "displacement_angstrom", "supercell_size", "molecule"}

    computer = _vasp_computer(tmp_path / "vasp", "vasp", "vasp.vasp")
    vasp = psteros.SurfaceWorkflowConfig(
        backend="vasp", execution=psteros.ExecutionPolicy(computer=computer),
        calculation=psteros.VaspCalculationConfig(f"vasp@{computer}", {"ENCUT": 400, "NSW": 0, "IBRION": -1}))
    dft = psteros.build_vibrations_workgraph({"bulk": structures["bulk"], "o2": structures["o2"]}, vasp,
                                             psteros.VibrationsConfig(molecules=("o2",)))
    names = {task.name for task in dft.tasks}
    assert not any("_vib_p" in name or "perfect" in name for name in names)
    assert sum(name.startswith("bulk_vib_s") and name.endswith("_vasp") for name in names) == 6 * len(structures["bulk"])
    assert set(dft.tasks["bulk_vibrations"].inputs.settings.value.get_dict()) == {
        "backend", "displaced_sites", "displacement_angstrom", "supercell_size", "molecule"}


# ------------------------------------------------------------------ the config


def test_phonopy_config_chooses_supercells_and_meshes_by_length() -> None:
    bulk, slab = cu3au_bulk(), cu_slab()
    config = psteros.PhonopyConfig(supercell_length_angstrom=10.0, mesh_length_angstrom=40.0, slabs=("slab",))
    assert config.supercell("bulk", bulk) == (3, 3, 3)  # ceil(10 / 3.75)
    assert config.mesh("bulk", bulk) == (11, 11, 11)  # round(40 / 3.75)
    width = 3.59 / np.sqrt(2)
    assert config.supercell("slab", slab) == (int(np.ceil(10 / width)),) * 2 + (1,)
    assert config.mesh("slab", slab) == (int(round(40 / width)),) * 2 + (1,)
    explicit = psteros.PhonopyConfig(supercells={"bulk": (2, 2, 2), "slab": (4, 4, 1)}, meshes={"slab": (8, 8, 1)},
                                     slabs=("slab",))
    assert explicit.supercell("bulk", bulk) == (2, 2, 2) and explicit.mesh("slab", slab) == (8, 8, 1)
    assert explicit.labels == {"bulk", "slab"}

    assert vacuum_along_c_angstrom(slab) > 10 and vacuum_along_c_angstrom(bulk) < 3
    config.check("slab", slab)
    with pytest.raises(ValueError, match="bulk: PhonopyConfig.slabs lists it, but it has no vacuum"):
        psteros.PhonopyConfig(slabs=("bulk",)).check("bulk", bulk)
    with pytest.warns(UserWarning, match="other: a gap of .* PhonopyConfig.slabs does not list it"):
        config.check("other", slab)

    with pytest.raises(ValueError, match="a slab has 1 along c"):
        psteros.PhonopyConfig(slabs=("slab",), meshes={"slab": (8, 8, 2)})
    with pytest.raises(ValueError, match="three positive integers"):
        psteros.PhonopyConfig(supercells={"bulk": (2, 2)})
    with pytest.raises(ValueError, match="supercell_length_angstrom must be a positive number"):
        psteros.PhonopyConfig(supercell_length_angstrom=0.0)
    with pytest.raises(ValueError, match="distinct labels"):
        psteros.PhonopyConfig(slabs=("a", "a"))


def test_vibrations_config_with_phonopy_is_strict() -> None:
    phonopy = psteros.PhonopyConfig(slabs=("slab",))
    config = psteros.VibrationsConfig(molecules=("o2",), phonopy=phonopy)
    assert config.labels == {"o2", "slab"}
    with pytest.raises(ValueError, match=r"leave displaced_sites and supercells out \(given for bulk, slab\)"):
        psteros.VibrationsConfig(supercells={"bulk": (2, 2, 2)}, displaced_sites={"slab": (0,)}, phonopy=phonopy)
    with pytest.raises(ValueError, match="o2 are molecules"):
        psteros.VibrationsConfig(molecules=("o2",), phonopy=psteros.PhonopyConfig(slabs=("o2",)))
    with pytest.raises(TypeError, match="must be a PhonopyConfig"):
        psteros.VibrationsConfig(phonopy={"slabs": ("slab",)})


# ------------------------------------------------------------------ the physics


def test_commensurate_mesh_gives_the_gamma_point_modes_of_the_supercell() -> None:
    """On the q-points commensurate with the supercell, Fourier interpolation is exact: phonopy's modes are those of
    the Gamma-point finite differences in that supercell (the existing route), per input cell."""

    from ase.calculators.emt import EMT
    from pymatgen.io.ase import AseAtomsAdaptor

    bulk = cu3au_bulk()
    modes = emt_modes(bulk, (2, 2, 2), (2, 2, 2), shifted_mesh=False)
    atoms = AseAtomsAdaptor.get_atoms(bulk).repeat((2, 2, 2))
    atoms.calc = EMT()
    reference = atoms.get_positions()
    plus, minus = [], []
    for site in range(len(atoms)):
        for axis in range(3):
            for sign, store in ((1.0, plus), (-1.0, minus)):
                positions = reference.copy()
                positions[site, axis] += sign * 0.01
                atoms.set_positions(positions)
                store.append(atoms.get_forces())
    gamma = psteros.harmonic_vibrations_from_forces(
        masses_amu=atoms.get_masses(), displaced_sites=list(range(len(atoms))), displacement_angstrom=0.01,
        forces_plus=plus, forces_minus=minus, composition={"Cu": 3, "Au": 1}, zero_modes=3, supercell_size=8)
    assert modes.zero_modes == 3 and modes.vibrations.supercell_size == 8
    assert len(modes.vibrations.frequencies_cm1) == len(gamma.frequencies_cm1) == 8 * 12 - 3
    assert np.allclose(np.sort(modes.vibrations.frequencies_cm1), np.sort(gamma.frequencies_cm1), atol=0.1)
    for temperature in (0.0, 300.0, 900.0):  # phonopy symmetrises the force constants and adds -d where needed only
        assert modes.vibrations.free_energy_ev(temperature) == pytest.approx(gamma.free_energy_ev(temperature), abs=1e-4)
    assert modes.vibrations.composition == {"Cu": 3, "Au": 1}


def test_free_energy_is_phonopys_on_the_same_mesh() -> None:
    """The mesh average of :class:`HarmonicVibrations` is phonopy's own thermal free energy (THz -> cm^-1, kJ/mol)."""

    phonopy = pytest.importorskip("phonopy")
    bulk = cu_bulk()
    plan = psteros.phonopy_displacements(bulk, (4, 4, 4))
    perfect, displaced = plan.structures(bulk)
    forces = emt_forces(displaced) - emt_forces([perfect])[0]
    modes = psteros.phonopy_modes(bulk, plan, forces, reference_forces_ev_per_angstrom=None, mesh=(8, 8, 8))
    unit = phonopy.structure.atoms.PhonopyAtoms(symbols=["Cu"], cell=bulk.lattice.matrix,
                                                scaled_positions=bulk.frac_coords, masses=[63.546])
    phonon = phonopy.Phonopy(unit, supercell_matrix=np.diag([4, 4, 4]), primitive_matrix=np.eye(3))
    phonon.dataset = {"natom": 64, "first_atoms": [{"number": a, "displacement": np.array(v), "forces": f}
                                                   for (a, v), f in zip(plan.displacements, forces)]}
    phonon.produce_force_constants()
    phonon.symmetrize_force_constants()
    phonon.run_mesh([8, 8, 8], is_gamma_center=False)  # an even mesh: shifted by half a step, as psteros's
    phonon.run_thermal_properties(temperatures=[300.0, 900.0])
    expected = np.asarray(phonon.thermal_properties.free_energy) / 96.485_332_12  # kJ/mol -> eV
    assert modes.zero_modes == 0 and modes.vibrations.supercell_size == 512
    assert [modes.vibrations.free_energy_ev(t) for t in (300.0, 900.0)] == pytest.approx(expected, abs=1e-6)
    assert THZ_TO_CM1 == pytest.approx(33.356_41, rel=1e-6)


def test_forces_of_the_undisplaced_supercell_are_subtracted() -> None:
    """A geometry away from the minimum (DFT at a model's geometry): its residual forces leave the modes unchanged."""

    pytest.importorskip("phonopy")
    bulk = cu3au_bulk()
    plan = psteros.phonopy_displacements(bulk, (2, 2, 2))
    perfect, displaced = plan.structures(bulk)
    forces, reference = emt_forces(displaced), emt_forces([perfect])[0]
    residual = np.random.default_rng(3).normal(scale=0.05, size=reference.shape)
    clean = psteros.phonopy_modes(bulk, plan, forces, reference_forces_ev_per_angstrom=reference, mesh=(4, 4, 4))
    offset = psteros.phonopy_modes(bulk, plan, forces + residual, reference_forces_ev_per_angstrom=reference + residual,
                                   mesh=(4, 4, 4))
    assert np.allclose(offset.vibrations.frequencies_cm1, clean.vibrations.frequencies_cm1, atol=1e-6)
    raw = psteros.phonopy_modes(bulk, plan, forces + residual, reference_forces_ev_per_angstrom=None, mesh=(4, 4, 4),
                                imaginary_modes="drop")
    assert not np.allclose(raw.vibrations.frequencies_cm1, clean.vibrations.frequencies_cm1, atol=1.0)
    with pytest.raises(ValueError, match="forces_ev_per_angstrom must have shape"):
        psteros.phonopy_modes(bulk, plan, forces[:-1], reference_forces_ev_per_angstrom=reference, mesh=(4, 4, 4))
    with pytest.raises(ValueError, match="not the one of the displacements"):
        psteros.phonopy_modes(cu_bulk(), plan, forces, reference_forces_ev_per_angstrom=reference, mesh=(4, 4, 4))


def test_stability_check_reports_imaginary_modes() -> None:
    pytest.importorskip("phonopy")
    bulk = cu_bulk()
    stable = emt_modes(bulk, (3, 3, 3), (6, 6, 6))
    # the shifted 6x6x6 mesh and the 27 commensurate points, Gamma left out (only acoustic modes in a one-atom cell)
    assert stable.stability.stable and stable.stability.q_points == 216 + 26
    assert "stable: no imaginary mode at 242 q-points" in stable.stability.summary()
    plan = psteros.phonopy_displacements(bulk, (3, 3, 3))
    perfect, displaced = plan.structures(bulk)
    reference = emt_forces([perfect])[0]
    flipped = 2 * reference - emt_forces(displaced)  # the opposite force constants: a maximum, not a minimum
    with pytest.raises(ValueError, match="imaginary mode"):
        psteros.phonopy_modes(bulk, plan, flipped, reference_forces_ev_per_angstrom=reference, mesh=(6, 6, 6))
    unstable = psteros.phonopy_modes(bulk, plan, flipped, reference_forces_ev_per_angstrom=reference, mesh=(6, 6, 6),
                                     imaginary_modes="drop")
    assert not unstable.stability.stable and unstable.stability.lowest_frequency_cm1 < -100
    assert len(unstable.stability.imaginary) == unstable.stability.q_points
    assert unstable.stability.summary().startswith("NOT STABLE: imaginary modes at 242 of 242 q-points")
    again = psteros.PhononStability.from_dict(unstable.to_dict())
    assert again == unstable.stability


def test_slab_vibrates_as_a_whole() -> None:
    """Full-slab phonons: in-plane supercell, 2D mesh, every atom displaced, no region counted as bulk."""

    pytest.importorskip("phonopy")
    slab, bulk = cu_slab(3), cu_bulk()
    config = psteros.PhonopyConfig(slabs=("slab",), supercells={"bulk": (4, 4, 4), "slab": (4, 4, 1)},
                                   meshes={"bulk": (12, 12, 12), "slab": (12, 12, 1)})
    slab_modes = emt_modes(slab, config.supercell("slab", slab), config.mesh("slab", slab))
    bulk_modes = emt_modes(bulk, config.supercell("bulk", bulk), config.mesh("bulk", bulk))
    vibrations = slab_modes.vibrations
    assert vibrations.composition == {"Cu": 3} and vibrations.frozen_composition == {}
    assert vibrations.supercell_size == 144 and len(vibrations.frequencies_cm1) == 144 * 9
    assert slab_modes.stability.stable
    area = slab.lattice.a * slab.lattice.b
    gamma_vib = (vibrations.free_energy_ev(900.0) - 3 * bulk_modes.vibrations.free_energy_ev(900.0)) / (2 * area)
    assert -0.02 < gamma_vib < -0.005  # eV/A^2: the softer surface lowers gamma
    # the frozen-centre bookkeeping is not needed: E + F_vib(slab)
    free = psteros.free_energies_with_vibrations({"slab": -10.0, "bulk": -3.5},
                                                 {"slab": vibrations, "bulk": bulk_modes.vibrations}, 900.0,
                                                 bulk_label="bulk")
    assert free["slab"] == pytest.approx(-10.0 + vibrations.free_energy_ev(900.0))


def test_displacements_keep_kinds_apart_and_round_trip() -> None:
    pytest.importorskip("phonopy")
    from pymatgen.core import Lattice, Structure

    plain = Structure(Lattice.cubic(3.0), ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    magnetic = plain.copy(site_properties={"kind_name": ["Fe1", "Fe2"]})
    assert len(psteros.phonopy_displacements(plain, (2, 2, 2))) == 1  # bcc: one displacement
    plan = psteros.phonopy_displacements(magnetic, (2, 2, 2), plus_minus=True)
    assert len(plan) == 4  # two inequivalent kinds, +d and -d
    assert sorted(plan.unit_sites) == [0] * 8 + [1] * 8
    assert psteros.PhonopyDisplacements.from_dict(json.loads(json.dumps(plan.to_dict()))) == plan
    perfect, displaced = plan.structures(magnetic)
    assert len(perfect) == 16 and perfect.site_properties["kind_name"].count("Fe2") == 8
    atom, vector = plan.displacements[0]
    moved = displaced[0].cart_coords - perfect.cart_coords
    assert np.allclose(moved[atom], vector) and np.allclose(np.delete(moved, atom, axis=0), 0.0)
    assert np.linalg.norm(vector) == pytest.approx(0.01)


# ------------------------------------------------------------------ the GPU job


def test_forces_job_gives_the_forces_of_the_displacements(tmp_path) -> None:
    import subprocess
    import sys

    pytest.importorskip("phonopy")
    from test_mlip_route import _stub_mace

    bulk = cu3au_bulk()
    plan = psteros.phonopy_displacements(bulk, (2, 1, 1))
    perfect, displaced = plan.structures(bulk)
    folder = psteros.prepare_forces_job(tmp_path / "job", {"cu3au": perfect}, ("finetune/models/x_seed1.model", "Default"),
                                        {"cu3au": plan.displacements}, device="cpu")
    job = json.loads((folder / "job.json").read_text())
    assert job["structures"][0]["displacements"] == [[atom, *vector] for atom, vector in plan.displacements]
    environment = {"PYTHONPATH": _stub_mace(tmp_path), "PATH": "/usr/bin:/bin"}
    completed = subprocess.run([sys.executable, "-u", "forces.py"], cwd=folder, env=environment, capture_output=True,
                               text=True, timeout=300)
    assert completed.returncode == 0, completed.stderr
    results = json.loads((folder / "results.json").read_text())
    assert results["status"] == "completed" and results["structures"]["cu3au"]["force_evaluations"] == len(plan) + 1
    data = np.load(folder / "forces" / "cu3au.npz")
    assert np.allclose(data["forces"], emt_forces(displaced), atol=1e-6)
    assert np.allclose(data["reference_forces"], emt_forces([perfect])[0], atol=1e-6)
    assert np.allclose(data["positions"], plan.positions_angstrom, atol=1e-8)

    with pytest.raises(ValueError, match="every structure"):
        psteros.prepare_forces_job(tmp_path / "x", {"cu3au": perfect}, ("m", "Default"), {})
    with pytest.raises(ValueError, match="cu3au: each displacement must be"):
        psteros.prepare_forces_job(tmp_path / "y", {"cu3au": perfect}, ("m", "Default"), {"cu3au": [(99, (0.01, 0, 0))]})
    with pytest.raises(ValueError, match="at most 0.1 A"):
        psteros.prepare_forces_job(tmp_path / "z", {"cu3au": perfect}, ("m", "Default"), {"cu3au": [(0, (0.5, 0, 0))]})


# ------------------------------------------------------------------ AiiDA: the DFT graph


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_dft_phonopy_graph_runs_the_symmetry_reduced_displacements(tmp_path) -> None:
    pytest.importorskip("aiida_vasp")
    pytest.importorskip("phonopy")
    from test_vibrations import _computer

    computer = _computer(tmp_path, "vasp", "vasp.vasp")
    incar = {"ENCUT": 520, "EDIFF": 1e-8, "NSW": 0, "IBRION": -1}
    slab, _ = psteros.sno2_110_slab(termination="o", triple_layers=3)
    fixed = tuple(psteros.central_sites(slab, half_width=1.5))
    static = psteros.SurfaceWorkflowConfig(
        backend="vasp", name="sno2", role_overrides={"slab": psteros.CalculationOverride(fixed_sites=fixed)},
        execution=psteros.ExecutionPolicy(computer=computer, max_concurrent_jobs=4),
        calculation=psteros.VaspCalculationConfig(f"vasp@{computer}", incar, potential_mapping={"Sn": "Sn_d"}),
    )
    bulk, o2 = psteros.rutile_sno2_bulk(), psteros.triplet_o2_cell()
    config = psteros.PhonopyConfig(supercell_length_angstrom=6.0, slabs=("slab",))
    graph = psteros.build_vibrations_workgraph(
        {"slab": slab, "bulk": bulk, "o2": o2}, static, psteros.VibrationsConfig(molecules=("o2",), phonopy=config))
    names = {task.name for task in graph.tasks}
    for label, structure in (("bulk", bulk), ("slab", slab)):
        settings = graph.tasks[f"{label}_vibrations"].inputs.settings.value.get_dict()
        plan = psteros.PhonopyDisplacements.from_dict(settings["displacements"])
        assert settings["method"] == "phonopy" and settings["backend"] == "vasp"
        assert settings["supercell"] == list(config.supercell(label, structure)) == list(plan.supercell)
        assert settings["mesh"] == list(config.mesh(label, structure)) and settings["shifted_mesh"] is True
        calculations = {name for name in names if name.startswith(f"{label}_vib_")}
        assert calculations == {f"{label}_vib_{key}_vasp" for key in ["perfect", *(f"p{i}" for i in range(len(plan)))]}
        cells = len(structure) * int(np.prod(plan.supercell))
        assert len(graph.tasks[f"{label}_vib_p0_vasp"].inputs.structure.value.sites) == cells
    assert graph.tasks["slab_vibrations"].inputs.settings.value.get_dict()["supercell"][2] == 1
    # the full slab is displaced (the frozen centre constrains only the relaxation), symmetry-reduced
    slab_plan = psteros.PhonopyDisplacements.from_dict(
        graph.tasks["slab_vibrations"].inputs.settings.value.get_dict()["displacements"])
    assert {slab_plan.unit_sites[atom] for atom, _ in slab_plan.displacements} - set(fixed)
    assert len(slab_plan) < 6 * len(slab)
    # the molecule keeps its Gamma-point finite differences
    assert sum(name.startswith("o2_vib_s") for name in names) == 12
    assert graph.tasks["o2_vibrations"].inputs.settings.value.get_dict()["molecule"] is True
    vasp = graph.tasks["slab_vib_p0_vasp"]
    assert vasp.inputs.parameters.value.get_dict() == {"incar": {key.lower(): value for key, value in incar.items()}}


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_dft_phonopy_modes_calcfunction() -> None:
    pytest.importorskip("phonopy")
    from test_vibrations import _folder, vasprun

    from psteros.backends.qe import as_aiida_structure
    from psteros.backends.vibration_tasks import phonopy_harmonic_modes
    from psteros.vibrations_workflow import phonopy_plan, phonopy_settings

    bulk = cu3au_bulk()
    node = as_aiida_structure(bulk)
    vibrations = psteros.VibrationsConfig(phonopy=psteros.PhonopyConfig(supercells={"bulk": (2, 2, 2)},
                                                                        meshes={"bulk": (4, 4, 4)}))
    plan = phonopy_plan(vibrations, "bulk", node.get_pymatgen_structure())
    displacements = plan["displacements"]
    perfect, displaced = displacements.structures(bulk)
    symbols = [site.specie.symbol for site in perfect]
    retrieved = {"perfect": _folder("vasprun.xml", vasprun(symbols, emt_forces([perfect])[0]))}
    for index, forces in enumerate(emt_forces(displaced)):
        retrieved[f"p{index}"] = _folder("vasprun.xml", vasprun(symbols, forces))
    from aiida import orm

    settings = orm.Dict({"backend": "vasp", **phonopy_settings(vibrations, plan)})
    result = phonopy_harmonic_modes._callable(structure=node, settings=settings, **retrieved).get_dict()
    expected = emt_modes(bulk, (2, 2, 2), (4, 4, 4))
    assert result["method"] == "phonopy" and result["supercell_size"] == 64 and result["frozen_composition"] == {}
    assert np.allclose(result["frequencies_cm1"], expected.vibrations.frequencies_cm1, atol=0.05)  # 8-digit vasprun
    vibrations_read = psteros.HarmonicVibrations.from_dict(result)
    assert vibrations_read.free_energy_ev(600.0) == pytest.approx(expected.vibrations.free_energy_ev(600.0), abs=1e-5)

    class Graph:
        outputs = {"bulk_vibrations": result, "o2_vibrations": {"frequencies_cm1": [1580.0], "composition": {"O": 2}}}

    assert psteros.read_phonon_stability(Graph(), ["bulk"])["bulk"].stable
    assert psteros.read_vibrations(Graph(), ["bulk"])["bulk"].supercell_size == 64
    with pytest.raises(ValueError, match="o2: Gamma-point modes, not phonopy"):
        psteros.read_phonon_stability(Graph(), ["bulk", "o2"])
    del retrieved["p0"]
    with pytest.raises(ValueError, match=r"1 of the \d+ retrieved folders are missing \(p0\)"):
        phonopy_harmonic_modes._callable(structure=node, settings=settings, **retrieved)


# ------------------------------------------------------------------ AiiDA: the MLIP route


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_mlip_graph_with_phonopy(tmp_path) -> None:
    pytest.importorskip("aiida_cp2k")
    pytest.importorskip("aiida_shell")
    pytest.importorskip("phonopy")
    from test_mlip_workflow import _computer, _outputs, _recipes, _structures

    computer = _computer(tmp_path)
    static, fine_tune, job = _recipes(tmp_path, computer)
    structures = _structures()
    fixed = {"slab": psteros.central_sites(structures["slab"])}
    config = psteros.PhonopyConfig(slabs=("slab",), supercells={"bulk": (2, 2, 3)})
    vibrations = psteros.VibrationsConfig(molecules=("o2",), phonopy=config)
    graph = psteros.build_mlip_surface_workgraph(structures, static, fine_tune, job=job, metallic_bulk="bulk",
                                                 fixed_sites=fixed, vibrations=vibrations)
    names = {task.name for task in graph.tasks}
    assert {"mlip_vibrations_job", "mlip_vibrations"} | {f"{label}_mlip_vibrations" for label in structures} <= names
    assert not any(name.endswith("_mlip_supercell") for name in names)
    settings = graph.tasks["mlip_vibrations_job"].inputs.settings.value.get_dict()
    assert settings["phonopy"] == {"bulk": [2, 2, 3], "slab": list(config.supercell("slab", structures["slab"]))}
    assert settings["phonopy"]["slab"][2] == 1
    assert settings["molecules"] == {"o2": [0, 1]}
    assert settings["model"] == ["finetune/models/ruo2_round1_seed11.model", "Default"]
    assert list(graph.tasks["mlip_vibrations"].inputs.arguments.value) == ["-u", "forces.py"]
    slab_modes = graph.tasks["slab_mlip_vibrations"].inputs.settings.value.get_dict()
    assert slab_modes == {"label": "slab", "model": "seed11", "method": "phonopy",
                          "supercell": settings["phonopy"]["slab"], "mesh": list(config.mesh("slab", structures["slab"])),
                          "shifted_mesh": True, "symprec": 1e-5}
    assert graph.tasks["o2_mlip_vibrations"].inputs.settings.value.get_dict()["method"] == "gamma"
    links = {(link.from_task.name, link.to_task.name) for link in graph.links}
    assert {("slab_mlip_round1_relaxed", "mlip_vibrations_job"), ("mlip_round1_finetune", "mlip_vibrations"),
            ("mlip_vibrations_job", "mlip_vibrations"), ("mlip_vibrations", "slab_mlip_vibrations"),
            ("mlip_vibrations_job", "slab_mlip_vibrations")} <= links
    assert {f"{label}_mlip_vibrations" for label in structures} | {"mlip_vibrations_results"} <= _outputs(graph)

    with pytest.raises(ValueError, match="bulk: PhonopyConfig.slabs lists it"):
        psteros.build_mlip_surface_workgraph(
            structures, static, fine_tune, job=job, metallic_bulk="bulk",
            vibrations=psteros.VibrationsConfig(molecules=("o2",), phonopy=psteros.PhonopyConfig(slabs=("bulk",))))


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_mlip_phonopy_calcfunctions_run_the_job(tmp_path) -> None:
    import subprocess
    import sys

    pytest.importorskip("phonopy")
    from aiida import orm
    from test_mlip_route import _stub_mace

    from psteros.backends.mlip_tasks import mlip_displacements_folder, mlip_phonon_modes
    from psteros.backends.qe import as_aiida_structure

    bulk, slab = cu3au_bulk(), cu_slab(3)
    from pymatgen.core import Lattice, Structure

    n2 = Structure(Lattice.cubic(9.0), ["N", "N"], [[0.5, 0.5, 0.5], [0.5 + 1.1 / 9.0, 0.5, 0.5]])
    nodes = {label: as_aiida_structure(structure).store() for label, structure in
             (("bulk", bulk), ("slab", slab), ("n2", n2))}
    settings = {"model": ["finetune/models/x_seed1.model", "Default"], "device": "cpu", "default_dtype": "float64",
                "displacement_angstrom": 0.01, "plus_minus": False, "symprec": 1e-5,
                "phonopy": {"bulk": [2, 2, 2], "slab": [3, 3, 1]}, "molecules": {"n2": [0, 1]}}
    folder = mlip_displacements_folder._callable(
        settings=orm.Dict(settings), **{f"structure_{label}": node for label, node in nodes.items()})
    assert {"job.json", "forces.py", "phonopy.json", "bulk.extxyz", "slab.extxyz", "n2.extxyz"} <= set(
        folder.list_object_names())
    work = tmp_path / "work"
    work.mkdir()
    for name in folder.list_object_names():
        (work / name).write_text(folder.get_object_content(name))
    completed = subprocess.run([sys.executable, "-u", "forces.py"], cwd=work, capture_output=True, text=True,
                               env={"PYTHONPATH": _stub_mace(tmp_path), "PATH": "/usr/bin:/bin"}, timeout=600)
    assert completed.returncode == 0, completed.stderr
    folder.store()
    forces = orm.FolderData(tree=str(work / "forces")).store()

    for label, supercell, mesh in (("bulk", (2, 2, 2), (4, 4, 4)), ("slab", (3, 3, 1), (6, 6, 1))):
        options = {"label": label, "model": "seed1", "method": "phonopy", "supercell": list(supercell),
                   "mesh": list(mesh), "shifted_mesh": True, "symprec": 1e-5}
        result = mlip_phonon_modes._callable(structure=nodes[label], settings=orm.Dict(options), job=folder,
                                             forces=forces).get_dict()
        structure = {"bulk": bulk, "slab": slab}[label]
        masses = [nodes[label].get_kind(site.kind_name).mass for site in nodes[label].sites]
        expected = emt_modes(structure, supercell, mesh, masses_amu=masses)
        assert result["model"] == "seed1" and result["method"] == "phonopy"
        assert np.allclose(result["frequencies_cm1"], expected.vibrations.frequencies_cm1, atol=1e-3)
    options = {"label": "n2", "model": "seed1", "method": "gamma", "displaced_sites": [0, 1],
               "displacement_angstrom": 0.01, "supercell_size": 1, "molecule": True}
    result = mlip_phonon_modes._callable(structure=nodes["n2"], settings=orm.Dict(options), job=folder,
                                         forces=forces).get_dict()
    assert result["zero_modes"] == 5 and len(result["frequencies_cm1"]) == 1 and result["frequencies_cm1"][0] > 500
    with pytest.raises(ValueError, match="wrote no ag.npz"):
        mlip_phonon_modes._callable(structure=nodes["n2"], settings=orm.Dict(dict(options, label="ag")), job=folder,
                                    forces=forces)


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_mlip_study_with_phonopy_lets_the_slabs_vibrate_as_a_whole(tmp_path) -> None:
    """The study adds its slabs to PhonopyConfig.slabs and no longer needs a bulk-like frozen centre."""

    pytest.importorskip("aiida_cp2k")
    pytest.importorskip("aiida_shell")
    pytest.importorskip("phonopy")
    from test_mlip_workflow import _computer, _recipes, _ruo2_study

    computer = _computer(tmp_path)
    static, fine_tune, job = _recipes(tmp_path, computer)
    study = _ruo2_study()
    off_centre = psteros.MlipSurfaceStudy(study.bulk, study.references, study.slabs, fixed_half_width_angstrom=0.6)
    with pytest.raises(ValueError, match="does not have the bulk stoichiometry"):
        off_centre.build(static, fine_tune, job=job, vibrations=psteros.VibrationsConfig())
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        graph = off_centre.build(static, fine_tune, job=job,
                                 vibrations=psteros.VibrationsConfig(phonopy=psteros.PhonopyConfig()))
    settings = graph.tasks["mlip_vibrations_job"].inputs.settings.value.get_dict()
    assert set(settings["phonopy"]) == set(study.structures) - {"ref_O"}
    assert settings["molecules"] == {"ref_O": [0, 1]}
    for label in study.slabs:
        assert settings["phonopy"][label][2] == 1


# ------------------------------------------------------------------ examples


def test_emt_example_gives_the_vibrational_surface_energy() -> None:
    import importlib.util
    from pathlib import Path

    pytest.importorskip("phonopy")
    path = Path(__file__).resolve().parents[2] / "examples" / "vibrations_phonopy" / "cu100_emt.py"
    spec = importlib.util.spec_from_file_location("example_cu100_emt", path)
    example = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(example)
    config = psteros.PhonopyConfig(slabs=("slab",), supercells={"bulk": (4, 4, 4), "slab": (4, 4, 1)},
                                   meshes={"bulk": (12, 12, 12), "slab": (12, 12, 1)})
    slab = example.relaxed_slab(3)
    bulk_modes, slab_modes = example.modes("bulk", cu_bulk(), config), example.modes("slab", slab, config)
    gamma = (slab_modes.free_energy_ev(900.0) - 3 * bulk_modes.free_energy_ev(900.0)) / (slab.lattice.a * slab.lattice.b * 2)
    assert gamma == pytest.approx(emt_gamma_vib_900k(), abs=1e-6)


def emt_gamma_vib_900k() -> float:
    """gamma_vib(900 K) of the 3-layer EMT Cu(100) slab with the supercells and meshes of the example test."""

    slab, bulk = cu_slab(3), cu_bulk()
    slab_modes = emt_modes(slab, (4, 4, 1), (12, 12, 1)).vibrations
    bulk_modes = emt_modes(bulk, (4, 4, 4), (12, 12, 12)).vibrations
    return (slab_modes.free_energy_ev(900.0) - 3 * bulk_modes.free_energy_ev(900.0)) / (2 * slab.lattice.a * slab.lattice.b)


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_sno2_mlip_example_builds_the_phonopy_graph(tmp_path) -> None:
    pytest.importorskip("aiida_vasp")
    pytest.importorskip("aiida_shell")
    pytest.importorskip("phonopy")
    import argparse

    from aiida import orm
    from ase import Atoms
    from ase.io import write
    from test_mlip_route import _example

    computer = orm.Computer(label=f"local-{tmp_path.name}", hostname="localhost", transport_type="core.local",
                            scheduler_type="core.direct", workdir=str(tmp_path)).store()
    orm.InstalledCode(computer=computer, filepath_executable="/bin/true", label="vasp",
                      default_calc_job_plugin="vasp.vasp").store()
    orm.InstalledCode(computer=computer, filepath_executable="/opt/mace/bin/python", label="mace",
                      default_calc_job_plugin="core.shell").store()
    frames = [Atoms("SnO", positions=[[0, 0, 0], [2.0 + 0.01 * i, 0, 0]], cell=[8, 8, 8], pbc=True) for i in range(6)]
    write(tmp_path / "replay.xyz", frames, format="extxyz")
    args = argparse.Namespace(code=f"vasp@{computer.label}", mace_code=f"mace@{computer.label}",
                              model="/models/mace-mh-1.model", replay=str(tmp_path / "replay.xyz"), queue="gpu",
                              vibrations=False, phonopy=True, submit=False)
    example = _example("sno2_charge_neutral_mlip")
    route = example.route()
    static, fine_tune, job = example.recipes(args)
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)  # every slab is known to PhonopyConfig: no vacuum warning
        graph = route.build(static, fine_tune, job=job, vibrations=example.phonons())
    settings = graph.tasks["mlip_vibrations_job"].inputs.settings.value.get_dict()
    assert set(settings["phonopy"]) == set(route.structures) - set(route.molecule_labels)
    assert all(settings["phonopy"][label][2] == 1 for label in route.slabs)
    assert set(settings["molecules"]) == set(route.molecule_labels)

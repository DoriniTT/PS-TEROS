"""Tests of the vibrational contributions: harmonic thermodynamics, finite differences, graphs."""

from __future__ import annotations

import math

import numpy as np
import pytest

import psteros
from psteros.vibrations import EV_PER_CM1, KB_EV_PER_K
from psteros.vibrations_workflow import displacement_key, qe_output_forces, vasprun_forces

# sqrt(k / mu) with k in eV/A^2 and mu in amu, as a wavenumber (cm^-1)
CM1_PER_SQRT_EV_PER_A2_AMU = 521.4709


def single_mode(frequency_cm1: float, temperature_k: float) -> tuple[float, float]:
    """Zero-point energy and Helmholtz free energy of one harmonic mode."""

    energy = frequency_cm1 * EV_PER_CM1
    if temperature_k == 0:
        return energy / 2, energy / 2
    return energy / 2, energy / 2 + KB_EV_PER_K * temperature_k * math.log(1 - math.exp(-energy / (KB_EV_PER_K * temperature_k)))


# ---------------------------------------------------------------- thermodynamics


def test_harmonic_thermodynamics_of_known_modes() -> None:
    vibrations = psteros.HarmonicVibrations((200.0, 500.0, 800.0))
    zpe = sum(single_mode(value, 0)[0] for value in (200.0, 500.0, 800.0))
    assert vibrations.zero_point_energy_ev == pytest.approx(zpe)
    assert vibrations.free_energy_ev(0.0) == pytest.approx(zpe)
    assert vibrations.entropy_ev_per_k(0.0) == 0.0
    expected = sum(single_mode(value, 600.0)[1] for value in (200.0, 500.0, 800.0))
    assert vibrations.free_energy_ev(600.0) == pytest.approx(expected, rel=1e-12)
    # F = U - TS, and F decreases with temperature (S > 0)
    assert vibrations.free_energy_ev(600.0) == pytest.approx(
        vibrations.internal_energy_ev(600.0) - 600.0 * vibrations.entropy_ev_per_k(600.0))
    assert vibrations.free_energy_ev(900.0) < vibrations.free_energy_ev(300.0) < zpe
    with pytest.raises(ValueError, match="temperature_k"):
        vibrations.free_energy_ev(-1.0)


def test_agrees_with_ase_harmonic_thermo() -> None:
    thermochemistry = pytest.importorskip("ase.thermochemistry")
    frequencies = (85.0, 230.0, 410.0, 655.0)
    reference = thermochemistry.HarmonicThermo([value * EV_PER_CM1 for value in frequencies])
    # ASE uses an older CODATA k_B: agreement to 1 micro-eV (and 1e-9 eV/K)
    vibrations = psteros.HarmonicVibrations(frequencies)
    for temperature in (100.0, 500.0, 1000.0):
        assert vibrations.free_energy_ev(temperature) == pytest.approx(
            reference.get_helmholtz_energy(temperature, verbose=False), abs=1e-6)
        assert vibrations.entropy_ev_per_k(temperature) == pytest.approx(
            reference.get_entropy(temperature, verbose=False), abs=1e-9)


def test_imaginary_modes_raise_unless_dropped() -> None:
    with pytest.raises(ValueError, match=r"1 imaginary mode\(s\) \(35.0i cm\^-1\)"):
        psteros.HarmonicVibrations((-35.0, 300.0))
    dropped = psteros.HarmonicVibrations((-35.0, 300.0), imaginary_modes="drop")
    assert dropped.imaginary_frequencies_cm1 == (35.0,)
    assert dropped.free_energy_ev(300.0) == pytest.approx(psteros.HarmonicVibrations((300.0,)).free_energy_ev(300.0))
    soft = psteros.HarmonicVibrations((5.0, 300.0), low_frequency_cutoff_cm1=50.0)
    assert soft.free_energy_ev(300.0) == pytest.approx(psteros.HarmonicVibrations((50.0, 300.0)).free_energy_ev(300.0))


def test_supercell_energies_are_per_input_cell() -> None:
    cell = psteros.HarmonicVibrations((100.0, 300.0) * 3, composition={"Sn": 2, "O": 4}, supercell_size=3)
    single = psteros.HarmonicVibrations((100.0, 300.0), composition={"Sn": 2, "O": 4})
    assert cell.free_energy_ev(500.0) == pytest.approx(single.free_energy_ev(500.0))
    assert cell.zero_point_energy_ev == pytest.approx(single.zero_point_energy_ev)
    restored = psteros.HarmonicVibrations.from_dict(cell.to_dict())
    assert restored == cell
    with pytest.raises(ValueError, match="supercell_size"):
        psteros.HarmonicVibrations((100.0,), supercell_size=0)


def test_solid_free_energy_counts_the_frozen_region_as_bulk() -> None:
    bulk = psteros.HarmonicVibrations((150.0, 400.0, 600.0), composition={"Sn": 2, "O": 4})
    slab = psteros.HarmonicVibrations(
        (120.0, 350.0), composition={"Sn": 6, "O": 14}, frozen_composition={"Sn": 4, "O": 8})
    temperature = 700.0
    assert psteros.solid_free_energy_ev(-10.0, bulk, temperature) == pytest.approx(-10.0 + bulk.free_energy_ev(temperature))
    assert psteros.solid_free_energy_ev(-50.0, slab, temperature, bulk=bulk) == pytest.approx(
        -50.0 + slab.free_energy_ev(temperature) + 2 * bulk.free_energy_ev(temperature))
    with pytest.raises(ValueError, match="pass bulk="):
        psteros.solid_free_energy_ev(-50.0, slab, temperature)
    off_stoichiometry = psteros.HarmonicVibrations((120.0,), composition={"Sn": 6, "O": 14},
                                                   frozen_composition={"Sn": 4, "O": 6})
    with pytest.raises(ValueError, match="bulk stoichiometry"):
        psteros.solid_free_energy_ev(-50.0, off_stoichiometry, temperature, bulk=bulk)
    passivated = psteros.HarmonicVibrations((120.0,), frozen_composition={"Ga": 1, "As": 1, "H": 1})
    with pytest.raises(ValueError, match="elements of the bulk"):
        psteros.solid_free_energy_ev(-50.0, passivated, temperature,
                                     bulk=psteros.HarmonicVibrations((100.0,), composition={"Ga": 1, "As": 1}))


def test_molecule_reference_adds_only_the_zero_point_energy() -> None:
    o2 = psteros.HarmonicVibrations((1580.0,))
    assert psteros.molecule_reference_energy_ev(-9.86, o2) == pytest.approx(-9.86 + 1580.0 * EV_PER_CM1 / 2)


def test_bulk_like_vibrations_leave_a_stoichiometric_surface_energy_unchanged() -> None:
    """Free energies feed the existing phase diagram; vibrations equal to the bulk ones cancel."""

    bulk = psteros.HarmonicVibrations((150.0, 400.0, 600.0), composition={"Sn": 2, "O": 4})
    # two bulk cells, of which one is frozen: free sites vibrate as one bulk cell
    slab = psteros.HarmonicVibrations((150.0, 400.0, 600.0), composition={"Sn": 4, "O": 8},
                                      frozen_composition={"Sn": 2, "O": 4})
    o2 = psteros.HarmonicVibrations((1580.0,))
    temperature = 800.0

    def diagram(slab_energy: float, bulk_energy: float, o2_energy: float) -> psteros.SurfacePhaseDiagram:
        references = psteros.BinaryOxideReferences(
            bulk_energy_ev=bulk_energy, bulk_composition={"Sn": 2, "O": 4}, oxygen_molecule_energy_ev=o2_energy)
        termination = psteros.SlabTermination("sto", slab_energy, {"Sn": 4, "O": 8}, 20.0)
        return psteros.surface_phase_diagram([termination], references, delta_mu_range=(-2.0, 0.0), points=5)

    at_zero = diagram(-78.0, -40.0, -9.86)
    with_vibrations = diagram(
        psteros.solid_free_energy_ev(-78.0, slab, temperature, bulk=bulk),
        psteros.solid_free_energy_ev(-40.0, bulk, temperature),
        psteros.molecule_reference_energy_ev(-9.86, o2),
    )
    for old, new in zip(at_zero.curves["sto"], with_vibrations.curves["sto"]):
        assert new.gamma_ev_per_angstrom2 == pytest.approx(old.gamma_ev_per_angstrom2, abs=1e-12)


# ---------------------------------------------------------------- finite differences


def displaced_forces(constants: np.ndarray, sites: list[int], step: float) -> tuple[list, list]:
    """Forces F = -K u of a harmonic model after +step and -step displacements of each site and axis."""

    number_of_sites = constants.shape[0] // 3
    plus, minus = [], []
    for site in sites:
        for axis in range(3):
            for sign, store in ((1.0, plus), (-1.0, minus)):
                displacement = np.zeros(3 * number_of_sites)
                displacement[3 * site + axis] = sign * step
                store.append((-constants @ displacement).reshape(number_of_sites, 3))
    return plus, minus


def test_finite_differences_recover_the_modes_of_a_harmonic_model() -> None:
    rng = np.random.default_rng(7)
    matrix = rng.normal(size=(9, 9))
    constants = matrix @ matrix.T + 9 * np.eye(9)  # positive definite, eV/A^2
    masses = [15.999, 118.71, 15.999]
    plus, minus = displaced_forces(constants, [0, 1, 2], 0.01)
    vibrations = psteros.harmonic_vibrations_from_forces(
        masses_amu=masses, displaced_sites=[0, 1, 2], displacement_angstrom=0.01,
        forces_plus=plus, forces_minus=minus)
    weights = 1 / np.sqrt(np.repeat(masses, 3))
    expected = np.sqrt(np.linalg.eigvalsh(constants * np.outer(weights, weights))) * CM1_PER_SQRT_EV_PER_A2_AMU
    assert np.allclose(sorted(vibrations.frequencies_cm1), sorted(expected), rtol=1e-4)


def test_partial_hessian_of_free_sites() -> None:
    rng = np.random.default_rng(3)
    matrix = rng.normal(size=(9, 9))
    constants = matrix @ matrix.T + 9 * np.eye(9)
    masses = [15.999, 118.71, 15.999]
    plus, minus = displaced_forces(constants, [1, 2], 0.01)
    vibrations = psteros.harmonic_vibrations_from_forces(
        masses_amu=masses, displaced_sites=[1, 2], displacement_angstrom=0.01,
        forces_plus=plus, forces_minus=minus, composition={"Sn": 1, "O": 2}, frozen_composition={"O": 1})
    block = constants[3:, 3:]
    weights = 1 / np.sqrt(np.repeat(masses[1:], 3))
    expected = np.sqrt(np.linalg.eigvalsh(block * np.outer(weights, weights))) * CM1_PER_SQRT_EV_PER_A2_AMU
    assert len(vibrations.frequencies_cm1) == 6
    assert np.allclose(sorted(vibrations.frequencies_cm1), sorted(expected), rtol=1e-4)
    assert vibrations.frozen_composition == {"O": 1}


def diatomic_constants(spring: float) -> np.ndarray:
    """Force constants of a diatomic along x, the 3 translations and 2 rotations at zero frequency."""

    constants = np.zeros((6, 6))
    constants[0, 0] = constants[3, 3] = spring
    constants[0, 3] = constants[3, 0] = -spring
    return constants


def test_molecule_keeps_only_its_stretch() -> None:
    spring, mass = 70.0, 15.999
    plus, minus = displaced_forces(diatomic_constants(spring), [0, 1], 0.01)
    vibrations = psteros.harmonic_vibrations_from_forces(
        masses_amu=[mass, mass], displaced_sites=[0, 1], displacement_angstrom=0.01,
        forces_plus=plus, forces_minus=minus, zero_modes=5)
    assert vibrations.frequencies_cm1 == pytest.approx((math.sqrt(spring / (mass / 2)) * CM1_PER_SQRT_EV_PER_A2_AMU,), rel=1e-4)
    with pytest.raises(ValueError, match="shape"):
        psteros.harmonic_vibrations_from_forces(
            masses_amu=[mass, mass], displaced_sites=[0], displacement_angstrom=0.01,
            forces_plus=plus, forces_minus=minus)


# ---------------------------------------------------------------- output parsers


def vasprun(symbols: list[str], forces) -> str:
    atoms = "".join(f"<rc><c>{symbol:<2}</c><c>   1</c></rc>" for symbol in symbols)
    rows = "".join(f"<v> {f[0]:.8f} {f[1]:.8f} {f[2]:.8f} </v>" for f in forces)
    zero = "".join("<v> 0.0 0.0 0.0 </v>" for _ in symbols)
    return (
        "<?xml version='1.0' encoding='ISO-8859-1'?><modeling>"
        f"<atominfo><atoms>{len(symbols)}</atoms><array name='atoms'><set>{atoms}</set></array></atominfo>"
        f"<calculation><varray name='forces'>{zero}</varray></calculation>"
        f"<calculation><varray name='forces'>{rows}</varray></calculation></modeling>"
    )


def pw_output(forces) -> str:
    lines = ["     Forces acting on atoms (cartesian axes, Ry/au):", ""]
    ry_per_bohr = 13.605693122994 / 0.529177210903
    lines += [
        f"     atom {i + 1:4d} type  1   force = {f[0] / ry_per_bohr:14.8f}{f[1] / ry_per_bohr:14.8f}{f[2] / ry_per_bohr:14.8f}"
        for i, f in enumerate(forces)
    ]
    lines += ["     The non-local contrib.  to forces"]
    lines += [f"     atom {i + 1:4d} type  1   force =      9.0     9.0     9.0" for i in range(len(forces))]
    return "\n".join(lines) + "\n     Total force =     0.001\n"


def test_force_parsers() -> None:
    forces = [[0.1, -0.2, 0.3], [-0.1, 0.2, -0.3]]
    symbols, parsed = vasprun_forces(vasprun(["Sn", "O"], forces))
    assert symbols == ["Sn", "O"] and np.allclose(parsed, forces)
    assert np.allclose(qe_output_forces(pw_output(forces), 2), forces, atol=1e-6)
    with pytest.raises(ValueError, match="tprnfor"):
        qe_output_forces("JOB DONE.", 2)
    with pytest.raises(ValueError, match="1 of 2"):
        qe_output_forces(pw_output(forces[:1]), 2)


# ---------------------------------------------------------------- recipe of the graph


def test_vibrations_config() -> None:
    config = psteros.VibrationsConfig(displaced_sites={"slab": (3, 1)}, supercells={"bulk": (2, 2, 1)},
                                      molecules=("o2",))
    assert config.displaced_sites == {"slab": (1, 3)}
    assert config.sites_to_displace("slab", 5) == (1, 3)
    assert config.sites_to_displace("other", 4, psteros.CalculationOverride(fixed_sites=(0, 2))) == (1, 3)
    assert config.sites_to_displace("bulk", 3) == (0, 1, 2)
    with pytest.raises(ValueError, match="outside"):
        config.sites_to_displace("slab", 2)
    with pytest.raises(ValueError, match="supercell displaces every site"):
        config.sites_to_displace("bulk", 3, psteros.CalculationOverride(fixed_sites=(0,)))
    with pytest.raises(ValueError, match="nothing to displace"):
        config.sites_to_displace("x", 1, psteros.CalculationOverride(fixed_sites=(0,)))
    with pytest.raises(ValueError, match="too large"):
        psteros.VibrationsConfig(displacement_angstrom=0.5)
    with pytest.raises(ValueError, match="three positive integers"):
        psteros.VibrationsConfig(supercells={"bulk": (2, 2)})
    with pytest.raises(ValueError, match="do not give displaced_sites"):
        psteros.VibrationsConfig(displaced_sites={"bulk": (0,)}, supercells={"bulk": (2, 2, 2)})


def test_read_vibrations_from_a_graph() -> None:
    class Graph:
        outputs = {"slab_vibrations": {"frequencies_cm1": [-20.0, 300.0], "composition": {"Sn": 1, "O": 2}}}

    with pytest.raises(ValueError, match="slab: 1 imaginary mode"):
        psteros.read_vibrations(Graph(), ["slab"])
    vibrations = psteros.read_vibrations(Graph(), ["slab"], imaginary_modes="drop")["slab"]
    assert vibrations.frequencies_cm1 == (-20.0, 300.0) and vibrations.composition == {"Sn": 1, "O": 2}
    with pytest.raises(ValueError, match="no vibrations output for bulk"):
        psteros.read_vibrations(Graph(), ["bulk"])


# ---------------------------------------------------------------- AiiDA


def _computer(tmp_path, label: str, plugin: str) -> str:
    from aiida import orm

    computer = orm.Computer(
        label=f"local-{tmp_path.name}", hostname="localhost", transport_type="core.local",
        scheduler_type="core.direct", workdir=str(tmp_path),
    ).store()
    orm.InstalledCode(computer=computer, filepath_executable="/bin/true", label=label,
                      default_calc_job_plugin=plugin).store()
    return computer.label


def _folder(name: str, content: str):
    import io

    from aiida import orm

    folder = orm.FolderData()
    folder.base.repository.put_object_from_filelike(io.BytesIO(content.encode()), name)
    return folder


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_vasp_vibrations_workgraph(tmp_path) -> None:
    pytest.importorskip("aiida_vasp")
    computer = _computer(tmp_path, "vasp", "vasp.vasp")
    incar = {"ENCUT": 520, "EDIFF": 1e-8, "NSW": 0, "IBRION": -1}
    slab, _ = psteros.sno2_110_slab(termination="o", triple_layers=3)
    fixed = tuple(psteros.central_sites(slab, half_width=1.5))
    static = psteros.SurfaceWorkflowConfig(
        backend="vasp", name="sno2", role_overrides={"slab": psteros.CalculationOverride(fixed_sites=fixed)},
        execution=psteros.ExecutionPolicy(computer=computer, max_concurrent_jobs=4),
        calculation=psteros.VaspCalculationConfig(f"vasp@{computer}", incar, potential_mapping={"Sn": "Sn_d"}),
    )
    o2 = psteros.triplet_o2_cell()
    graph = psteros.build_vibrations_workgraph(
        {"slab": slab, "bulk": psteros.rutile_sno2_bulk(), "o2": o2}, static,
        psteros.VibrationsConfig(supercells={"bulk": (1, 1, 2)}, molecules=("o2",)),
    )
    assert graph.max_number_jobs == 4
    names = {task.name for task in graph.tasks}
    free = len(slab) - len(fixed)
    assert sum(name.startswith("slab_vib_") for name in names) == 6 * free
    assert sum(name.startswith("bulk_vib_") for name in names) == 6 * 12
    assert sum(name.startswith("o2_vib_") for name in names) == 12
    assert {"slab_vibrations", "bulk_vibrations", "o2_vibrations", "bulk_supercell"} <= names
    free_site = next(index for index in range(len(slab)) if index not in fixed)
    vasp = graph.tasks[f"slab_vib_{displacement_key(free_site, 0, 'plus')}_vasp"]
    # the displaced static calculations do not inherit the selective dynamics of the relaxation
    assert vasp.inputs.parameters.value.get_dict() == {"incar": incar}
    settings = graph.tasks["slab_vibrations"].inputs.settings.value.get_dict()
    assert len(settings["displaced_sites"]) == free and not set(settings["displaced_sites"]) & set(fixed)
    assert settings["supercell_size"] == 1 and settings["molecule"] is False
    assert graph.tasks["bulk_vibrations"].inputs.settings.value.get_dict()["supercell_size"] == 2
    assert graph.tasks["o2_vibrations"].inputs.settings.value.get_dict()["molecule"] is True
    outputs = set(graph.outputs._get_keys()) if hasattr(graph.outputs, "_get_keys") else set(graph.outputs)
    assert {"slab_vibrations", "bulk_vibrations", "o2_vibrations"} <= outputs

    relaxing = psteros.SurfaceWorkflowConfig(
        backend="vasp", execution=static.execution,
        calculation=psteros.VaspCalculationConfig(f"vasp@{computer}", {"ENCUT": 520, "IBRION": 2, "NSW": 50}))
    with pytest.raises(ValueError, match="static INCAR must not relax"):
        psteros.build_vibrations_workgraph({"o2": o2}, relaxing)
    with pytest.raises(ValueError, match="not in structures"):
        psteros.build_vibrations_workgraph({"o2": o2}, static, psteros.VibrationsConfig(molecules=("h2o",)))


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_qe_vibrations_workgraph_prints_forces(tmp_path) -> None:
    pytest.importorskip("aiida_quantumespresso")
    from tests.unit.test_surface_study import _pseudo_family

    computer = _computer(tmp_path, "pw", "quantumespresso.pw")
    static = psteros.SurfaceWorkflowConfig(
        backend="qe",
        calculation=psteros.QeCalculationConfig(
            f"pw@{computer}", _pseudo_family(tmp_path, ["Sn", "O"]),
            {"CONTROL": {"calculation": "scf"}, "SYSTEM": {"ecutwfc": 40.0}, "ELECTRONS": {"conv_thr": 1e-10}}),
        execution=psteros.ExecutionPolicy(computer=computer),
        role_overrides={"o2": psteros.CalculationOverride(parameters={"SYSTEM": {"nspin": 2}}, fixed_sites=(0,))},
    )
    graph = psteros.build_vibrations_workgraph({"o2": psteros.triplet_o2_cell()}, static)
    assert graph.max_number_jobs == 1
    task = next(task for task in graph.tasks if task.name.endswith("_qe"))
    parameters = task.inputs.pw.parameters.value.get_dict()
    assert parameters["CONTROL"] == {"calculation": "scf", "tprnfor": True}
    assert parameters["SYSTEM"] == {"ecutwfc": 40.0, "nspin": 2}
    assert "settings" not in task.inputs.pw or task.inputs.pw.settings.value is None
    assert sum(task.name.startswith("o2_vib_") for task in graph.tasks) == 6  # site 0 is fixed


@pytest.mark.tier2
@pytest.mark.requires_aiida
@pytest.mark.parametrize("backend", ["vasp", "qe"])
def test_harmonic_modes_calcfunction(backend) -> None:
    from aiida import orm

    from psteros.backends.vibration_tasks import displace_site, harmonic_modes, repeat_structure

    spring, mass, step = 70.0, 15.999, 0.01
    molecule = orm.StructureData(cell=[[10.0, 0, 0], [0, 10.0, 0], [0, 0, 10.0]])
    molecule.append_atom(position=(5.0, 5.0, 5.0), symbols="O")
    molecule.append_atom(position=(6.2, 5.0, 5.0), symbols="O")
    plus, minus = displaced_forces(diatomic_constants(spring), [0, 1], step)
    write = (lambda f: _folder("vasprun.xml", vasprun(["O", "O"], f))) if backend == "vasp" else (
        lambda f: _folder("aiida.out", pw_output(f)))
    retrieved = {}
    for index, (site, axis) in enumerate((s, a) for s in (0, 1) for a in range(3)):
        retrieved[displacement_key(site, axis, "plus")] = write(plus[index])
        retrieved[displacement_key(site, axis, "minus")] = write(minus[index])
    settings = orm.Dict({"backend": backend, "displaced_sites": [0, 1], "displacement_angstrom": step,
                         "supercell_size": 1, "molecule": True})
    result = harmonic_modes._callable(structure=molecule, settings=settings, **retrieved)
    vibrations = psteros.HarmonicVibrations.from_dict(result)
    assert result["zero_modes"] == 5 and vibrations.composition == {"O": 2}
    assert vibrations.frequencies_cm1 == pytest.approx(
        (math.sqrt(spring / (mass / 2)) * CM1_PER_SQRT_EV_PER_A2_AMU,), rel=1e-3)

    moved = displace_site._callable(structure=molecule, site=orm.Int(1), vector=orm.List([0.0, step, 0.0]))
    assert np.allclose(moved.sites[1].position, (6.2, 5.0 + step, 5.0)) and moved.sites[0].position == (5.0, 5.0, 5.0)
    supercell = repeat_structure._callable(structure=molecule, size=orm.List([1, 1, 2]))
    assert len(supercell.sites) == 4 and supercell.cell[2][2] == 20.0


# ---------------------------------------------------------------- example


def _load_example(name: str):
    import importlib.util
    import sys
    from pathlib import Path

    directory = Path(__file__).resolve().parents[2] / "examples" / "vasp_surface_phase_diagram"
    sys.path.insert(0, str(directory))
    try:
        spec = importlib.util.spec_from_file_location(f"vibrations_example_{name}", directory / f"{name}.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(directory))
    return module


def test_example_phase_diagram_with_vibrations() -> None:
    analysis = _load_example("phase_diagram")
    bulk, metal = psteros.rutile_sno2_bulk(), psteros.alpha_sn_bulk()

    class Structure:
        def __init__(self, structure):
            self.structure = structure

        def get_pymatgen_structure(self):
            return self.structure

    class Graph:
        def __init__(self, outputs):
            self.outputs = outputs

    misc = lambda energy: {"total_energies": {"energy_extrapolated": energy}}  # noqa: E731
    refs = Graph({
        "sno2_bulk_static_misc": misc(-40.0), "sno2_bulk_relaxed_structure": Structure(bulk),
        "alpha_sn_static_misc": misc(-32.0), "alpha_sn_relaxed_structure": Structure(metal),
        "o2_static_misc": misc(-9.8), "o2_relaxed_structure": Structure(psteros.triplet_o2_cell()),
    })
    outputs, bulk_modes = {}, (150.0, 400.0, 600.0)
    vibrations = {
        "sno2_bulk_vibrations": {"frequencies_cm1": bulk_modes, "composition": {"Sn": 2, "O": 4}},
        "alpha_sn_vibrations": {"frequencies_cm1": (120.0,), "composition": {"Sn": 8}},
        "o2_vibrations": {"frequencies_cm1": (1580.0,), "composition": {"O": 2}},
    }
    for termination, label in zip(("o", "sno", "sn2o"), analysis.SLABS):
        slab, _ = psteros.sno2_110_slab(termination=termination, triple_layers=3)
        n_sn, n_o = int(slab.composition["Sn"]), int(slab.composition["O"])
        outputs[f"{label}_static_misc"] = misc(-20.0 * n_sn / 2 - 4.9 * (n_o - 2 * n_sn) + 3.0)
        outputs[f"{label}_relaxed_structure"] = slab
        # free sites vibrate like the bulk cells they contain; the frozen centre is one bulk cell
        vibrations[f"{label}_vibrations"] = {
            "frequencies_cm1": bulk_modes * (n_sn // 2 - 1), "composition": dict(slab.composition.get_el_amt_dict()),
            "frozen_composition": {"Sn": 2, "O": 4},
        }
    at_zero = analysis.build_diagram(refs, Graph(outputs))
    hot = analysis.build_diagram(refs, Graph(outputs), Graph(vibrations), 800.0)
    # Sn6O12 is stoichiometric: with bulk-like vibrations its gamma does not change
    assert hot.curves["slab_o"][-1].gamma_ev_per_angstrom2 == pytest.approx(
        at_zero.curves["slab_o"][-1].gamma_ev_per_angstrom2, abs=1e-12)
    assert hot.references.oxygen_poor_limit_ev != at_zero.references.oxygen_poor_limit_ev
    with pytest.raises(SystemExit):
        analysis.main(["--profile", "p", "--refs-pk", "1", "--slabs-pk", "2", "--vib-pk", "3"])


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_example_vibrations_graph(tmp_path) -> None:
    pytest.importorskip("aiida_vasp")
    import argparse

    example = _load_example("vibrations")
    example.relaxed_bulk_lattice = lambda pk: (4.737, 3.186)
    code = f"vasp@{_computer(tmp_path, 'vasp', 'vasp.vasp')}"
    args = argparse.Namespace(code=code, potential_family="PBE", computer="local", queue=None, ranks=1,
                              walltime=600, parallel_jobs=2)
    fixed = example.fixed_sites(refs_pk=0)
    assert set(fixed) == {"slab_o", "slab_sno", "slab_sn2o"} and all(len(sites) == 6 for sites in fixed.values())
    static = example.static_recipe(args, fixed)
    assert static.calculation.incar["EDIFF"] == 1e-7 and static.execution.max_concurrent_jobs == 2
    slab, _ = psteros.sno2_110_slab(termination="sn2o", triple_layers=3, a=4.737, c=3.186)
    graph = psteros.build_vibrations_workgraph(
        {"slab_sn2o": slab, "o2": psteros.triplet_o2_cell()}, static, psteros.VibrationsConfig(molecules=("o2",)))
    assert sum(task.name.endswith("_vasp") for task in graph.tasks) == 6 * (len(slab) - 6) + 12
    o2 = next(task for task in graph.tasks if task.name.startswith("o2_vib_") and task.name.endswith("_vasp"))
    assert o2.inputs.parameters.value.get_dict()["incar"]["ISPIN"] == 2

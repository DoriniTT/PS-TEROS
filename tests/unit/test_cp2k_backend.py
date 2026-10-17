"""CP2K backend: recipe validation, input assembly, parsers, readers and graphs (submit=False)."""

from __future__ import annotations

import math

import numpy as np
import pytest

import psteros
from psteros.backends.cp2k import (
    BOHR_TO_ANGSTROM,
    DOS_FILE,
    FORCES_FILE,
    HARTREE_TO_EV,
    cp2k_kpoints_mesh,
    cp2k_parameters,
    parse_cp2k_dos,
    parse_cp2k_forces,
)

DFT = {
    "QS": {"METHOD": "GPW", "EPS_DEFAULT": 1.0e-12},
    "MGRID": {"CUTOFF": 600, "REL_CUTOFF": 60},
    "XC": {"XC_FUNCTIONAL": {"_": "PBE"}},
    "SCF": {"ADDED_MOS": 20, "SMEAR": {"METHOD": "FERMI_DIRAC", "ELECTRONIC_TEMPERATURE": 300}},
}
BASIS = {"Ru": "DZVP-MOLOPT-PBE-GTH-q16", "O": "DZVP-MOLOPT-PBE-GTH-q6", "Sn": "DZVP-MOLOPT-SR-GTH-q4"}
POTENTIALS = {"Ru": "GTH-PBE-q16", "O": "GTH-PBE-q6", "Sn": "GTH-PBE-q4"}


def cp2k_config(code_label: str = "cp2k@cluster", run_type: str | None = None, **changes):
    parameters = {"FORCE_EVAL": {"METHOD": "Quickstep", "DFT": DFT}}
    if run_type:
        parameters["GLOBAL"] = {"RUN_TYPE": run_type}
        if run_type in ("GEO_OPT", "CELL_OPT"):
            parameters["MOTION"] = {run_type: {"MAX_FORCE": 1.94e-4}}
    values = dict(code_label=code_label, parameters=parameters, basis_sets=BASIS, potentials=POTENTIALS,
                  basis_set_file="BASIS_MOLOPT_UZH", potential_file="POTENTIAL_UZH")
    values.update(changes)
    return psteros.Cp2kCalculationConfig(**values)


def test_cp2k_recipe_validation() -> None:
    with pytest.raises(ValueError, match="FORCE_EVAL/DFT"):
        cp2k_config(parameters={"FORCE_EVAL": {}})
    with pytest.raises(ValueError, match="kpoints_spacing"):
        cp2k_config(parameters={"FORCE_EVAL": {"DFT": {"KPOINTS": {}}}})
    with pytest.raises(ValueError, match="CELL"):
        cp2k_config(parameters={"FORCE_EVAL": {"DFT": {}, "SUBSYS": {"CELL": {"ABC": "1 1 1"}}}})
    with pytest.raises(ValueError, match="same elements"):
        cp2k_config(potentials={"Ru": "GTH-PBE-q16"})
    with pytest.raises(ValueError, match="kpoints_spacing must be positive"):
        cp2k_config(kpoints_spacing=0.0)
    assert cp2k_config(kpoints_spacing=None).kpoints_spacing is None
    execution = psteros.ExecutionPolicy()
    config = psteros.SurfaceWorkflowConfig(backend="cp2k", calculation=cp2k_config(), execution=execution)
    assert config.backend == "cp2k"
    with pytest.raises(TypeError, match="Cp2kCalculationConfig"):
        psteros.SurfaceWorkflowConfig(
            backend="cp2k", calculation=psteros.VaspCalculationConfig("vasp@cluster"), execution=execution)
    with pytest.raises(TypeError, match="QeCalculationConfig"):
        psteros.SurfaceWorkflowConfig(backend="qe", calculation=cp2k_config(), execution=execution)
    with pytest.raises(ValueError, match=r"or 'cp2k' with Cp2kCalculationConfig\), got 'siesta'"):
        psteros.SurfaceWorkflowConfig(backend="siesta", calculation=cp2k_config(), execution=execution)


def test_static_input_gets_kinds_files_and_prints() -> None:
    config = cp2k_config()
    override = psteros.CalculationOverride(
        parameters={"FORCE_EVAL": {"DFT": {"UKS": True, "MULTIPLICITY": 3}}}, fixed_sites=(0,))
    parameters = cp2k_parameters(config, override, kinds={"O": "O"}, label="o2")
    assert parameters["GLOBAL"] == {"RUN_TYPE": "ENERGY_FORCE", "PRINT_LEVEL": "MEDIUM"}
    dft = parameters["FORCE_EVAL"]["DFT"]
    assert dft["UKS"] is True and dft["MULTIPLICITY"] == 3 and dft["MGRID"]["CUTOFF"] == 600
    assert dft["BASIS_SET_FILE_NAME"] == "BASIS_MOLOPT_UZH" and dft["POTENTIAL_FILE_NAME"] == "POTENTIAL_UZH"
    assert parameters["FORCE_EVAL"]["SUBSYS"]["KIND"] == [
        {"_": "O", "ELEMENT": "O", "BASIS_SET": "DZVP-MOLOPT-PBE-GTH-q6", "POTENTIAL": "GTH-PBE-q6"}]
    assert parameters["FORCE_EVAL"]["PRINT"]["FORCES"]["FILENAME"] == f"={FORCES_FILE}"
    assert dft["PRINT"]["DOS"]["FILENAME"] == f"={DOS_FILE}"
    assert "MOTION" not in parameters  # fixed sites only matter for a relaxation
    assert "UKS" not in config.parameters["FORCE_EVAL"]["DFT"]  # the recipe is not changed


def test_relaxation_input_fixes_sites_one_based() -> None:
    config = cp2k_config(run_type="GEO_OPT")
    parameters = cp2k_parameters(config, psteros.CalculationOverride(fixed_sites=(2, 0)),
                                 kinds={"Ru": "Ru", "O": "O"}, label="slab")
    assert parameters["MOTION"]["CONSTRAINT"]["FIXED_ATOMS"] == {"LIST": "1 3", "COMPONENTS_TO_FIX": "XYZ"}
    assert parameters["MOTION"]["GEO_OPT"] == {"MAX_FORCE": 1.94e-4}
    assert "PRINT" not in parameters["FORCE_EVAL"]  # no force or DOS file during a relaxation
    both = psteros.CalculationOverride(
        fixed_sites=(0,), parameters={"MOTION": {"CONSTRAINT": {"FIXED_ATOMS": {"LIST": "2"}}}})
    with pytest.raises(ValueError, match="not both"):
        cp2k_parameters(config, both, kinds={"Ru": "Ru"}, label="slab")


def test_kinds_must_be_elements_with_basis_sets() -> None:
    config = cp2k_config()
    assert cp2k_parameters(config, None, kinds={"Ru1": "Ru"})["FORCE_EVAL"]["SUBSYS"]["KIND"][0]["_"] == "Ru1"
    with pytest.raises(ValueError, match="slab: CP2K kind 'H0p75'"):
        cp2k_parameters(config, None, kinds={"H0p75": "H"}, label="slab")
    with pytest.raises(ValueError, match=r"slab: no CP2K basis set or potential for \['Ti'\]"):
        cp2k_parameters(config, None, kinds={"Ti": "Ti"}, label="slab")
    given = {"FORCE_EVAL": {"SUBSYS": {"KIND": [{"_": "Ru", "BASIS_SET": "mine"}]}}}
    kinds = cp2k_parameters(config, psteros.CalculationOverride(parameters=given), kinds={"Ru": "Ru"})
    assert kinds["FORCE_EVAL"]["SUBSYS"]["KIND"] == [{"_": "Ru", "BASIS_SET": "mine"}]


def test_kpoint_mesh_uses_two_pi_spacing() -> None:
    assert cp2k_kpoints_mesh(np.diag([4.0, 4.0, 4.0]), 0.2) == (8, 8, 8)
    assert cp2k_kpoints_mesh(np.diag([6.36, 3.11, 40.0]), 0.2) == (5, 11, 1)
    assert cp2k_kpoints_mesh(np.diag([4.0, 4.0, 4.0]), None) == (1, 1, 1)
    assert cp2k_kpoints_mesh(np.diag([2 * math.pi / 2.0] * 3), 1.0) == (2, 2, 2)


def test_one_kpoint_across_vacuum() -> None:
    from psteros.backends.cp2k import vacuum_gaps

    slab, _ = psteros.rutile_110_slab(metal="Ru", termination="o_bridge", a=4.492, c=3.107, vacuum_angstrom=22.0)
    gaps = vacuum_gaps(slab.lattice.matrix, slab.frac_coords)
    assert gaps[2] >= 22.0 and max(gaps[:2]) < 3.0  # pymatgen pads the vacuum to whole layers
    assert cp2k_kpoints_mesh(slab.lattice.matrix, 0.15) == (14, 7, 2)
    assert cp2k_kpoints_mesh(slab.lattice.matrix, 0.15, frac_coords=slab.frac_coords, vacuum_gap_angstrom=8.0) == (14, 7, 1)
    o2 = psteros.triplet_o2_cell(cell_length=12.0)
    assert cp2k_kpoints_mesh(o2.lattice.matrix, 0.2, frac_coords=o2.frac_coords, vacuum_gap_angstrom=8.0) == (1, 1, 1)
    bulk = psteros.rutile_bulk("Ru", a=4.492, c=3.107)
    assert cp2k_kpoints_mesh(bulk.lattice.matrix, 0.2, frac_coords=bulk.frac_coords, vacuum_gap_angstrom=8.0) == (7, 7, 11)
    with pytest.raises(ValueError, match="vacuum_gap_angstrom must be positive"):
        cp2k_config(vacuum_gap_angstrom=0.0)


FORCES_AU = """
 ATOMIC FORCES in [a.u.]

 # Atom   Kind   Element          X                      Y                      Z
      1      1      Ru          0.0100000000000000     0.0000000000000000    -0.0200000000000000
      2      2      O          -0.0100000000000000     0.0050000000000000     0.0200000000000000
 SUM OF ATOMIC FORCES           0.0000000000000000     0.0050000000000000     0.0000000000000000     0.0050000
"""


def test_parse_cp2k_forces_reads_the_unit_from_the_header() -> None:
    forces = parse_cp2k_forces(FORCES_AU)
    scale = HARTREE_TO_EV / BOHR_TO_ANGSTROM
    assert forces.shape == (2, 3)
    assert forces[0] == pytest.approx([0.01 * scale, 0.0, -0.02 * scale])
    assert forces[0, 0] == pytest.approx(0.514220674763, rel=1e-9)
    in_ev = FORCES_AU.replace("[a.u.]", "[eV/Angstrom]")
    assert parse_cp2k_forces(in_ev)[1] == pytest.approx([-0.01, 0.005, 0.02])
    with pytest.raises(ValueError, match="one ATOMIC FORCES block"):
        parse_cp2k_forces(FORCES_AU + FORCES_AU)
    with pytest.raises(ValueError, match="in order"):
        parse_cp2k_forces(FORCES_AU.replace("      2      2      O", "      3      2      O"))
    with pytest.raises(ValueError, match="unknown unit"):
        parse_cp2k_forces(FORCES_AU.replace("[a.u.]", "[kcal]"))


def _dos(rows, header=" # DOS at iteration step i = 12", spins=1) -> str:
    title = (" #  Energy[a.u.]       Density     Occupation" if spins == 1
             else " #  Energy[a.u.]  Alpha_Density    Occupation   Beta_Density     Occupation")
    lines = [header, title] + ["".join(f"{value:15.8f}" for value in row) for row in rows]
    return "\n".join(lines) + "\n"


def test_parse_cp2k_dos_counts_partially_filled_bins() -> None:
    # k-point file: weights 0.25 per state, 4 states per spin in all; full occupation 2 (RKS).
    metal = _dos([(-0.30, 0.25, 0.50), (-0.20, 0.25, 0.25), (-0.10, 0.50, 0.0)])
    dos = parse_cp2k_dos(metal)
    assert dos.energies_ev[0] == pytest.approx(-0.30 * HARTREE_TO_EV)
    assert dos.fermi_energy_ev is None
    assert [ratio.tolist() for ratio in dos.filling()] == [pytest.approx([1.0, 0.5, 0.0])]
    assert dos.partially_occupied_bins() == 1
    insulator = _dos([(-0.30, 0.50, 1.0), (-0.20, 0.0, 0.0), (0.10, 0.50, 0.0)])
    assert parse_cp2k_dos(insulator).partially_occupied_bins() == 0
    gamma = _dos([(-0.30, 0.5, 1.0), (0.1, 0.5, 0.0)],
                 header=" # DOS at iteration step i = 3, E_Fermi[a.u.] =   -0.100000")
    assert parse_cp2k_dos(gamma).fermi_energy_ev == pytest.approx(-0.1 * HARTREE_TO_EV)
    spin = _dos([(-0.3, 0.5, 0.5, 0.5, 0.5), (-0.1, 0.5, 0.5, 0.5, 0.0)], spins=2)
    assert parse_cp2k_dos(spin).partially_occupied_bins() == 0
    assert parse_cp2k_dos(spin + _dos([(-0.3, 0.5, 0.5, 0.5, 0.25)], spins=2)).partially_occupied_bins() == 1
    with pytest.raises(ValueError, match="no 'DOS at iteration step'"):
        parse_cp2k_dos("nothing here")


class _Folder:
    def __init__(self, files):
        self.base = self
        self.repository = self
        self._files = files

    def get_object_content(self, name):
        if name not in self._files:
            raise FileNotFoundError(name)
        return self._files[name]


def test_readers_convert_units_and_prefer_the_static_stage() -> None:
    class Graph:
        outputs = {
            "a_relax_parameters": {"energy": -9.0, "energy_units": "a.u."},
            "a_static_parameters": {"energy": -10.0, "energy_units": "a.u."},
            "a_relaxed_structure": "relaxed-a",
            "a_static_retrieved": _Folder({FORCES_FILE: FORCES_AU, DOS_FILE: _dos([(-0.3, 0.5, 1.0)])}),
            "b_parameters": {"energy": -1.0},
            "b_retrieved": _Folder({}),
        }

    energies, structures = psteros.read_cp2k_results(Graph(), ["a", "b"])
    assert energies == {"a": pytest.approx(-10.0 * HARTREE_TO_EV), "b": pytest.approx(-HARTREE_TO_EV)}
    assert structures == {"a": "relaxed-a"}
    assert psteros.read_cp2k_forces(Graph(), ["a"])["a"].shape == (2, 3)
    assert psteros.read_cp2k_dos(Graph(), ["a"])["a"].partially_occupied_bins() == 0
    with pytest.raises(ValueError, match="b: psteros_forces.dat was not retrieved"):
        psteros.read_cp2k_forces(Graph(), ["b"])
    with pytest.raises(ValueError, match="no output parameters for c"):
        psteros.read_cp2k_results(Graph(), ["c"])


def _structures():
    pytest.importorskip("pymatgen")
    from pymatgen.core import Lattice, Structure

    ruo2 = Structure.from_spacegroup("P4_2/mnm", Lattice.tetragonal(4.49, 3.11), ["Ru", "O"],
                                     [[0, 0, 0], [0.3056, 0.3056, 0]])
    slab, _ = psteros.sno2_110_slab(termination="o", triple_layers=3)
    return {"bulk": ruo2, "slab": slab, "o2": psteros.triplet_o2_cell()}


def test_qe_and_vasp_builders_refuse_cp2k_recipes() -> None:
    structures = _structures()
    execution = psteros.ExecutionPolicy()
    static = psteros.SurfaceWorkflowConfig(backend="cp2k", calculation=cp2k_config(), execution=execution)
    relax = psteros.SurfaceWorkflowConfig(backend="cp2k", calculation=cp2k_config(run_type="GEO_OPT"),
                                          execution=execution)
    with pytest.raises(ValueError, match="build_cp2k_surface_workgraph"):
        psteros.build_surface_workgraph(structures, static)
    with pytest.raises(ValueError, match="build_cp2k_relax_static_workgraph"):
        psteros.build_relax_static_workgraph(structures, relax, static)
    with pytest.raises(ValueError, match="CP2K is not supported yet"):
        psteros.build_vibrations_workgraph(structures, static)


def test_cp2k_relax_static_checks_run_types() -> None:
    structures = _structures()
    execution = psteros.ExecutionPolicy()
    relax = psteros.SurfaceWorkflowConfig(backend="cp2k", calculation=cp2k_config(run_type="GEO_OPT"),
                                          execution=execution)
    static = psteros.SurfaceWorkflowConfig(backend="cp2k", calculation=cp2k_config(), execution=execution)
    with pytest.raises(ValueError, match="bulk: the CP2K relaxation needs GLOBAL/RUN_TYPE GEO_OPT or CELL_OPT"):
        psteros.build_cp2k_relax_static_workgraph(structures, static, static)
    with pytest.raises(ValueError, match="bulk: the CP2K static calculation must not move atoms"):
        psteros.build_cp2k_relax_static_workgraph(structures, relax, relax)
    with pytest.raises(ValueError, match="share one execution policy"):
        psteros.build_cp2k_relax_static_workgraph(
            structures, relax, psteros.SurfaceWorkflowConfig(
                backend="cp2k", calculation=cp2k_config(), execution=psteros.ExecutionPolicy(queue="other")))
    qe = psteros.QeCalculationConfig("pw@cluster", "sssp", {"CONTROL": {}, "SYSTEM": {}, "ELECTRONS": {}})
    with pytest.raises(TypeError, match="requires CP2K recipes"):
        psteros.build_cp2k_relax_static_workgraph(
            structures, psteros.SurfaceWorkflowConfig(backend="qe", calculation=qe, execution=execution), static)
    with pytest.raises(TypeError, match="requires CP2K recipes"):
        psteros.build_cp2k_surface_workgraph(
            structures, psteros.SurfaceWorkflowConfig(backend="qe", calculation=qe, execution=execution))


def _cp2k_code(tmp_path):
    from aiida import orm

    computer = orm.Computer(
        label=f"local-{tmp_path.name}", hostname="localhost", transport_type="core.local",
        scheduler_type="core.direct", workdir=str(tmp_path),
    ).store()
    orm.InstalledCode(computer=computer, filepath_executable="/bin/true", label="cp2k",
                      default_calc_job_plugin="cp2k").store()
    return computer.label


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_cp2k_relax_static_workgraph(tmp_path) -> None:
    pytest.importorskip("aiida_cp2k")
    computer = _cp2k_code(tmp_path)
    basis = tmp_path / "BASIS_MOLOPT_UZH"
    basis.write_text("Ru DZVP-MOLOPT-PBE-GTH-q16\n")
    execution = psteros.ExecutionPolicy(computer=computer, queue="par128", max_concurrent_jobs=3,
                                        extra_options={"import_sys_environment": False})
    code = f"cp2k@{computer}"
    o2 = psteros.CalculationOverride(kpoints_distance=100.0,
                                     parameters={"FORCE_EVAL": {"DFT": {"UKS": True, "MULTIPLICITY": 3}}})
    relax = psteros.SurfaceWorkflowConfig(
        backend="cp2k", name="ruo2", execution=execution,
        calculation=cp2k_config(code, run_type="GEO_OPT", input_files={"BASIS_MOLOPT_UZH": str(basis)}),
        role_overrides={
            "bulk": psteros.CalculationOverride(parameters={"GLOBAL": {"RUN_TYPE": "CELL_OPT"},
                                                            "MOTION": {"CELL_OPT": {"MAX_FORCE": 1.94e-4}}}),
            "slab": psteros.CalculationOverride(fixed_sites=(0, 1)),
            "o2": o2,
        },
    )
    static = psteros.SurfaceWorkflowConfig(
        backend="cp2k", name="ruo2_static", execution=execution, calculation=cp2k_config(code),
        role_overrides={"o2": o2},
    )
    structures = _structures()
    graph = psteros.build_cp2k_relax_static_workgraph(structures, relax, static)
    names = {task.name for task in graph.tasks}
    assert {f"{label}_{stage}_cp2k" for label in structures for stage in ("relax", "static")} <= names

    bulk = graph.tasks["bulk_relax_cp2k"].inputs.cp2k
    assert bulk.parameters.value.get_dict()["GLOBAL"]["RUN_TYPE"] == "CELL_OPT"
    assert bulk.kpoints.value.get_kpoints_mesh()[0] == list(cp2k_kpoints_mesh(structures["bulk"].lattice.matrix, 0.2))
    assert bulk.metadata.options.parser_name.value == "cp2k_advanced_parser"
    assert bulk.metadata.options.queue_name.value == "par128"
    assert bulk.metadata.options.import_sys_environment.value is False
    assert bulk.file.BASIS_MOLOPT_UZH.value.filename == "BASIS_MOLOPT_UZH"
    assert graph.tasks["bulk_relax_cp2k"].inputs.handler_overrides.value.get_dict() == {
        "restart_incomplete_calculation": {"enabled": True}}

    slab_relax = graph.tasks["slab_relax_cp2k"].inputs.cp2k.parameters.value.get_dict()
    assert slab_relax["MOTION"]["CONSTRAINT"]["FIXED_ATOMS"]["LIST"] == "1 2"
    slab_static = graph.tasks["slab_static_cp2k"].inputs.cp2k
    static_parameters = slab_static.parameters.value.get_dict()
    assert static_parameters["GLOBAL"]["RUN_TYPE"] == "ENERGY_FORCE" and "MOTION" not in static_parameters
    assert slab_static.settings.value.get_dict() == {"additional_retrieve_list": [FORCES_FILE, DOS_FILE]}
    assert {kind["_"] for kind in static_parameters["FORCE_EVAL"]["SUBSYS"]["KIND"]} == {"Sn", "O"}

    o2_static = graph.tasks["o2_static_cp2k"].inputs.cp2k
    assert o2_static.kpoints.value is None  # Gamma only for the molecule
    slab_mesh = slab_static.kpoints.value.get_kpoints_mesh()[0]
    assert slab_mesh[2] == 1  # one point across the vacuum of the slab
    assert o2_static.parameters.value.get_dict()["FORCE_EVAL"]["DFT"]["MULTIPLICITY"] == 3

    links = {(link.from_task.name, link.to_task.name) for link in graph.links}
    assert ("slab_relax_cp2k", "slab_static_cp2k") in links
    outputs = set(graph.outputs._get_keys()) if hasattr(graph.outputs, "_get_keys") else set(graph.outputs)
    assert {"slab_relaxed_structure", "slab_relax_parameters", "slab_static_parameters",
            "slab_relax_retrieved", "slab_static_retrieved"} <= outputs


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_cp2k_surface_workgraph_one_static_per_label(tmp_path) -> None:
    pytest.importorskip("aiida_cp2k")
    computer = _cp2k_code(tmp_path)
    config = psteros.SurfaceWorkflowConfig(
        backend="cp2k", name="labels", calculation=cp2k_config(f"cp2k@{computer}", kpoints_spacing=None),
        execution=psteros.ExecutionPolicy(computer=computer),
    )
    graph = psteros.build_cp2k_surface_workgraph(_structures(), config)
    assert {task.name for task in graph.tasks} >= {"bulk_cp2k", "slab_cp2k", "o2_cp2k"}
    assert graph.tasks["bulk_cp2k"].inputs.cp2k.kpoints.value is None
    outputs = set(graph.outputs._get_keys()) if hasattr(graph.outputs, "_get_keys") else set(graph.outputs)
    assert {"bulk_parameters", "bulk_structure", "bulk_retrieved"} <= outputs

"""Tier-1/2 tests for the VASP backend, the central psteros engine."""

from __future__ import annotations

import pytest

import psteros
from psteros.backends.vasp import vasp_incar, vasp_parameters, vasp_positions_dof, vasp_potential_mapping
from psteros.polar_workflow import read_vasp_results


def vasp_config(**options) -> psteros.VaspCalculationConfig:
    values = dict(code_label="vasp@cluster", incar={"ENCUT": 520, "IBRION": 2, "NSW": 100, "ISIF": 2},
                  potential_mapping={"Sn": "Sn_d"})
    values.update(options)
    return psteros.VaspCalculationConfig(**values)


def test_parameters_use_the_aiida_vasp_incar_namespace() -> None:
    config = vasp_config()
    override = psteros.CalculationOverride(parameters={"INCAR": {"ISIF": 3}})
    assert vasp_incar(config, override) == {"ENCUT": 520, "IBRION": 2, "NSW": 100, "ISIF": 3}
    assert vasp_parameters(config, override) == {"incar": {"ENCUT": 520, "IBRION": 2, "NSW": 100, "ISIF": 3}}
    assert config.incar["ISIF"] == 2
    with pytest.raises(ValueError, match="'INCAR' namespace"):
        vasp_incar(config, psteros.CalculationOverride(parameters={"SYSTEM": {"ecutwfc": 40}}))


def test_fixed_sites_become_selective_dynamics() -> None:
    override = psteros.CalculationOverride(fixed_sites=(2, 0, 2))
    assert override.fixed_sites == (0, 2)
    parameters = vasp_parameters(vasp_config(), override, number_of_sites=3)
    assert parameters["dynamics"] == {"positions_dof": [[False] * 3, [True] * 3, [False] * 3]}
    assert vasp_positions_dof(2, [1]) == [[True] * 3, [False] * 3]
    with pytest.raises(ValueError, match="number of sites"):
        vasp_parameters(vasp_config(), override)
    with pytest.raises(ValueError, match="outside"):
        vasp_parameters(vasp_config(), override, number_of_sites=2)
    with pytest.raises(ValueError, match="non-negative integer"):
        psteros.CalculationOverride(fixed_sites=(-1,))


def test_potential_mapping_defaults_to_the_element() -> None:
    config = vasp_config()
    assert vasp_potential_mapping(config, {"Sn": "Sn", "O": "O"}) == {"Sn": "Sn_d", "O": "O"}
    with pytest.raises(ValueError, match="H0p75"):
        vasp_potential_mapping(config, {"Ga": "Ga", "H0p75": "H"})
    assert vasp_potential_mapping(vasp_config(potential_mapping={"H0p75": "H.75"}), {"H0p75": "H"}) == {
        "H0p75": "H.75"}


def test_config_validates_the_incar() -> None:
    with pytest.raises(ValueError, match="flat INCAR"):
        vasp_config(incar={"incar": {"ENCUT": 520}})
    with pytest.raises(ValueError, match="max_iterations"):
        vasp_config(max_iterations=0)
    assert vasp_config(max_iterations=3).max_iterations == 3


def test_central_sites_of_a_slab() -> None:
    slab, _ = psteros.sno2_110_slab(termination="o", triple_layers=5)
    fixed = psteros.central_sites(slab, half_width=1.5)
    import numpy as np

    normal = np.cross(slab.lattice.matrix[0], slab.lattice.matrix[1])
    normal /= np.linalg.norm(normal)
    heights = np.array([site.coords @ normal for site in slab])
    middle = (heights.max() + heights.min()) / 2
    assert fixed and all(abs(heights[i] - middle) < 1.5 for i in fixed)
    assert len(fixed) < len(slab)
    with pytest.raises(ValueError, match="half_width"):
        psteros.central_sites(slab, half_width=0)


def test_relax_static_recipes_are_checked_before_the_graph() -> None:
    execution = psteros.ExecutionPolicy(computer="cluster", queue="normal")
    relax = psteros.SurfaceWorkflowConfig(backend="vasp", calculation=vasp_config(), execution=execution)
    static = psteros.SurfaceWorkflowConfig(backend="vasp", calculation=vasp_config(incar={"ENCUT": 520, "NSW": 0}),
                                           execution=execution)
    structures = {"bulk": psteros.rutile_sno2_bulk()}
    with pytest.raises(ValueError, match="static INCAR must not relax"):
        psteros.build_relax_static_workgraph(structures, relax, relax)
    with pytest.raises(ValueError, match="relaxation INCAR needs"):
        psteros.build_relax_static_workgraph(structures, static, static)
    single_point = psteros.CalculationOverride(parameters={"INCAR": {"NSW": 0}})
    with pytest.raises(ValueError, match="bulk: the relaxation INCAR"):
        psteros.build_relax_static_workgraph(
            structures, psteros.SurfaceWorkflowConfig(
                backend="vasp", calculation=vasp_config(), execution=execution,
                role_overrides={"bulk": single_point}),
            static,
        )
    qe = psteros.QeCalculationConfig("pw@cluster", "sssp", {"CONTROL": {}, "SYSTEM": {}, "ELECTRONS": {}})
    with pytest.raises(ValueError, match="same backend"):
        psteros.build_relax_static_workgraph(
            structures, relax, psteros.SurfaceWorkflowConfig(backend="qe", calculation=qe, execution=execution))
    with pytest.raises(TypeError, match="requires QE recipes"):
        psteros.build_qe_relax_static_workgraph(structures, relax, static)


def test_read_vasp_results_from_a_relax_static_graph() -> None:
    class Graph:
        outputs = {
            "a_relax_misc": {"total_energies": {"energy_extrapolated": -12.0}},
            "a_static_misc": {"total_energies": {"energy_extrapolated": -12.5}},
            "a_relaxed_structure": "relaxed-a",
        }

    energies, structures = read_vasp_results(Graph(), ["a"])
    assert energies == {"a": -12.5} and structures == {"a": "relaxed-a"}
    with pytest.raises(ValueError, match="no misc output"):
        read_vasp_results(Graph(), ["b"])


def test_study_overrides_by_stage() -> None:
    pytest.importorskip("pymatgen")
    from pymatgen.core import Lattice, Structure

    mgo = Structure.from_spacegroup("Fm-3m", Lattice.cubic(4.21), ["Mg", "O"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    mg = Structure.from_spacegroup("P6_3/mmc", Lattice.hexagonal(3.21, 5.21), ["Mg"], [[1 / 3, 2 / 3, 0.25]])
    study = psteros.ChargeNeutralSurfaceStudy(mgo, [(1, 0, 0)], {"Mg": mg, "O": psteros.triplet_o2_cell()},
                                              min_slab_thickness=8)
    static = study.vasp_overrides("static")
    assert dict(static["bulk"].parameters["INCAR"]) == {}
    assert dict(static["ref_Mg"].parameters["INCAR"]) == {}
    assert dict(static["ref_O"].parameters["INCAR"]) == {"ISPIN": 2, "MAGMOM": [1.0, 1.0]}
    assert static["ref_O"].kpoints_distance == 10.0
    fixed = study.vasp_overrides(extra={"MgO_100_term_0": psteros.CalculationOverride(fixed_sites=(0,))})
    assert fixed["MgO_100_term_0"].fixed_sites == (0,)
    with pytest.raises(ValueError, match="stage"):
        study.vasp_overrides("scf")


def _vasp_code(tmp_path):
    from aiida import orm

    computer = orm.Computer(
        label=f"local-{tmp_path.name}", hostname="localhost", transport_type="core.local",
        scheduler_type="core.direct", workdir=str(tmp_path),
    ).store()
    orm.InstalledCode(computer=computer, filepath_executable="/bin/true", label="vasp",
                      default_calc_job_plugin="vasp.vasp").store()
    return computer.label


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_vasp_relax_static_workgraph(tmp_path) -> None:
    pytest.importorskip("aiida_vasp")
    computer = _vasp_code(tmp_path)
    execution = psteros.ExecutionPolicy(computer=computer, queue="debug", max_concurrent_jobs=1)
    relax = psteros.SurfaceWorkflowConfig(
        backend="vasp", name="sno2", execution=execution,
        calculation=vasp_config(code_label=f"vasp@{computer}", max_iterations=2),
        role_overrides={
            "bulk": psteros.CalculationOverride(parameters={"INCAR": {"ISIF": 3}}),
            "slab": psteros.CalculationOverride(fixed_sites=(0, 1)),
        },
    )
    static = psteros.SurfaceWorkflowConfig(
        backend="vasp", name="sno2_static", execution=execution,
        calculation=vasp_config(code_label=f"vasp@{computer}", incar={"ENCUT": 520, "NSW": 0, "IBRION": -1}),
    )
    slab, _ = psteros.sno2_110_slab(termination="o", triple_layers=3)
    graph = psteros.build_relax_static_workgraph(
        {"bulk": psteros.rutile_sno2_bulk(), "slab": slab}, relax, static)
    names = {task.name for task in graph.tasks}
    assert {"bulk_relax_vasp", "bulk_static_vasp", "slab_relax_vasp", "slab_static_vasp"} <= names
    bulk = graph.tasks["bulk_relax_vasp"].inputs
    assert bulk.parameters.value.get_dict() == {"incar": {"ENCUT": 520, "IBRION": 2, "NSW": 100, "ISIF": 3}}
    assert bulk.potential_mapping.value.get_dict() == {"Sn": "Sn_d", "O": "O"}
    assert bulk.max_iterations.value.value == 2
    dynamics = graph.tasks["slab_relax_vasp"].inputs.parameters.value.get_dict()["dynamics"]["positions_dof"]
    assert dynamics[:3] == [[False] * 3, [False] * 3, [True] * 3] and len(dynamics) == len(slab)
    assert graph.tasks["slab_static_vasp"].inputs.parameters.value.get_dict() == {
        "incar": {"ENCUT": 520, "NSW": 0, "IBRION": -1}}
    links = {(link.from_task.name, link.to_task.name) for link in graph.links}
    assert ("slab_relax_vasp", "slab_static_vasp") in links
    outputs = set(graph.outputs._get_keys()) if hasattr(graph.outputs, "_get_keys") else set(graph.outputs)
    assert {"slab_relaxed_structure", "slab_static_misc", "slab_relax_misc"} <= outputs


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_qe_fixed_sites_become_fixed_coords(tmp_path) -> None:
    pytest.importorskip("aiida_quantumespresso")
    from aiida import orm

    from tests.unit.test_surface_study import _pseudo_family

    computer = orm.Computer(
        label=f"local-{tmp_path.name}", hostname="localhost", transport_type="core.local",
        scheduler_type="core.direct", workdir=str(tmp_path),
    ).store()
    orm.InstalledCode(computer=computer, filepath_executable="/bin/true", label="pw",
                      default_calc_job_plugin="quantumespresso.pw").store()
    slab, _ = psteros.sno2_110_slab(termination="o", triple_layers=3)
    config = psteros.SurfaceWorkflowConfig(
        backend="qe",
        calculation=psteros.QeCalculationConfig(
            f"pw@{computer.label}", _pseudo_family(tmp_path, ["Sn", "O"]),
            {"CONTROL": {"calculation": "relax"}, "SYSTEM": {"ecutwfc": 40.0}, "ELECTRONS": {}}),
        execution=psteros.ExecutionPolicy(computer=computer.label, queue="debug"),
        role_overrides={"slab": psteros.CalculationOverride(fixed_sites=(1,))},
    )
    graph = psteros.build_surface_workgraph({"slab": slab}, config)
    flags = graph.tasks["slab_qe"].inputs.pw.settings.value.get_dict()["FIXED_COORDS"]
    assert flags[1] == [True] * 3 and flags[0] == [False] * 3 and len(flags) == len(slab)

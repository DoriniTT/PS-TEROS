"""The VASP adapter hands aiida-vasp the k-point spacing in the unit aiida-vasp expects."""

from __future__ import annotations

import math

import pytest

from psteros.backends.vasp import aiida_vasp_kpoints_spacing

pytest.importorskip("pymatgen")


def _mesh(lattice, spacing_for_aiida_vasp):
    """Mesh aiida-vasp builds: AiiDA's ``set_kpoints_mesh_from_density(spacing * 2*pi)``."""
    distance = spacing_for_aiida_vasp * 2.0 * math.pi
    return [math.ceil(norm / distance) for norm in lattice.reciprocal_lattice.abc]


def test_conversion_removes_the_two_pi_that_aiida_vasp_adds():
    assert aiida_vasp_kpoints_spacing(0.2) == pytest.approx(0.2 / (2.0 * math.pi))
    assert aiida_vasp_kpoints_spacing(0.2) * 2.0 * math.pi == pytest.approx(0.2)


def test_documented_spacing_gives_a_real_mesh_for_gaas():
    from pymatgen.core import Lattice

    gaas = Lattice.cubic(5.75)
    assert _mesh(gaas, 0.2) == [1, 1, 1]  # what an unconverted 0.2 gave: Gamma only
    assert _mesh(gaas, aiida_vasp_kpoints_spacing(0.2)) == [6, 6, 6]


def test_molecule_spacing_stays_gamma_only():
    from pymatgen.core import Lattice

    assert _mesh(Lattice.cubic(12.0), aiida_vasp_kpoints_spacing(10.0)) == [1, 1, 1]


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_the_workgraph_input_is_converted(tmp_path):
    pytest.importorskip("aiida_vasp")
    from aiida import orm
    from aiida_workgraph import WorkGraph
    from pymatgen.core import Lattice, Structure

    import psteros
    from psteros.backends.vasp import add_vasp_task

    computer = orm.Computer(
        label=f"local-{tmp_path.name}", hostname="localhost", transport_type="core.local",
        scheduler_type="core.direct", workdir=str(tmp_path),
    ).store()
    orm.InstalledCode(computer=computer, filepath_executable="/bin/true", label="vasp",
                      default_calc_job_plugin="vasp.vasp").store()
    config = psteros.VaspCalculationConfig(
        code_label=f"vasp@{computer.label}", incar={"ENCUT": 400}, potential_mapping={"Ga": "Ga_d", "As": "As"},
        kpoints_spacing=0.2,
    )
    gaas = Structure.from_spacegroup("F-43m", Lattice.cubic(5.75), ["Ga", "As"], [[0, 0, 0], [0.25, 0.25, 0.25]])
    graph = WorkGraph("kpoints")
    task = add_vasp_task(graph, label="bulk", structure=gaas, config=config,
                         execution=psteros.ExecutionPolicy(computer=computer.label, queue="debug"))
    assert task.inputs.kpoints_spacing.value.value == pytest.approx(0.2 / (2.0 * math.pi))
    assert config.kpoints_spacing == 0.2  # the recipe itself keeps the documented unit

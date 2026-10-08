"""Tier-1 tests for the psteros 2.0 public surface API."""

from __future__ import annotations

import pytest

import psteros


def qe_parameters() -> dict:
    return {
        "CONTROL": {"calculation": "relax"},
        "SYSTEM": {"ecutwfc": 80.0, "ecutrho": 640.0},
        "ELECTRONS": {"conv_thr": 1.0e-8},
    }


def test_public_api_is_vasp_first_and_versioned() -> None:
    assert psteros.__version__ == "2.0.0"
    assert "VASP" in psteros.__doc__.splitlines()[2]
    assert callable(psteros.build_surface_workgraph)
    assert callable(psteros.build_relax_static_workgraph)
    assert callable(psteros.build_qe_relax_static_workgraph)
    assert psteros.VaspCalculationConfig is not None
    assert psteros.QeCalculationConfig is not None


def test_execution_policy_assumes_no_computer() -> None:
    policy = psteros.ExecutionPolicy()
    assert policy.computer is None and policy.queue is None
    assert policy.max_concurrent_jobs == 1
    # Only what was asked for reaches AiiDA: no queue, wall time or scheduler directive.
    assert policy.scheduler_options() == {"resources": {"num_machines": 1}, "withmpi": True}


def test_execution_policy_passes_the_users_scheduler_settings() -> None:
    policy = psteros.ExecutionPolicy(
        computer="my-cluster", queue="my-queue", account="my-project",
        resources={"num_machines": 2, "num_mpiprocs_per_machine": 32},
        max_wallclock_seconds=3600, max_concurrent_jobs=4,
        custom_scheduler_commands="#SBATCH --constraint=cpu", prepend_text="module load vasp",
    )
    assert policy.scheduler_options() == {
        "resources": {"num_machines": 2, "num_mpiprocs_per_machine": 32},
        "withmpi": True,
        "max_wallclock_seconds": 3600,
        "queue_name": "my-queue",
        "account": "my-project",
        "custom_scheduler_commands": "#SBATCH --constraint=cpu",
        "prepend_text": "module load vasp",
    }
    assert psteros.ExecutionPolicy(max_concurrent_jobs=None).max_concurrent_jobs is None


def test_execution_policy_rejects_invalid_values() -> None:
    with pytest.raises(ValueError, match="max_concurrent_jobs"):
        psteros.ExecutionPolicy(max_concurrent_jobs=0)
    with pytest.raises(ValueError, match="max_wallclock_seconds"):
        psteros.ExecutionPolicy(max_wallclock_seconds=0)
    with pytest.raises(ValueError, match="queue must not be empty"):
        psteros.ExecutionPolicy(queue="")
    with pytest.raises(ValueError, match="resources"):
        psteros.ExecutionPolicy(resources={})


def test_recipe_requires_an_execution_policy() -> None:
    qe = psteros.QeCalculationConfig("pw@my-cluster", "sssp", qe_parameters())
    with pytest.raises(TypeError):
        psteros.SurfaceWorkflowConfig(backend="qe", calculation=qe)


def test_qe_config_requires_complete_namelists() -> None:
    config = psteros.QeCalculationConfig(
        code_label="pw@my-cluster",
        pseudo_family="SSSP/1.3/PBE/precision",
        parameters=qe_parameters(),
    )
    assert config.kpoints_distance == 0.20
    assert config.max_iterations == 1
    with pytest.raises(ValueError, match="missing namelists"):
        psteros.QeCalculationConfig(
            code_label="pw@my-cluster",
            pseudo_family="sssp",
            parameters={"CONTROL": {}, "SYSTEM": {}},
        )


def test_qe_relaxation_controls_must_use_control_namelist() -> None:
    parameters = qe_parameters()
    parameters["IONS"] = {"ion_dynamics": "bfgs", "forc_conv_thr": 1.0e-3}
    with pytest.raises(ValueError, match="must be in CONTROL"):
        psteros.QeCalculationConfig("pw@my-cluster", "sssp", parameters)


def test_backend_type_must_match_recipe() -> None:
    qe = psteros.QeCalculationConfig("pw@my-cluster", "sssp", qe_parameters())
    execution = psteros.ExecutionPolicy()
    with pytest.raises(TypeError, match="requires VaspCalculationConfig"):
        psteros.SurfaceWorkflowConfig(backend="vasp", calculation=qe, execution=execution)
    with pytest.raises(ValueError, match="backend must be 'qe' or 'vasp'"):
        psteros.SurfaceWorkflowConfig(backend="Qe", calculation=qe, execution=execution)


def test_sno2_slab_rejects_unknown_termination_and_layer_count() -> None:
    with pytest.raises(ValueError, match="termination must be"):
        psteros.sno2_110_slab(termination="O-rich")
    with pytest.raises(TypeError, match="triple_layers must be an integer"):
        psteros.sno2_110_slab(triple_layers=5.0)


def test_calculation_override_keeps_qe_settings_per_structure() -> None:
    fixed = [[False, False, False], [True, True, True]]
    override = psteros.CalculationOverride(settings={"FIXED_COORDS": fixed})
    assert override.settings["FIXED_COORDS"] == fixed


def test_qe_fixed_coordinate_flags_use_aiida_qe_fixed_semantics() -> None:
    flags = psteros.qe_fixed_coordinate_flags(4, {1, 3})
    assert flags == [
        [False, False, False],
        [True, True, True],
        [False, False, False],
        [True, True, True],
    ]
    with pytest.raises(ValueError, match="outside"):
        psteros.qe_fixed_coordinate_flags(2, {2})


def test_sno2_published_formula_series_is_symmetric() -> None:
    expected = {
        "o": "Sn18 O36",
        "sno": "Sn18 O34",
        "sn2o": "Sn18 O32",
    }
    for termination, formula in expected.items():
        slab, identity = psteros.sno2_110_slab(
            termination=termination, triple_layers=9, vacuum_angstrom=20.0
        )
        assert slab.composition.formula == formula
        assert identity.miller_index == (1, 1, 0)
        assert identity.termination == termination


def test_sno2_slab_is_convertible_to_plain_structure() -> None:
    slab, _ = psteros.sno2_110_slab(termination="o", triple_layers=5)
    from pymatgen.core import Structure

    plain = Structure(
        lattice=slab.lattice,
        species=slab.species_and_occu,
        coords=slab.frac_coords,
        coords_are_cartesian=False,
    )
    assert plain.composition == slab.composition


def test_reference_starting_structures_have_expected_compositions() -> None:
    assert psteros.alpha_sn_bulk().composition.reduced_formula == "Sn"
    assert psteros.litharge_sno_bulk().composition.reduced_formula == "SnO"
    oxygen = psteros.triplet_o2_cell()
    assert oxygen.composition.formula == "O2"
    assert oxygen.get_distance(0, 1) == pytest.approx(1.208)


def test_reference_crystal_starting_models_have_expected_symmetry() -> None:
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    assert SpacegroupAnalyzer(psteros.rutile_sno2_bulk()).get_space_group_number() == 136
    assert SpacegroupAnalyzer(psteros.alpha_sn_bulk()).get_space_group_number() == 227
    assert SpacegroupAnalyzer(psteros.litharge_sno_bulk()).get_space_group_number() == 129


def test_oxide_surface_energy_and_stability_are_deterministic() -> None:
    point = psteros.surface_energy_oxide_equilibrium(
        slab_energy_ev=-100.0,
        n_metal=10,
        n_oxygen=20,
        bulk_formula_energy_ev=-9.0,
        oxygen_reference_energy_ev=-10.0,
        delta_mu_oxygen_ev=-1.0,
        surface_area_angstrom2=20.0,
    )
    assert point.gamma_ev_per_angstrom2 == pytest.approx(-0.25)
    assert point.gamma_j_per_m2 == pytest.approx(-4.005441585)
    stable = psteros.stable_termination(
        {
            "o": [psteros.SurfaceEnergyPoint(-2.0, 0.10), psteros.SurfaceEnergyPoint(0.0, 0.20)],
            "sno": [psteros.SurfaceEnergyPoint(-2.0, 0.20), psteros.SurfaceEnergyPoint(0.0, 0.10)],
        }
    )
    assert stable == [(-2.0, "o", 0.10), (0.0, "sno", 0.10)]

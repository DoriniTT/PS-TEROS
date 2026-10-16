"""Precision levels of the VASP recipe: what they set, how an explicit value wins, and that recipes without them are unchanged."""

from __future__ import annotations

import math

import pytest

import psteros
from psteros.backends.vasp import aiida_vasp_kpoints_spacing, vasp_incar, vasp_parameters, vasp_potential_mapping
from psteros.precision import COMMON_INCAR, PRECISION_LEVELS, level_incar, precision_level

pytest.importorskip("aiida_vasp")  # the POTCAR and magnetic-moment tables are read from its protocol file

TWO_PI = 2.0 * math.pi


def recipe(**options) -> psteros.VaspCalculationConfig:
    values = {"code_label": "vasp@cluster", "incar": {"IBRION": -1, "NSW": 0}, "precision": "moderate"}
    values.update(options)
    return psteros.VaspCalculationConfig(**values)


# ------------------------------------------------------------------ recipes without a level do not change


def test_a_recipe_without_precision_is_what_it_was() -> None:
    incar = {"ENCUT": 520, "IBRION": 2, "NSW": 100, "ISIF": 2}
    config = psteros.VaspCalculationConfig("vasp@cluster", incar, potential_mapping={"Sn": "Sn_d"})  # positional, as before
    assert config.precision is None and config.incar == incar
    assert config.kpoints_spacing == 0.20 and config.effective_kpoints_spacing == 0.20
    assert vasp_incar(config) == incar
    assert vasp_incar(config, site_symbols=["Sn", "O"]) == incar  # nothing is added, not even the spin
    assert vasp_parameters(config, site_symbols=["Sn", "O"]) == {"incar": {key.lower(): value for key, value in incar.items()}}
    assert vasp_potential_mapping(config, {"Sn": "Sn", "O": "O", "Sr": "Sr"}) == {"Sn": "Sn_d", "O": "O", "Sr": "Sr"}
    with pytest.raises(ValueError, match="INCAR must not be empty"):
        psteros.VaspCalculationConfig("vasp@cluster", {})
    with pytest.raises(ValueError, match="flat INCAR"):
        psteros.VaspCalculationConfig("vasp@cluster", {"ENCUT": 520, "nested": {"b": 2}})


def test_a_spacing_written_by_the_user_is_kept_even_when_it_equals_the_default() -> None:
    assert recipe(precision=None, kpoints_spacing=0.2).effective_kpoints_spacing == 0.2
    assert recipe(kpoints_spacing=0.2).effective_kpoints_spacing == 0.2  # explicit 0.2 beats the level
    assert recipe(kpoints_spacing=0.5, precision="fast").effective_kpoints_spacing == 0.5


# ------------------------------------------------------------------ the levels


def test_the_three_levels_are_ordered_from_cheap_to_precise() -> None:
    assert list(PRECISION_LEVELS) == ["fast", "moderate", "precise"]
    fast, moderate, precise = PRECISION_LEVELS.values()
    assert fast.encut_ev < moderate.encut_ev < precise.encut_ev
    assert fast.kpoints_spacing > moderate.kpoints_spacing > precise.kpoints_spacing  # a smaller spacing is a denser mesh
    assert fast.ediff_ev > moderate.ediff_ev > precise.ediff_ev
    assert fast.force_cutoff_ev_per_angstrom > moderate.force_cutoff_ev_per_angstrom > precise.force_cutoff_ev_per_angstrom
    assert moderate.kpoints_spacing == pytest.approx(0.05 * TWO_PI)
    assert moderate.ediffg_ev_per_angstrom == -0.03


def test_an_unknown_level_names_the_known_ones() -> None:
    with pytest.raises(ValueError, match=r"\['fast', 'moderate', 'precise'\].*'accurate'"):
        precision_level("accurate")
    with pytest.raises(ValueError, match="precision must be one of"):
        psteros.VaspCalculationConfig("vasp@cluster", {"NSW": 0}, precision="stringent")


def test_leaving_the_incar_out_selects_the_moderate_level() -> None:
    config = psteros.VaspCalculationConfig("vasp@cluster")
    assert config.precision == "moderate" and config.incar == {}
    incar = vasp_incar(config)
    assert incar["ENCUT"] == 550.0 and incar["EDIFF"] == 1e-6 and incar["GGA"] == "PE" and incar["ISMEAR"] == 0
    assert "EDIFFG" not in incar and "IBRION" not in incar  # nothing says that it relaxes
    assert config.effective_kpoints_spacing == pytest.approx(0.05 * TWO_PI)
    assert psteros.VaspCalculationConfig("vasp@cluster", precision="fast").incar == {}  # a level alone is complete too


def test_every_level_sets_its_cutoff_convergence_and_spacing() -> None:
    for name, level in PRECISION_LEVELS.items():
        config = recipe(precision=name)
        incar = vasp_incar(config)
        assert incar["ENCUT"] == level.encut_ev and incar["EDIFF"] == level.ediff_ev
        assert config.effective_kpoints_spacing == pytest.approx(level.kpoints_spacing)
        assert {key: incar[key] for key in COMMON_INCAR} == dict(COMMON_INCAR)


# ------------------------------------------------------------------ what the user writes wins


def test_explicit_tags_win_whatever_their_case() -> None:
    config = recipe(incar={"encut": 600, "Sigma": 0.1, "IBRION": -1, "NSW": 0}, precision="fast")
    incar = vasp_incar(config)
    assert incar["encut"] == 600 and "ENCUT" not in incar
    assert incar["Sigma"] == 0.1 and "SIGMA" not in incar
    assert incar["EDIFF"] == 1e-5  # what the user did not write comes from the level
    override = psteros.CalculationOverride(parameters={"INCAR": {"ENCUT": 700, "EDIFF": 1e-8}})
    stored = vasp_parameters(config, override)["incar"]  # what aiida-vasp receives: the override is the last word
    assert stored["encut"] == 700 and stored["ediff"] == 1e-8 and stored["sigma"] == 0.1


def test_the_stored_parameters_are_lower_case_for_aiida_vasp() -> None:
    parameters = vasp_parameters(recipe(precision="precise"))
    assert set(parameters) == {"incar"} and all(key == key.lower() for key in parameters["incar"])
    assert parameters["incar"]["encut"] == 700.0 and parameters["incar"]["ediff"] == 1e-7


def test_the_force_criterion_is_added_only_to_a_relaxation() -> None:
    relax = recipe(incar={"IBRION": 2, "NSW": 100, "ISIF": 3}, precision="precise")
    assert vasp_incar(relax)["EDIFFG"] == -0.01
    assert "EDIFFG" not in vasp_incar(recipe(incar={"IBRION": -1, "NSW": 0}))  # a static calculation
    assert "EDIFFG" not in vasp_incar(recipe(incar={"IBRION": 5, "NSW": 1}))  # finite differences are no relaxation
    assert vasp_incar(recipe(incar={"IBRION": 2, "NSW": 10, "EDIFFG": -0.001}))["EDIFFG"] == -0.001
    turned_into_relaxation = psteros.CalculationOverride(parameters={"INCAR": {"IBRION": 2, "NSW": 50}})
    assert vasp_incar(recipe(), turned_into_relaxation)["EDIFFG"] == -0.03  # decided on the final INCAR of the structure


# ------------------------------------------------------------------ spin and POTCARs


def test_spin_is_off_unless_the_structure_holds_a_magnetic_element() -> None:
    config = recipe()
    assert vasp_incar(config, site_symbols=["Sn", "Sn", "O", "O"])["ISPIN"] == 1
    assert "MAGMOM" not in vasp_incar(config, site_symbols=["Sn", "O"])
    magnetic = vasp_incar(config, site_symbols=["Fe", "O", "O"])
    assert magnetic["ISPIN"] == 2 and magnetic["MAGMOM"] == [5.0, 1.0, 1.0]  # the moments of the aiida-vasp protocol
    assert "ISPIN" not in vasp_incar(config)  # without the structure the level cannot tell
    given = vasp_incar(recipe(incar={"IBRION": -1, "NSW": 0, "ispin": 1}), site_symbols=["Fe", "O"])
    assert given["ispin"] == 1 and "ISPIN" not in given and "MAGMOM" not in given
    assert "ISPIN" not in vasp_incar(recipe(incar={"IBRION": -1, "NSW": 0, "MAGMOM": [3.0]}), site_symbols=["Fe"])
    assert level_incar("fast", {}, site_symbols=["Co"])["MAGMOM"] == [5.0]


def test_potcars_come_from_the_aiida_vasp_protocol_for_the_elements_the_recipe_does_not_name() -> None:
    config = recipe(potential_mapping={"O": "O_s"})
    assert vasp_potential_mapping(config, {"Sr": "Sr", "Ti": "Ti", "O": "O"}) == {"Sr": "Sr_sv", "Ti": "Ti_pv", "O": "O_s"}
    assert vasp_potential_mapping(recipe(), {"Sn": "Sn", "O": "O"}) == {"Sn": "Sn_d", "O": "O"}
    with pytest.raises(ValueError, match="H0p75"):
        vasp_potential_mapping(recipe(), {"H0p75": "H"})  # a custom kind still needs its own entry
    assert vasp_potential_mapping(recipe(precision=None), {"Sr": "Sr"}) == {"Sr": "Sr"}  # no level: as before


# ------------------------------------------------------------------ the link to aiida-vasp


def test_moderate_is_the_balanced_protocol_of_aiida_vasp() -> None:
    pytest.importorskip("aiida")
    from aiida_vasp.workchains.v2.vasp import VaspWorkChain

    protocol = VaspWorkChain.get_protocol_inputs("balanced")
    tags = {key.upper(): value for key, value in protocol["parameters"]["incar"].items()}
    moderate = PRECISION_LEVELS["moderate"]
    assert moderate.encut_ev == tags["ENCUT"]
    assert moderate.kpoints_spacing_aiida_vasp == protocol["kpoints_spacing"]
    assert moderate.kpoints_spacing == pytest.approx(0.05 * TWO_PI)
    for key in ("ISMEAR", "SIGMA", "LREAL", "LASPH", "NELM"):
        assert COMMON_INCAR[key] == tags[key]
    assert str(COMMON_INCAR["PREC"]).lower() == str(tags["PREC"]).lower() and str(COMMON_INCAR["ALGO"]).lower() == tags["ALGO"]
    assert moderate.ediff_ev == protocol["meta_parameters"]["ediff_per_atom"]  # per atom there, per cell here
    from psteros.precision import protocol_tables

    assert protocol_tables()["potential_mapping"] == dict(protocol["potential_mapping"])
    assert protocol_tables()["magmom_mapping"] == dict(protocol["magmom_mapping"])
    assert VaspWorkChain.get_protocol_inputs("moderate")["kpoints_spacing"] == protocol["kpoints_spacing"]  # alias
    assert tags["GGA"] == "ps" and COMMON_INCAR["GGA"] == "PE"  # the one deliberate difference with the functional


def test_every_tag_that_a_level_writes_is_accepted_by_aiida_vasp() -> None:
    from aiida_vasp.assistant.parameters import ParametersMassage

    for name in PRECISION_LEVELS:
        for incar in ({"IBRION": 2, "NSW": 100, "ISIF": 3}, {"IBRION": -1, "NSW": 0}):
            parameters = vasp_parameters(recipe(incar=incar, precision=name), site_symbols=["Fe", "O"])
            assert ParametersMassage(parameters).parameters.incar["encut"] == PRECISION_LEVELS[name].encut_ev


# ------------------------------------------------------------------ graphs


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_the_graph_builders_use_the_level(tmp_path) -> None:
    from aiida import orm

    computer = orm.Computer(
        label=f"local-{tmp_path.name}", hostname="localhost", transport_type="core.local",
        scheduler_type="core.direct", workdir=str(tmp_path),
    ).store()
    orm.InstalledCode(computer=computer, filepath_executable="/bin/true", label="vasp", default_calc_job_plugin="vasp.vasp").store()
    code = f"vasp@{computer.label}"
    execution = psteros.ExecutionPolicy(computer=computer.label)
    relax = psteros.SurfaceWorkflowConfig(
        backend="vasp", name="level", execution=execution,
        calculation=psteros.VaspCalculationConfig(code, {"IBRION": 2, "NSW": 100, "ISIF": 3}, precision="moderate"),
    )
    assert psteros.VaspCalculationConfig(code, {"IBRION": -1, "NSW": 0}).precision is None  # an INCAR alone selects no level
    both = psteros.SurfaceWorkflowConfig(
        backend="vasp", name="level_static", execution=execution,
        calculation=psteros.VaspCalculationConfig(code, {"IBRION": -1, "NSW": 0}, precision="moderate"),
    )
    graph = psteros.build_relax_static_workgraph({"bulk": psteros.rutile_sno2_bulk()}, relax, both)
    relax_task, static_task = graph.tasks["bulk_relax_vasp"], graph.tasks["bulk_static_vasp"]
    relax_incar = relax_task.inputs.parameters.value.get_dict()["incar"]
    static_incar = static_task.inputs.parameters.value.get_dict()["incar"]
    assert relax_incar["encut"] == 550.0 and relax_incar["ediffg"] == -0.03 and relax_incar["isif"] == 3
    assert static_incar["encut"] == 550.0 and "ediffg" not in static_incar and static_incar["ispin"] == 1
    assert relax_task.inputs.kpoints_spacing.value.value == pytest.approx(aiida_vasp_kpoints_spacing(0.05 * TWO_PI))
    assert relax_task.inputs.potential_mapping.value.get_dict() == {"Sn": "Sn_d", "O": "O"}


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_the_reference_blocks_use_the_level(tmp_path) -> None:
    from aiida import orm

    computer = orm.Computer(
        label=f"local-{tmp_path.name}", hostname="localhost", transport_type="core.local",
        scheduler_type="core.direct", workdir=str(tmp_path),
    ).store()
    orm.InstalledCode(computer=computer, filepath_executable="/bin/true", label="vasp", default_calc_job_plugin="vasp.vasp").store()
    config = psteros.SurfaceWorkflowConfig(
        backend="vasp", name="level", execution=psteros.ExecutionPolicy(computer=computer.label),
        calculation=psteros.VaspCalculationConfig(f"vasp@{computer.label}", precision="fast"),
    )
    references = {"sno2": psteros.ReferenceSystem(psteros.rutile_sno2_bulk(), "solid")}
    graph = psteros.build_vasp_reference_workgraph(references, config, blocks=(psteros.Relax(incar={"ibrion": 2, "nsw": 100}), psteros.Static()))
    relax = graph.tasks["sno2_relax_vasp"].inputs
    static = graph.tasks["sno2_static_vasp"].inputs
    assert relax.parameters.value.get_dict()["incar"]["encut"] == 400.0
    assert relax.parameters.value.get_dict()["incar"]["ediffg"] == -0.05
    assert "ediffg" not in static.parameters.value.get_dict()["incar"]
    assert relax.kpoints_spacing.value.value == pytest.approx(0.07)
    assert relax.potential_mapping.value.get_dict() == {"Sn": "Sn_d", "O": "O"}


# ------------------------------------------------------------------ the example


def _campaign():
    import importlib.util
    import sys
    from pathlib import Path

    directory = Path(__file__).resolve().parents[2] / "examples" / "vasp_surface_phase_diagram"
    sys.path.insert(0, str(directory))
    try:
        spec = importlib.util.spec_from_file_location("precision_example_campaign", directory / "campaign.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(directory))
    return module


def test_the_example_campaign_takes_a_level_or_keeps_its_own_settings() -> None:
    import types

    campaign = _campaign()
    args = types.SimpleNamespace(
        code="vasp@cluster", potential_family="PBE", computer=None, queue=None, ranks=4, walltime=3600, precision=None
    )
    relax, static = campaign.recipes(args, "sno2")  # without --precision: the settings written in the script
    assert relax.calculation.precision is None and relax.calculation.incar["ENCUT"] == 400
    assert relax.calculation.kpoints_spacing == 0.3 and relax.calculation.incar["EDIFFG"] == -0.02
    args.precision = "fast"
    relax, static = campaign.recipes(args, "sno2")
    assert relax.calculation.precision == "fast"
    assert dict(relax.calculation.incar) == {"IBRION": 2, "NSW": 100, "ISIF": 2}
    assert dict(static.calculation.incar) == {"IBRION": -1, "NSW": 0}
    assert relax.calculation.effective_kpoints_spacing == pytest.approx(0.07 * TWO_PI)

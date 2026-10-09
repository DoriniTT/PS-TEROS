"""Tier-1 tests for the campaign API: slab systems and public names (no AiiDA needed)."""

from __future__ import annotations

import dataclasses
import inspect

import pytest

import psteros

CAMPAIGN_NAMES = (
    "SlabSystem",
    "build_vasp_campaign_workgraph",
    "campaign_results",
    "campaign_terminations",
)

EXISTING_NAMES = (
    "ReferenceSystem",
    "build_vasp_reference_workgraph",
    "build_surface_workgraph",
    "build_qe_relax_static_workgraph",
    "reference_results",
    "reference_thermochemistry",
    "Relax",
    "Static",
    "Vibrations",
)


def bulk_structure():
    # SlabSystem keeps whatever structure it is given; a bulk cell is a realistic stand-in.
    return psteros.alpha_sn_bulk()


def test_slab_system_is_a_frozen_dataclass_for_a_solid() -> None:
    system = psteros.SlabSystem(bulk_structure())
    assert dataclasses.is_dataclass(system)
    assert {field.name for field in dataclasses.fields(psteros.SlabSystem)} == {
        "structure",
        "override",
        "block_overrides",
    }
    assert system.phase == "solid"  # a class constant, not a field
    assert system.override is None
    assert dict(system.block_overrides) == {}
    with pytest.raises(dataclasses.FrozenInstanceError):
        system.override = psteros.CalculationOverride()


def test_slab_system_keeps_its_structure_and_does_not_share_its_default_mapping() -> None:
    structure = bulk_structure()
    first, second = psteros.SlabSystem(structure), psteros.SlabSystem(structure)
    assert first.structure is structure
    assert first.block_overrides is not second.block_overrides


def test_slab_overrides_must_be_calculation_overrides() -> None:
    structure = bulk_structure()
    with pytest.raises(TypeError, match="CalculationOverride"):
        psteros.SlabSystem(structure, override={"INCAR": {"ENCUT": 600}})
    with pytest.raises(TypeError, match="CalculationOverride"):
        psteros.SlabSystem(structure, block_overrides={"static": {"INCAR": {"ENCUT": 600}}})


def test_vasp_slab_overrides_take_incar_tags_only() -> None:
    structure = bulk_structure()
    other = psteros.CalculationOverride(parameters={"SYSTEM": {"ecutwfc": 40}})
    with pytest.raises(ValueError, match="INCAR"):
        psteros.SlabSystem(structure, override=other)
    mixed = psteros.CalculationOverride(
        parameters={"INCAR": {"ENCUT": 600}, "DYNAMICS": {"positions_dof": [[True, True, True]]}}
    )
    with pytest.raises(ValueError, match="INCAR"):
        psteros.SlabSystem(structure, override=mixed)


def test_slab_override_errors_name_the_block_they_come_from() -> None:
    # AGENTS.md: error messages name the block (a SlabSystem has no label of its own).
    structure = bulk_structure()
    with pytest.raises(TypeError) as type_error:
        psteros.SlabSystem(structure, block_overrides={"static": {"INCAR": {"ENCUT": 600}}})
    assert "static" in str(type_error.value)
    other = psteros.CalculationOverride(parameters={"SYSTEM": {"ecutwfc": 40}})
    with pytest.raises(ValueError) as namespace_error:
        psteros.SlabSystem(structure, block_overrides={"static": other})
    assert "static" in str(namespace_error.value)


def test_slab_overrides_accept_incar_tags_and_the_other_override_fields() -> None:
    structure = bulk_structure()
    full = psteros.CalculationOverride(
        parameters={"INCAR": {"ENCUT": 600}},
        kpoints_distance=0.3,
        settings={"parser_settings": {"add_dos": True}},
        metadata={"max_wallclock_seconds": 600},
    )
    system = psteros.SlabSystem(
        structure,
        override=full,
        block_overrides={"static": psteros.CalculationOverride(kpoints_distance=0.4)},
    )
    assert system.override == full
    assert system.block_overrides["static"].kpoints_distance == 0.4
    empty = psteros.SlabSystem(structure, override=psteros.CalculationOverride())
    assert empty.override == psteros.CalculationOverride()


def test_build_vasp_campaign_workgraph_takes_keyword_only_options_after_the_systems() -> None:
    parameters = inspect.signature(psteros.build_vasp_campaign_workgraph).parameters
    assert list(parameters) == ["references", "slabs", "config", "reference_blocks", "slab_blocks", "submit"]
    for name in ("reference_blocks", "slab_blocks", "submit"):
        assert parameters[name].kind is inspect.Parameter.KEYWORD_ONLY, name
    assert parameters["submit"].default is False
    assert parameters["reference_blocks"].default == (psteros.Relax(), psteros.Static(), psteros.Vibrations())
    assert parameters["slab_blocks"].default == (psteros.Relax(), psteros.Static())


def test_new_campaign_names_are_exported() -> None:
    for name in CAMPAIGN_NAMES:
        assert name in psteros.__all__, name
        assert callable(getattr(psteros, name)), name


def test_existing_public_names_are_still_exported() -> None:
    for name in EXISTING_NAMES:
        assert name in psteros.__all__, name
        assert callable(getattr(psteros, name)), name
    assert len(psteros.__all__) == len(set(psteros.__all__))
    for name in psteros.__all__:
        assert hasattr(psteros, name), name

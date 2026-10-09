"""Tier-2 tests: the VASP campaign WorkGraph, references and slabs in one graph.

Graphs are built on the throwaway profile of ``conftest.py`` and never submitted.
The reader tests store bare nodes that carry the extras a submitted graph sets.
"""

from __future__ import annotations

import dataclasses
import math
import re

import pytest

pytest.importorskip("aiida_vasp")

import psteros  # noqa: E402

pytestmark = pytest.mark.requires_aiida

INCAR = {"ENCUT": 520, "EDIFF": 1e-7, "ISMEAR": 0, "SIGMA": 0.05, "IBRION": 2, "NSW": 100, "EDIFFG": -0.005}
BUILT_IN_TASKS = {"graph_ctx", "graph_inputs", "graph_outputs"}
BAD_LABELS = ["110_o", "a__b", "slab-o", "_x", "slab_"]
DEFAULT_REFERENCE_BLOCKS = [("relax", "relax"), ("static", "static"), ("vibrations", "vibrations")]
DEFAULT_SLAB_BLOCKS = [("relax", "relax"), ("static", "static")]
OUTPUTS_BY_KIND = {
    "relax": {"misc", "retrieved", "remote", "energy", "structure"},
    "static": {"misc", "retrieved", "remote", "energy"},
    "vibrations": {"misc", "retrieved", "remote", "frequencies"},
}
REFERENCE_BLOCKS_EXTRA = [
    {"name": "relax", "kind": "relax"},
    {"name": "static", "kind": "static"},
    {"name": "vibrations", "kind": "vibrations"},
]
SLAB_BLOCKS_EXTRA = [{"name": "relax", "kind": "relax"}, {"name": "static", "kind": "static"}]


@pytest.fixture(scope="module")
def code_label() -> str:
    from aiida import orm

    label = "vasp-campaign-test"
    try:
        computer = orm.load_computer("localhost-test")
    except Exception:
        computer = orm.Computer(
            label="localhost-test", hostname="localhost", transport_type="core.local", scheduler_type="core.direct"
        ).store()
    try:
        orm.load_code(f"{label}@localhost-test")
    except Exception:
        orm.InstalledCode(
            computer=computer, filepath_executable="/bin/true", default_calc_job_plugin="vasp.vasp", label=label
        ).store()
    return f"{label}@localhost-test"


@pytest.fixture(scope="module")
def slab_o():
    return psteros.sno2_110_slab(termination="o", triple_layers=3, vacuum_angstrom=15.0)[0]


def recipe(code_label: str, incar=None) -> psteros.SurfaceWorkflowConfig:
    return psteros.SurfaceWorkflowConfig(
        backend="vasp",
        calculation=psteros.VaspCalculationConfig(
            code_label=code_label,
            incar=INCAR if incar is None else incar,
            potential_mapping={"Sn": "Sn_d", "O": "O"},
            kpoints_spacing=0.25,
        ),
        execution=psteros.ExecutionPolicy(computer="localhost-test", queue="debug"),
        name="sno2",
    )


def references() -> dict:
    cell_relax = psteros.CalculationOverride(parameters={"INCAR": {"ISIF": 3}})
    return {
        "o2": psteros.ReferenceSystem(
            psteros.triplet_o2_cell(),
            "gas",
            symmetry_number=2,
            spin=1.0,
            override=psteros.CalculationOverride(
                parameters={"INCAR": {"ISPIN": 2, "NUPDOWN": 2}},
                kpoints_distance=5.0,
                metadata={"max_wallclock_seconds": 600},
            ),
        ),
        "sno2": psteros.ReferenceSystem(
            psteros.rutile_sno2_bulk(), "solid", supercell=(2, 2, 3), block_overrides={"relax": cell_relax}
        ),
    }


def slab_systems(slab) -> dict:
    return {"slab_o": psteros.SlabSystem(slab)}


def task_names(workgraph) -> set[str]:
    return {task.name for task in workgraph.tasks} - BUILT_IN_TASKS


def links(workgraph) -> set[tuple[str, str, str, str]]:
    return {
        (link["from_task"], link["from_socket"], link["to_task"], link["to_socket"])
        for link in (link.to_dict() for link in workgraph.links)
    }


def incar(workgraph, task: str) -> dict:
    return workgraph.tasks[task].inputs.parameters.value.get_dict()["incar"]


def task_inputs(workgraph, task: str) -> dict:
    """Inputs of a ``vasp.v2.vasp`` task that the campaign and the reference builder must share."""

    inputs = workgraph.tasks[task].inputs
    return {
        "code": inputs.code.value.full_label,
        "parameters": inputs.parameters.value.get_dict(),
        "kpoints_spacing": inputs.kpoints_spacing.value.value,
        "potential_family": inputs.potential_family.value.value,
        "potential_mapping": inputs.potential_mapping.value.get_dict(),
        "options": inputs.options.value.get_dict(),
        "settings": inputs.settings.value.get_dict(),
        "clean_workdir": inputs.clean_workdir.value.value,
    }


def output_names(namespace, prefix: str = "") -> set[str]:
    """Dotted names of the leaf outputs below a WorkGraph output namespace."""

    names: set[str] = set()
    for name, socket in namespace._sockets.items():
        if hasattr(socket, "_sockets"):  # a nested namespace, not a leaf output
            names |= output_names(socket, f"{prefix}{name}.")
        else:
            names.add(f"{prefix}{name}")
    return names


def nested_outputs(group: str, label: str, blocks: list[tuple[str, str]]) -> set[str]:
    """Expected dotted output names of one label; ``blocks`` holds (block name, kind) pairs."""

    return {f"{group}.{label}.{block}.{socket}" for block, kind in blocks for socket in OUTPUTS_BY_KIND[kind]}


def flat_reference_outputs(label: str) -> set[str]:
    return {
        f"{label}_relax_misc", f"{label}_relax_retrieved", f"{label}_relax_structure", f"{label}_relax_energy",
        f"{label}_static_misc", f"{label}_static_retrieved", f"{label}_static_energy",
        f"{label}_vibrations_misc", f"{label}_vibrations_retrieved", f"{label}_vibrations_frequencies",
    }


# =============================================================================
# Graph structure: task names, outputs, links
# =============================================================================


def test_default_blocks_name_every_task_by_label_and_block(code_label, slab_o) -> None:
    workgraph = psteros.build_vasp_campaign_workgraph(references(), slab_systems(slab_o), recipe(code_label))
    assert task_names(workgraph) == {
        # o2, a gas: relaxation, static and vibrations (no supercell)
        "o2_relax_vasp", "o2_relax_energy",
        "o2_static_vasp", "o2_static_energy",
        "o2_vibrations_vasp", "o2_vibrations_frequencies",
        # sno2, a solid: its vibrations are displaced in a 2 x 2 x 3 supercell
        "sno2_relax_vasp", "sno2_relax_energy",
        "sno2_static_vasp", "sno2_static_energy",
        "sno2_vibrations_supercell", "sno2_vibrations_vasp", "sno2_vibrations_frequencies",
        # slab_o: the default slab blocks are relaxation and static, no vibrations
        "slab_o_relax_vasp", "slab_o_relax_energy",
        "slab_o_static_vasp", "slab_o_static_energy",
    }


def test_a_solid_in_its_own_cell_has_no_supercell_task(code_label) -> None:
    bulk = {"sn_bulk": psteros.ReferenceSystem(psteros.alpha_sn_bulk(), "solid")}
    workgraph = psteros.build_vasp_campaign_workgraph(bulk, {}, recipe(code_label))
    assert task_names(workgraph) == {
        "sn_bulk_relax_vasp", "sn_bulk_relax_energy",
        "sn_bulk_static_vasp", "sn_bulk_static_energy",
        "sn_bulk_vibrations_vasp", "sn_bulk_vibrations_frequencies",
    }


def test_outputs_are_nested_by_structure_and_block(code_label, slab_o) -> None:
    workgraph = psteros.build_vasp_campaign_workgraph(references(), slab_systems(slab_o), recipe(code_label))
    outputs = output_names(workgraph.outputs)
    assert outputs == (
        nested_outputs("references", "o2", DEFAULT_REFERENCE_BLOCKS)
        | nested_outputs("references", "sno2", DEFAULT_REFERENCE_BLOCKS)
        | nested_outputs("slabs", "slab_o", DEFAULT_SLAB_BLOCKS)
    )
    for name in (
        "references.o2.static.energy",
        "references.o2.relax.structure",
        "references.sno2.vibrations.frequencies",
        "references.o2.relax.remote",
        "slabs.slab_o.relax.structure",
        "slabs.slab_o.static.energy",
        "slabs.slab_o.static.remote",
    ):
        assert name in outputs, name
    # A static block has no structure, and vibrations have no energy.
    assert "references.o2.static.structure" not in outputs
    assert "slabs.slab_o.static.structure" not in outputs
    assert "references.sno2.vibrations.energy" not in outputs


def test_each_block_takes_the_structure_it_is_linked_to(code_label, slab_o) -> None:
    workgraph = psteros.build_vasp_campaign_workgraph(references(), slab_systems(slab_o), recipe(code_label))
    graph_links = links(workgraph)
    # references: the gas is displaced in its own cell, the solid in a supercell
    assert ("o2_relax_vasp", "structure", "o2_static_vasp", "structure") in graph_links
    assert ("o2_relax_vasp", "structure", "o2_vibrations_vasp", "structure") in graph_links
    assert ("sno2_relax_vasp", "structure", "sno2_static_vasp", "structure") in graph_links
    assert ("sno2_relax_vasp", "structure", "sno2_vibrations_supercell", "structure") in graph_links
    assert ("sno2_vibrations_supercell", "result", "sno2_vibrations_vasp", "structure") in graph_links
    # slabs: static on the relaxed slab
    assert ("slab_o_relax_vasp", "structure", "slab_o_static_vasp", "structure") in graph_links
    # energies and frequencies read the outputs of their own calculation
    assert ("o2_static_vasp", "misc", "o2_static_energy", "misc") in graph_links
    assert ("o2_static_vasp", "retrieved", "o2_static_energy", "retrieved") in graph_links
    assert ("o2_vibrations_vasp", "retrieved", "o2_vibrations_frequencies", "retrieved") in graph_links
    assert ("slab_o_static_vasp", "misc", "slab_o_static_energy", "misc") in graph_links


def test_slab_blocks_take_their_names_and_structure_from_the_blocks(code_label, slab_o) -> None:
    blocks = (
        psteros.Relax(name="coarse"),
        psteros.Relax(name="fine"),
        psteros.Static(name="scf", structure_from="coarse"),
    )
    workgraph = psteros.build_vasp_campaign_workgraph({}, slab_systems(slab_o), recipe(code_label), slab_blocks=blocks)
    assert task_names(workgraph) == {
        "slab_o_coarse_vasp", "slab_o_coarse_energy",
        "slab_o_fine_vasp", "slab_o_fine_energy",
        "slab_o_scf_vasp", "slab_o_scf_energy",
    }
    graph_links = links(workgraph)
    assert ("slab_o_coarse_vasp", "structure", "slab_o_fine_vasp", "structure") in graph_links
    assert ("slab_o_coarse_vasp", "structure", "slab_o_scf_vasp", "structure") in graph_links
    assert ("slab_o_fine_vasp", "structure", "slab_o_scf_vasp", "structure") not in graph_links
    assert output_names(workgraph.outputs) == nested_outputs(
        "slabs", "slab_o", [("coarse", "relax"), ("fine", "relax"), ("scf", "static")]
    )


def test_either_group_may_be_left_empty(code_label, slab_o) -> None:
    slabs_only = psteros.build_vasp_campaign_workgraph({}, slab_systems(slab_o), recipe(code_label))
    assert task_names(slabs_only) == {
        "slab_o_relax_vasp", "slab_o_relax_energy", "slab_o_static_vasp", "slab_o_static_energy",
    }
    references_only = psteros.build_vasp_campaign_workgraph(references(), {}, recipe(code_label))
    assert output_names(references_only.outputs) == (
        nested_outputs("references", "o2", DEFAULT_REFERENCE_BLOCKS)
        | nested_outputs("references", "sno2", DEFAULT_REFERENCE_BLOCKS)
    )


def test_graph_is_named_after_the_recipe_and_runs_one_job_at_a_time(code_label, slab_o) -> None:
    workgraph = psteros.build_vasp_campaign_workgraph(references(), slab_systems(slab_o), recipe(code_label))
    assert workgraph.name == "sno2_campaign"
    assert workgraph.max_number_jobs == 1
    renamed = dataclasses.replace(recipe(code_label), name="psteros_test")
    named = psteros.build_vasp_campaign_workgraph({}, slab_systems(slab_o), renamed)
    assert named.name == "psteros_test_campaign"


# =============================================================================
# INCAR, k-points, options and settings
# =============================================================================


def test_reference_tasks_match_the_reference_builder_task_by_task(code_label, slab_o) -> None:
    reference_graph = psteros.build_vasp_reference_workgraph(references(), recipe(code_label))
    campaign = psteros.build_vasp_campaign_workgraph(references(), slab_systems(slab_o), recipe(code_label))
    vasp_tasks = {name for name in task_names(reference_graph) if name.endswith("_vasp")}
    assert len(vasp_tasks) == 6
    assert vasp_tasks <= task_names(campaign)
    for name in sorted(vasp_tasks):
        assert task_inputs(campaign, name) == task_inputs(reference_graph, name), name
    # The comparison is not vacuous: the reference overrides reach the campaign too.
    assert incar(campaign, "sno2_relax_vasp")["isif"] == 3
    assert "isif" not in incar(campaign, "sno2_static_vasp")
    assert campaign.tasks["sno2_relax_vasp"].inputs.structure.value.get_pymatgen_structure().composition == (
        psteros.rutile_sno2_bulk().composition
    )
    supercell = "sno2_vibrations_supercell"
    assert campaign.tasks[supercell].inputs.size.value.get_list() == [2, 2, 3]
    assert reference_graph.tasks[supercell].inputs.size.value.get_list() == [2, 2, 3]


def test_slab_incar_follows_the_priority_chain(code_label, slab_o) -> None:
    # recipe < block defaults < Block.incar < slab override < block_overrides < required tags
    system = psteros.SlabSystem(
        slab_o,
        override=psteros.CalculationOverride(parameters={"INCAR": {"ENCUT": 700, "NSW": 50}}),
        block_overrides={"static": psteros.CalculationOverride(parameters={"INCAR": {"ENCUT": 800}})},
    )
    workgraph = psteros.build_vasp_campaign_workgraph(
        {},
        {"slab_o": system},
        recipe(code_label),
        slab_blocks=(psteros.Relax(incar={"encut": 600, "ediffg": -0.002}), psteros.Static()),
    )
    relax = incar(workgraph, "slab_o_relax_vasp")
    static = incar(workgraph, "slab_o_static_vasp")
    assert relax["encut"] == 700  # slab override beats Block.incar (600) and the recipe (520)
    assert relax["ediffg"] == -0.002  # Block.incar beats the recipe (-0.005)
    assert relax["nsw"] == 50  # slab override beats the recipe (100)
    assert static["encut"] == 800  # block_overrides beat the slab override
    assert static["ibrion"] == -1  # block defaults beat the recipe (2)
    assert static["nsw"] == 0  # required tags beat the slab override (50)
    assert static["ediffg"] == -0.005  # Block.incar of the relax block does not reach the static block


def test_block_incar_beats_the_block_defaults(code_label, slab_o) -> None:
    workgraph = psteros.build_vasp_campaign_workgraph(
        {}, slab_systems(slab_o), recipe(code_label), slab_blocks=(psteros.Relax(), psteros.Static(incar={"ibrion": 1}))
    )
    assert incar(workgraph, "slab_o_static_vasp")["ibrion"] == 1


def test_slab_kpoints_settings_and_metadata_follow_the_same_layers(code_label, slab_o) -> None:
    system = psteros.SlabSystem(
        slab_o,
        override=psteros.CalculationOverride(
            kpoints_distance=0.3,
            settings={"parser_settings": {"add_dos": True}},
            metadata={"max_wallclock_seconds": 600},
        ),
        block_overrides={"static": psteros.CalculationOverride(kpoints_distance=0.4)},
    )
    workgraph = psteros.build_vasp_campaign_workgraph({}, {"slab_o": system}, recipe(code_label))
    relax = workgraph.tasks["slab_o_relax_vasp"].inputs
    static = workgraph.tasks["slab_o_static_vasp"].inputs
    assert relax.kpoints_spacing.value.value == pytest.approx(0.3 / (2 * math.pi))
    assert static.kpoints_spacing.value.value == pytest.approx(0.4 / (2 * math.pi))
    for inputs in (relax, static):
        assert inputs.options.value.get_dict()["max_wallclock_seconds"] == 600
        assert inputs.settings.value.get_dict()["parser_settings"] == {"add_dos": True}


def test_a_bare_slab_structure_is_wrapped_like_a_slab_system(code_label, slab_o) -> None:
    bare = psteros.build_vasp_campaign_workgraph({}, {"slab_o": slab_o}, recipe(code_label))
    wrapped = psteros.build_vasp_campaign_workgraph({}, slab_systems(slab_o), recipe(code_label))
    assert task_names(bare) == task_names(wrapped)
    assert output_names(bare.outputs) == output_names(wrapped.outputs)
    for name in sorted(task for task in task_names(bare) if task.endswith("_vasp")):
        assert task_inputs(bare, name) == task_inputs(wrapped, name), name
    structure = bare.tasks["slab_o_relax_vasp"].inputs.structure.value
    assert structure.get_pymatgen_structure().composition == slab_o.composition


# =============================================================================
# Validation: every error names the offending label or block
# =============================================================================


def test_a_qe_recipe_is_rejected() -> None:
    qe = psteros.SurfaceWorkflowConfig(
        backend="qe",
        calculation=psteros.QeCalculationConfig("qe@x", "sssp", {"CONTROL": {}, "SYSTEM": {}, "ELECTRONS": {}}),
    )
    with pytest.raises(TypeError, match="(?i)vasp"):
        psteros.build_vasp_campaign_workgraph(references(), {}, qe)


def test_both_groups_empty_is_an_error(code_label) -> None:
    with pytest.raises(ValueError, match=r"(?i)empty|at least one"):
        psteros.build_vasp_campaign_workgraph({}, {}, recipe(code_label))


def test_role_overrides_are_not_used_by_the_campaign(code_label, slab_o) -> None:
    config = dataclasses.replace(
        recipe(code_label),
        role_overrides={"o2": psteros.CalculationOverride(parameters={"INCAR": {"ISPIN": 2}})},
    )
    with pytest.raises(ValueError, match="role_overrides"):
        psteros.build_vasp_campaign_workgraph(references(), slab_systems(slab_o), config)


@pytest.mark.parametrize("label", BAD_LABELS)
def test_a_bad_reference_label_is_named(code_label, label) -> None:
    bulk = {label: psteros.ReferenceSystem(psteros.alpha_sn_bulk(), "solid")}
    with pytest.raises(ValueError, match=re.escape(label)):
        psteros.build_vasp_campaign_workgraph(bulk, {}, recipe(code_label))


@pytest.mark.parametrize("label", BAD_LABELS)
def test_a_bad_slab_label_is_named(code_label, slab_o, label) -> None:
    with pytest.raises(ValueError, match=re.escape(label)):
        psteros.build_vasp_campaign_workgraph({}, {label: slab_o}, recipe(code_label))


def test_every_bad_label_of_a_group_is_listed(code_label) -> None:
    bulk = psteros.ReferenceSystem(psteros.alpha_sn_bulk(), "solid")
    with pytest.raises(ValueError) as info:
        psteros.build_vasp_campaign_workgraph({"110_o": bulk, "a__b": bulk}, {}, recipe(code_label))
    assert "110_o" in str(info.value) and "a__b" in str(info.value)


@pytest.mark.parametrize("label", ["o2", "O2", "slab_o", "a1_b2", "Sn2O", "s_1_x"])
def test_valid_labels_are_accepted(code_label, label) -> None:
    bulk = {label: psteros.ReferenceSystem(psteros.alpha_sn_bulk(), "solid")}
    workgraph = psteros.build_vasp_campaign_workgraph(bulk, {}, recipe(code_label))
    assert f"{label}_relax_vasp" in task_names(workgraph)


def test_a_label_cannot_be_both_a_reference_and_a_slab(code_label, slab_o) -> None:
    o2 = psteros.ReferenceSystem(psteros.triplet_o2_cell(), "gas", symmetry_number=2, spin=1.0)
    with pytest.raises(ValueError, match="o2"):
        psteros.build_vasp_campaign_workgraph({"o2": o2}, {"o2": slab_o}, recipe(code_label))


def test_a_reference_must_be_a_reference_system(code_label) -> None:
    with pytest.raises(TypeError, match="o2"):
        psteros.build_vasp_campaign_workgraph({"o2": psteros.triplet_o2_cell()}, {}, recipe(code_label))


def test_reference_block_overrides_must_name_a_block_of_the_group(code_label) -> None:
    bulk = psteros.ReferenceSystem(
        psteros.alpha_sn_bulk(), "solid", block_overrides={"scf": psteros.CalculationOverride()}
    )
    with pytest.raises(ValueError) as info:
        psteros.build_vasp_campaign_workgraph({"sn_bulk": bulk}, {}, recipe(code_label))
    assert "sn_bulk" in str(info.value) and "scf" in str(info.value)


def test_reference_block_overrides_must_be_in_reference_blocks(code_label) -> None:
    bulk = psteros.ReferenceSystem(
        psteros.alpha_sn_bulk(), "solid", block_overrides={"vibrations": psteros.CalculationOverride()}
    )
    with pytest.raises(ValueError) as info:
        psteros.build_vasp_campaign_workgraph(
            {"sn_bulk": bulk}, {}, recipe(code_label), reference_blocks=(psteros.Relax(), psteros.Static())
        )
    assert "sn_bulk" in str(info.value) and "vibrations" in str(info.value)


def test_slab_block_overrides_must_name_a_block_of_the_group(code_label, slab_o) -> None:
    system = psteros.SlabSystem(slab_o, block_overrides={"scf": psteros.CalculationOverride()})
    with pytest.raises(ValueError) as info:
        psteros.build_vasp_campaign_workgraph({}, {"slab_o": system}, recipe(code_label))
    assert "slab_o" in str(info.value) and "scf" in str(info.value)


@pytest.mark.parametrize("nsw", [0, -1])
def test_a_reference_relaxation_needs_nsw_above_zero(code_label, nsw) -> None:
    bulk = psteros.ReferenceSystem(
        psteros.alpha_sn_bulk(),
        "solid",
        block_overrides={"relax": psteros.CalculationOverride(parameters={"INCAR": {"NSW": nsw}})},
    )
    with pytest.raises(ValueError) as info:
        psteros.build_vasp_campaign_workgraph({"sn_bulk": bulk}, {}, recipe(code_label))
    assert "sn_bulk" in str(info.value) and "relax" in str(info.value)


@pytest.mark.parametrize("nsw", [0, -1])
def test_a_slab_relaxation_needs_nsw_above_zero(code_label, slab_o, nsw) -> None:
    system = psteros.SlabSystem(
        slab_o, block_overrides={"relax": psteros.CalculationOverride(parameters={"INCAR": {"NSW": nsw}})}
    )
    with pytest.raises(ValueError) as info:
        psteros.build_vasp_campaign_workgraph({}, {"slab_o": system}, recipe(code_label))
    assert "slab_o" in str(info.value) and "relax" in str(info.value)


def test_a_recipe_without_nsw_cannot_relax(code_label) -> None:
    # Missing NSW means no ionic steps, the same rule as build_vasp_reference_workgraph.
    no_nsw = {key: value for key, value in INCAR.items() if key != "NSW"}
    bulk = {"sn_bulk": psteros.ReferenceSystem(psteros.alpha_sn_bulk(), "solid")}
    with pytest.raises(ValueError) as info:
        psteros.build_vasp_campaign_workgraph(bulk, {}, recipe(code_label, no_nsw))
    assert "sn_bulk" in str(info.value) and "relax" in str(info.value)


def test_vibrations_are_not_a_slab_block(code_label, slab_o) -> None:
    with pytest.raises(ValueError, match="(?i)vibrations"):
        psteros.build_vasp_campaign_workgraph(
            {}, slab_systems(slab_o), recipe(code_label), slab_blocks=(psteros.Relax(), psteros.Vibrations())
        )


@pytest.mark.parametrize(
    ("option", "blocks", "error", "message"),
    [
        ("reference_blocks", (), ValueError, "at least one block"),
        ("slab_blocks", (psteros.Relax(), psteros.Relax()), ValueError, "duplicate block name"),
        ("slab_blocks", ("relax",), TypeError, "Relax, Static or Vibrations"),
    ],
)
def test_block_lists_are_checked_like_the_reference_builder(code_label, slab_o, option, blocks, error, message) -> None:
    with pytest.raises(error, match=message):
        psteros.build_vasp_campaign_workgraph(
            references(), slab_systems(slab_o), recipe(code_label), **{option: blocks}
        )


# =============================================================================
# Regression: the existing builders keep their flat outputs and task names
# =============================================================================


def test_reference_builder_keeps_its_flat_outputs_and_task_names(code_label) -> None:
    workgraph = psteros.build_vasp_reference_workgraph(references(), recipe(code_label))
    assert workgraph.name == "sno2_references"
    assert workgraph.max_number_jobs == 1
    assert task_names(workgraph) == {
        "o2_relax_vasp", "o2_relax_energy",
        "o2_static_vasp", "o2_static_energy",
        "o2_vibrations_vasp", "o2_vibrations_frequencies",
        "sno2_relax_vasp", "sno2_relax_energy",
        "sno2_static_vasp", "sno2_static_energy",
        "sno2_vibrations_supercell", "sno2_vibrations_vasp", "sno2_vibrations_frequencies",
    }
    assert set(workgraph.outputs._sockets) == flat_reference_outputs("o2") | flat_reference_outputs("sno2")
    graph_links = links(workgraph)
    assert ("o2_relax_vasp", "structure", "o2_static_vasp", "structure") in graph_links
    assert ("sno2_relax_vasp", "structure", "sno2_vibrations_supercell", "structure") in graph_links


def test_surface_builder_keeps_its_flat_outputs(code_label) -> None:
    workgraph = psteros.build_surface_workgraph({"bulk": psteros.rutile_sno2_bulk()}, recipe(code_label))
    assert workgraph.name == "sno2"
    assert task_names(workgraph) == {"bulk_vasp"}
    assert set(workgraph.outputs._sockets) == {"bulk_misc", "bulk_structure", "bulk_retrieved"}


# =============================================================================
# Submission: the extras the readers need (submit is stubbed, nothing runs)
# =============================================================================


@pytest.fixture
def submit_without_daemon(monkeypatch):
    """Replace WorkGraph.submit by one that stores a bare process node and runs nothing."""

    from aiida import orm
    from aiida_workgraph import WorkGraph

    def submit(self):
        self.process = orm.WorkflowNode().store()
        return self.process

    monkeypatch.setattr(WorkGraph, "submit", submit)


def test_building_does_not_submit_by_default(code_label, slab_o, monkeypatch) -> None:
    from aiida_workgraph import WorkGraph

    def refuse(self):
        raise AssertionError("a campaign graph is only submitted with submit=True")

    monkeypatch.setattr(WorkGraph, "submit", refuse)
    workgraph = psteros.build_vasp_campaign_workgraph(references(), slab_systems(slab_o), recipe(code_label))
    assert workgraph.process is None


def test_submit_stores_the_descriptions_the_readers_need(code_label, slab_o, submit_without_daemon) -> None:
    workgraph = psteros.build_vasp_campaign_workgraph(
        references(), slab_systems(slab_o), recipe(code_label), submit=True
    )
    node = workgraph.process
    slab_composition = {element: int(count) for element, count in slab_o.composition.get_el_amt_dict().items()}
    assert node.base.extras.get("psteros_campaign") == {
        "version": 1,
        "references": {
            "blocks": REFERENCE_BLOCKS_EXTRA,
            "labels": ["o2", "sno2"],
            "systems": {
                "o2": {"phase": "gas", "composition": {"O": 2}},
                "sno2": {"phase": "solid", "composition": {"Sn": 2, "O": 4}},
            },
        },
        "slabs": {
            "blocks": SLAB_BLOCKS_EXTRA,
            "labels": ["slab_o"],
            "systems": {"slab_o": {"phase": "solid", "composition": slab_composition}},
        },
    }
    # The reference description is the one reference_results and reference_thermochemistry read.
    assert node.base.extras.get("psteros_references") == {
        "blocks": REFERENCE_BLOCKS_EXTRA,
        "references": {
            "o2": {"phase": "gas", "supercell": [1, 1, 1], "symmetry_number": 2, "spin": 1.0},
            "sno2": {"phase": "solid", "supercell": [2, 2, 3], "symmetry_number": None, "spin": None},
        },
    }
    assert set(psteros.campaign_results(node.pk)["slabs"]) == {"slab_o"}


def test_a_slab_only_submission_stores_no_reference_description(code_label, slab_o, submit_without_daemon) -> None:
    workgraph = psteros.build_vasp_campaign_workgraph({}, slab_systems(slab_o), recipe(code_label), submit=True)
    extras = workgraph.process.base.extras
    assert extras.get("psteros_references", None) is None
    assert extras.get("psteros_campaign")["references"]["labels"] == []
    assert extras.get("psteros_campaign")["slabs"]["labels"] == ["slab_o"]


# =============================================================================
# Readers (no daemon): a bare node carries the extras a submitted graph would set
# =============================================================================

CAMPAIGN_EXTRA = {
    "version": 1,
    "references": {"blocks": REFERENCE_BLOCKS_EXTRA, "labels": ["o2"]},
    "slabs": {"blocks": SLAB_BLOCKS_EXTRA, "labels": ["slab_o"]},
}
REFERENCE_EXTRA = {
    "blocks": REFERENCE_BLOCKS_EXTRA,
    "references": {"o2": {"phase": "gas", "supercell": [1, 1, 1], "symmetry_number": 2, "spin": 1.0}},
}


def unstarted_campaign() -> int:
    """PK of a stored node carrying the campaign extras and no child calculation."""

    from aiida import orm

    node = orm.WorkflowNode().store()
    node.base.extras.set("psteros_campaign", CAMPAIGN_EXTRA)
    node.base.extras.set("psteros_references", REFERENCE_EXTRA)
    return node.pk


def test_readers_reject_a_node_without_the_campaign_extra() -> None:
    from aiida import orm

    pk = orm.Int(1).store().pk
    with pytest.raises(ValueError, match="psteros_campaign"):
        psteros.campaign_results(pk)
    with pytest.raises(ValueError, match="psteros_campaign"):
        psteros.campaign_terminations(pk)


def test_campaign_results_of_an_unstarted_graph_are_complete_but_empty() -> None:
    results = psteros.campaign_results(unstarted_campaign())
    assert set(results) == {"references", "slabs"}
    static = results["slabs"]["slab_o"]["static"]
    assert set(static) == {"kind", "state", "pk", "structure", "misc", "remote", "retrieved", "energy"}
    assert (static["kind"], static["state"], static["pk"], static["energy"]) == ("static", "not started", None, None)
    vibrations = results["references"]["o2"]["vibrations"]
    assert set(vibrations) == {"kind", "state", "pk", "structure", "misc", "remote", "retrieved", "frequencies"}
    assert vibrations["frequencies"] is None


def test_terminations_wait_for_the_slab_block_to_finish() -> None:
    with pytest.raises(ValueError, match="slab_o") as info:
        psteros.campaign_terminations(unstarted_campaign())
    assert "not started" in str(info.value)


def test_terminations_take_the_energy_from_the_last_static_block_by_default() -> None:
    from aiida import orm

    with pytest.raises(ValueError, match="block 'static'"):
        psteros.campaign_terminations(unstarted_campaign())
    with pytest.raises(ValueError, match="block 'relax'"):
        psteros.campaign_terminations(unstarted_campaign(), energy_block="relax")
    with pytest.raises(ValueError, match="energy_block must name a relax or static slab block"):
        psteros.campaign_terminations(unstarted_campaign(), energy_block="vibrations")
    relax_only = orm.WorkflowNode().store()
    relax_only.base.extras.set(
        "psteros_campaign",
        {**CAMPAIGN_EXTRA, "slabs": {"blocks": [{"name": "relax", "kind": "relax"}], "labels": ["slab_o"]}},
    )
    with pytest.raises(ValueError, match="block 'relax'"):
        psteros.campaign_terminations(relax_only.pk)


def test_reference_readers_work_unchanged_on_a_campaign_graph() -> None:
    pk = unstarted_campaign()
    results = psteros.reference_results(pk)
    assert results["o2"]["static"]["state"] == "not started"
    assert results["o2"]["static"]["energy"] is None
    assert results["o2"]["vibrations"]["frequencies"] is None
    with pytest.raises(ValueError, match=r"o2: blocks 'static' \(not started\)"):
        psteros.reference_thermochemistry(pk)


def _finished_child(graph, link_label: str, inputs: dict, outputs: dict):
    """Store a finished child process of ``graph`` with the given input and output nodes."""

    from aiida import orm
    from aiida.common.links import LinkType
    from plumpy import ProcessState

    child = orm.CalcFunctionNode()
    child.set_process_state(ProcessState.FINISHED)
    child.set_exit_status(0)
    child.base.links.add_incoming(graph, LinkType.CALL_CALC, link_label)
    for name, node in inputs.items():
        child.base.links.add_incoming(node, LinkType.INPUT_CALC, name)
    child.store()
    for name, node in outputs.items():
        node.base.links.add_incoming(child, LinkType.CREATE, name)
        node.store()
    return child


def test_readers_return_the_values_of_a_finished_slab() -> None:
    from aiida import orm

    from psteros.backends.qe import as_aiida_structure

    graph = orm.WorkflowNode().store()
    graph.base.extras.set("psteros_campaign", CAMPAIGN_EXTRA)
    initial = as_aiida_structure(psteros.sno2_110_slab(termination="o", triple_layers=3)[0]).store()
    relaxed = as_aiida_structure(psteros.sno2_110_slab(termination="o", triple_layers=3, a=4.83, c=3.24)[0])
    relax = _finished_child(
        graph, "slab_o_relax_vasp", {"structure": initial}, {"structure": relaxed, "misc": orm.Dict({"step": "relax"})}
    )
    static = _finished_child(graph, "slab_o_static_vasp", {"structure": relaxed}, {"misc": orm.Dict({"step": "static"})})
    _finished_child(graph, "slab_o_relax_energy", {}, {"result": orm.Float(-99.0)})
    _finished_child(graph, "slab_o_static_energy", {}, {"result": orm.Float(-100.0)})

    results = psteros.campaign_results(graph.pk)["slabs"]["slab_o"]
    assert (results["relax"]["state"], results["relax"]["pk"], results["relax"]["energy"]) == ("finished", relax.pk, -99.0)
    assert results["relax"]["structure"].uuid == relaxed.uuid
    assert (results["static"]["pk"], results["static"]["energy"], results["static"]["misc"]) == (
        static.pk, -100.0, {"step": "static"},
    )
    assert results["static"]["structure"].uuid == relaxed.uuid  # the input of the static block

    (termination,) = psteros.campaign_terminations(graph.pk)
    expected = psteros.SlabTermination.from_structure("slab_o", -100.0, relaxed)
    assert termination == expected
    (from_relax,) = psteros.campaign_terminations(graph.pk, energy_block="relax")
    assert from_relax.slab_energy_ev == -99.0


# =============================================================================
# Campaign entries and the analysis readers: fake finished children, no daemon
# =============================================================================

SYSTEMS_EXTRA = {
    "version": 1,
    "references": {
        "blocks": REFERENCE_BLOCKS_EXTRA,
        "labels": ["o2", "sno2"],
        "systems": {
            "o2": {"phase": "gas", "composition": {"O": 2}},
            "sno2": {"phase": "solid", "composition": {"Sn": 2, "O": 4}},
        },
    },
    "slabs": {
        "blocks": SLAB_BLOCKS_EXTRA,
        "labels": ["slab_o"],
        "systems": {"slab_o": {"phase": "solid", "composition": {"Sn": 4, "O": 8}}},
    },
}
ELEMENTAL_EXTRA = {
    "version": 1,
    "references": {
        "blocks": REFERENCE_BLOCKS_EXTRA,
        "labels": ["o2", "sn_bulk", "sno2"],
        "systems": {
            "o2": {"phase": "gas", "composition": {"O": 2}},
            "sn_bulk": {"phase": "solid", "composition": {"Sn": 8}},
            "sno2": {"phase": "solid", "composition": {"Sn": 2, "O": 4}},
        },
    },
    "slabs": {
        "blocks": SLAB_BLOCKS_EXTRA,
        "labels": ["slab_o"],
        "systems": {"slab_o": {"phase": "solid", "composition": {"Sn": 4, "O": 8}}},
    },
}


def campaign_graph(extra: dict):
    """A stored graph node that carries a campaign description and no child yet."""

    from aiida import orm

    graph = orm.WorkflowNode().store()
    graph.base.extras.set("psteros_campaign", extra)
    return graph


def as_aiida(structure):
    from psteros.backends.qe import as_aiida_structure

    return as_aiida_structure(structure)


def as_stored(structure):
    return as_aiida(structure).store()


def slab_cell(a_angstrom: float):
    """Sn4O8 in an orthorhombic box a x 4 x 15 A: its exposed face has |a x b| = 4 a A^2."""

    from pymatgen.core import Lattice, Structure

    coords = [
        [0.0, 0.0, 0.10], [0.5, 0.5, 0.10], [0.0, 0.5, 0.20], [0.5, 0.0, 0.20],
        [0.25, 0.25, 0.30], [0.75, 0.75, 0.30], [0.25, 0.75, 0.40], [0.75, 0.25, 0.40],
        [0.0, 0.0, 0.50], [0.5, 0.5, 0.50], [0.0, 0.5, 0.60], [0.5, 0.0, 0.60],
    ]
    return Structure(Lattice.orthorhombic(a_angstrom, 4.0, 15.0), ["Sn"] * 4 + ["O"] * 8, coords)


def finished_campaign():
    """o2 and sno2 with static energies; slab_o relaxed from a = 3.0 to a = 4.0 A, then static."""

    from aiida import orm

    graph = campaign_graph(SYSTEMS_EXTRA)
    o2 = as_stored(psteros.triplet_o2_cell())
    _finished_child(graph, "o2_relax_vasp", {"structure": o2}, {"structure": as_aiida(psteros.triplet_o2_cell())})
    _finished_child(graph, "o2_relax_energy", {}, {"result": orm.Float(-9.1)})
    _finished_child(graph, "o2_static_vasp", {"structure": o2}, {})
    _finished_child(graph, "o2_static_energy", {}, {"result": orm.Float(-9.0)})
    _finished_child(graph, "sno2_static_vasp", {"structure": as_stored(psteros.rutile_sno2_bulk())}, {})
    _finished_child(graph, "sno2_static_energy", {}, {"result": orm.Float(-38.0)})
    relaxed = as_aiida(slab_cell(4.0))
    _finished_child(graph, "slab_o_relax_vasp", {"structure": as_stored(slab_cell(3.0))}, {"structure": relaxed})
    _finished_child(graph, "slab_o_relax_energy", {}, {"result": orm.Float(-99.0)})
    _finished_child(graph, "slab_o_static_vasp", {"structure": relaxed}, {})
    _finished_child(graph, "slab_o_static_energy", {}, {"result": orm.Float(-100.0)})
    return graph


def test_submit_stores_the_phase_and_integer_composition_of_every_system(
    code_label, slab_o, submit_without_daemon
) -> None:
    workgraph = psteros.build_vasp_campaign_workgraph(references(), slab_systems(slab_o), recipe(code_label), submit=True)
    campaign = workgraph.process.base.extras.get("psteros_campaign")
    assert campaign["version"] == 1
    # The composition of each input structure: the 2 x 2 x 3 supercell of sno2 is not counted.
    assert campaign["references"]["systems"] == {
        "o2": {"phase": "gas", "composition": {"O": 2}},
        "sno2": {"phase": "solid", "composition": {"Sn": 2, "O": 4}},
    }
    slab_composition = {element: int(count) for element, count in slab_o.composition.get_el_amt_dict().items()}
    assert slab_composition["O"] == 2 * slab_composition["Sn"]  # the stoichiometric SnO2 slab
    assert campaign["slabs"]["systems"] == {"slab_o": {"phase": "solid", "composition": slab_composition}}
    for group in ("references", "slabs"):
        for label, system in campaign[group]["systems"].items():
            assert all(type(count) is int for count in system["composition"].values()), label


def test_a_slab_only_submission_stores_no_reference_systems(code_label, slab_o, submit_without_daemon) -> None:
    workgraph = psteros.build_vasp_campaign_workgraph({}, slab_systems(slab_o), recipe(code_label), submit=True)
    campaign = workgraph.process.base.extras.get("psteros_campaign")
    assert (campaign["references"]["labels"], campaign["references"]["systems"]) == ([], {})
    assert list(campaign["slabs"]["systems"]) == ["slab_o"]


def test_campaign_entries_read_the_energies_and_compositions_of_each_group() -> None:
    entries = psteros.campaign_entries(finished_campaign().pk)
    # References first, in graph order; the default energy block is the last static one.
    assert [entry.label for entry in entries] == ["o2", "sno2", "slab_o"]
    o2, sno2, slab_entry = entries
    assert (o2.group, o2.phase, dict(o2.composition), o2.energy_ev, o2.block, o2.state) == (
        "references", "gas", {"O": 2}, -9.0, "static", "finished",
    )
    assert (sno2.group, sno2.phase, dict(sno2.composition), sno2.energy_ev) == (
        "references", "solid", {"Sn": 2, "O": 4}, -38.0,
    )
    assert (slab_entry.group, slab_entry.phase, dict(slab_entry.composition), slab_entry.energy_ev) == (
        "slabs", "solid", {"Sn": 4, "O": 8}, -100.0,
    )
    assert slab_entry.block == "static"
    assert o2.surface_area_angstrom2 is None and sno2.surface_area_angstrom2 is None
    # The static block starts from the relaxed slab, a = 4.0 A and b = 4.0 A: |a x b| = 16.0 A^2.
    assert slab_entry.surface_area_angstrom2 == pytest.approx(4.0 * 4.0)


def test_the_energy_blocks_of_campaign_entries_can_be_chosen_by_name() -> None:
    entries = {
        entry.label: entry
        for entry in psteros.campaign_entries(
            finished_campaign().pk, reference_energy_block="relax", slab_energy_block="relax"
        )
    }
    assert (entries["o2"].block, entries["o2"].energy_ev, entries["o2"].state) == ("relax", -9.1, "finished")
    # sno2 has no relaxation yet: no energy, its state, and the composition of its input structure.
    assert (entries["sno2"].block, entries["sno2"].energy_ev, entries["sno2"].state) == ("relax", None, "not started")
    assert dict(entries["sno2"].composition) == {"Sn": 2, "O": 4}
    # The area is that of the relaxed slab (4.0 x 4.0 A^2), not of the initial one (3.0 x 4.0 A^2).
    assert (entries["slab_o"].block, entries["slab_o"].energy_ev) == ("relax", -99.0)
    assert entries["slab_o"].surface_area_angstrom2 == pytest.approx(16.0)


def test_unfinished_blocks_give_no_energy_and_the_composition_of_the_description() -> None:
    entries = psteros.campaign_entries(campaign_graph(SYSTEMS_EXTRA).pk)
    assert [(entry.label, entry.phase, dict(entry.composition), entry.energy_ev, entry.state) for entry in entries] == [
        ("o2", "gas", {"O": 2}, None, "not started"),
        ("sno2", "solid", {"Sn": 2, "O": 4}, None, "not started"),
        ("slab_o", "solid", {"Sn": 4, "O": 8}, None, "not started"),
    ]
    assert entries[2].surface_area_angstrom2 is None  # the static block has no structure yet


def test_a_running_block_gives_its_composition_and_area_but_no_energy() -> None:
    from plumpy import ProcessState

    graph = campaign_graph(SYSTEMS_EXTRA)
    sno2_child = _finished_child(graph, "sno2_static_vasp", {"structure": as_stored(psteros.rutile_sno2_bulk())}, {})
    sno2_child.set_process_state(ProcessState.RUNNING)
    slab_child = _finished_child(graph, "slab_o_static_vasp", {"structure": as_stored(slab_cell(4.0))}, {})
    slab_child.set_process_state(ProcessState.RUNNING)
    entries = {entry.label: entry for entry in psteros.campaign_entries(graph.pk)}
    assert (entries["sno2"].energy_ev, entries["sno2"].state) == (None, "running")
    assert dict(entries["sno2"].composition) == {"Sn": 2, "O": 4}
    assert (entries["slab_o"].energy_ev, entries["slab_o"].state) == (None, "running")
    assert entries["slab_o"].surface_area_angstrom2 == pytest.approx(16.0)


def test_campaign_entries_reject_a_node_without_the_campaign_extra() -> None:
    from aiida import orm

    with pytest.raises(ValueError, match="psteros_campaign"):
        psteros.campaign_entries(orm.Int(1).store().pk)


def test_the_analysis_readers_take_the_pk_of_a_campaign_graph() -> None:
    from aiida import orm

    graph = campaign_graph(ELEMENTAL_EXTRA)
    _finished_child(graph, "o2_static_vasp", {"structure": as_stored(psteros.triplet_o2_cell())}, {})
    _finished_child(graph, "o2_static_energy", {}, {"result": orm.Float(-9.0)})
    _finished_child(graph, "sn_bulk_static_vasp", {"structure": as_stored(psteros.alpha_sn_bulk())}, {})
    _finished_child(graph, "sn_bulk_static_energy", {}, {"result": orm.Float(-32.0)})
    _finished_child(graph, "sno2_static_vasp", {"structure": as_stored(psteros.rutile_sno2_bulk())}, {})
    _finished_child(graph, "sno2_static_energy", {}, {"result": orm.Float(-38.0)})
    _finished_child(graph, "slab_o_static_vasp", {"structure": as_stored(slab_cell(4.0))}, {})
    _finished_child(graph, "slab_o_static_energy", {}, {"result": orm.Float(-40.0)})

    # mu_O = -9.0 / 2 = -4.5 and mu_Sn = -32.0 / 8 = -4.0; the SnO2 cell is not an elemental limit.
    assert psteros.campaign_chemical_potentials(graph.pk) == pytest.approx({"O": -4.5, "Sn": -4.0})
    assert psteros.campaign_chemical_potentials(graph.pk) == psteros.campaign_chemical_potentials(
        psteros.campaign_entries(graph.pk)
    )
    # A slab of two elements needs its chemical potentials given: the elemental limits are not its equilibrium.
    with pytest.raises(ValueError, match="chemical_potentials_ev"):
        psteros.campaign_surface_energies(graph.pk)
    # gamma = (E - N_Sn mu_Sn - N_O mu_O) / (2 A) = (-40.0 - (4 * -4.0 + 8 * -4.5)) / (2 * 16.0) = 12.0 / 32.0
    limits = {"O": -9.0 / 2, "Sn": -32.0 / 8}
    assert psteros.campaign_surface_energies(graph.pk, chemical_potentials_ev=limits) == pytest.approx({"slab_o": 0.375})
    assert psteros.campaign_references(graph.pk, host="sno2") == psteros.BinaryOxideReferences(
        -38.0, {"Sn": 2, "O": 4}, -9.0, -4.0
    )


# =============================================================================
# Readers: a finished VASP task whose energy (or frequencies) task has not delivered
# =============================================================================


def test_a_finished_vasp_task_without_its_energy_task_names_it_in_the_state() -> None:
    graph = campaign_graph(SYSTEMS_EXTRA)
    _finished_child(graph, "slab_o_static_vasp", {"structure": as_stored(slab_cell(4.0))}, {})
    static = psteros.campaign_results(graph.pk)["slabs"]["slab_o"]["static"]
    assert (static["state"], static["energy"]) == ("finished, slab_o_static_energy not started", None)
    (slab_entry,) = [entry for entry in psteros.campaign_entries(graph.pk) if entry.group == "slabs"]
    assert (slab_entry.state, slab_entry.energy_ev) == ("finished, slab_o_static_energy not started", None)
    assert slab_entry.surface_area_angstrom2 == pytest.approx(4.0 * 4.0)  # |a x b| of the input, 4.0 x 4.0 A
    with pytest.raises(ValueError, match="slab_o") as info:
        psteros.campaign_terminations(graph.pk)
    assert "finished, slab_o_static_energy not started" in str(info.value)


def test_an_energy_task_that_is_running_is_named_in_the_state() -> None:
    from plumpy import ProcessState

    graph = campaign_graph(SYSTEMS_EXTRA)
    _finished_child(graph, "slab_o_static_vasp", {"structure": as_stored(slab_cell(4.0))}, {})
    energy = _finished_child(graph, "slab_o_static_energy", {}, {})
    energy.set_process_state(ProcessState.RUNNING)
    static = psteros.campaign_results(graph.pk)["slabs"]["slab_o"]["static"]
    assert (static["state"], static["energy"]) == ("finished, slab_o_static_energy running", None)


def test_the_frequencies_task_of_a_vibrations_block_is_named_in_the_state() -> None:
    graph = campaign_graph(SYSTEMS_EXTRA)
    _finished_child(graph, "o2_vibrations_vasp", {"structure": as_stored(psteros.triplet_o2_cell())}, {})
    vibrations = psteros.campaign_results(graph.pk)["references"]["o2"]["vibrations"]
    assert (vibrations["state"], vibrations["frequencies"]) == ("finished, o2_vibrations_frequencies not started", None)


def test_an_unfinished_vasp_task_keeps_its_plain_state() -> None:
    from plumpy import ProcessState

    graph = campaign_graph(SYSTEMS_EXTRA)
    vasp = _finished_child(graph, "slab_o_static_vasp", {"structure": as_stored(slab_cell(4.0))}, {})
    vasp.set_process_state(ProcessState.RUNNING)
    assert psteros.campaign_results(graph.pk)["slabs"]["slab_o"]["static"]["state"] == "running"


def test_campaign_entries_checks_an_explicit_block_even_for_a_group_without_labels() -> None:
    empty_references = {
        "version": 1,
        "references": {"blocks": REFERENCE_BLOCKS_EXTRA, "labels": [], "systems": {}},
        "slabs": SYSTEMS_EXTRA["slabs"],
    }
    empty_slabs = {
        "version": 1,
        "references": SYSTEMS_EXTRA["references"],
        "slabs": {"blocks": SLAB_BLOCKS_EXTRA, "labels": [], "systems": {}},
    }
    with pytest.raises(ValueError, match="bogus"):
        psteros.campaign_entries(campaign_graph(empty_references).pk, reference_energy_block="bogus")
    with pytest.raises(ValueError, match="bogus"):
        psteros.campaign_entries(campaign_graph(empty_slabs).pk, slab_energy_block="bogus")
    # A valid block of an empty group is accepted, and the empty group contributes no entries.
    entries = psteros.campaign_entries(campaign_graph(empty_references).pk, reference_energy_block="relax")
    assert [entry.label for entry in entries] == ["slab_o"]
    entries = psteros.campaign_entries(campaign_graph(empty_slabs).pk, slab_energy_block="relax")
    assert [entry.label for entry in entries] == ["o2", "sno2"]


# =============================================================================
# Builder: label/block pairs with the same task names, and the structures given
# =============================================================================


def test_label_block_pairs_with_the_same_task_names_are_named(code_label) -> None:
    bulk = psteros.ReferenceSystem(psteros.alpha_sn_bulk(), "solid")
    # sn/o_relax and sn_o/relax both give the task sn_o_relax_vasp.
    with pytest.raises(ValueError) as info:
        psteros.build_vasp_campaign_workgraph(
            {"sn": bulk, "sn_o": bulk},
            {},
            recipe(code_label),
            reference_blocks=(psteros.Relax(name="relax"), psteros.Relax(name="o_relax")),
        )
    assert "sn/o_relax" in str(info.value) and "sn_o/relax" in str(info.value)
    # Without the clash the same labels build: sn_relax_vasp and sn_o_relax_vasp.
    workgraph = psteros.build_vasp_campaign_workgraph(
        {"sn": bulk, "sn_o": bulk}, {}, recipe(code_label), reference_blocks=(psteros.Relax(name="relax"),)
    )
    assert task_names(workgraph) >= {"sn_relax_vasp", "sn_o_relax_vasp"}


def test_label_block_pairs_clash_across_the_reference_and_slab_groups(code_label, slab_o) -> None:
    bulk = psteros.ReferenceSystem(psteros.alpha_sn_bulk(), "solid")
    with pytest.raises(ValueError) as info:
        psteros.build_vasp_campaign_workgraph(
            {"sn": bulk},
            {"sn_o": slab_o},
            recipe(code_label),
            reference_blocks=(psteros.Relax(name="o_relax"),),
            slab_blocks=(psteros.Relax(name="relax"),),
        )
    assert "sn/o_relax" in str(info.value) and "sn_o/relax" in str(info.value)


@pytest.mark.parametrize("structure", ["not a structure", 3.5, [1, 2, 3]])
def test_a_structure_that_is_no_structure_is_a_type_error_naming_its_label(code_label, structure) -> None:
    with pytest.raises(TypeError, match="slab_o"):
        psteros.build_vasp_campaign_workgraph({}, {"slab_o": structure}, recipe(code_label))
    bulk = psteros.ReferenceSystem(structure, "solid")
    with pytest.raises(TypeError, match="sn_bulk"):
        psteros.build_vasp_campaign_workgraph({"sn_bulk": bulk}, {}, recipe(code_label))


def test_a_pk_of_a_node_that_is_no_structure_is_a_type_error_naming_its_label(code_label) -> None:
    from aiida import orm

    not_a_structure = orm.Int(3).store().pk
    with pytest.raises(TypeError, match="slab_o"):
        psteros.build_vasp_campaign_workgraph({}, {"slab_o": not_a_structure}, recipe(code_label))


def test_the_pk_of_a_structure_is_accepted(code_label, slab_o) -> None:
    pk = as_stored(slab_o).pk
    workgraph = psteros.build_vasp_campaign_workgraph({}, {"slab_o": pk}, recipe(code_label))
    assert "slab_o_relax_vasp" in task_names(workgraph)


def test_a_structure_aiida_cannot_convert_keeps_its_own_error_and_names_the_label(code_label) -> None:
    from pymatgen.core import DummySpecies, Lattice, Structure

    dummy = Structure(Lattice.cubic(5.0), [DummySpecies("X")], [[0.0, 0.0, 0.0]])
    with pytest.raises(Exception, match="slab_x: cannot convert the structure") as info:
        psteros.build_vasp_campaign_workgraph({}, {"slab_x": dummy}, recipe(code_label))
    assert info.value.__cause__ is not None

"""aiida-vasp returns relaxed structures with element names only; psteros restores the kinds."""

from __future__ import annotations

import numpy as np
import pytest

import psteros
from psteros.backends.vasp import assign_kinds, has_custom_kinds


def test_kinds_follow_the_atoms_not_the_site_order():
    cell = np.diag([10.0, 10.0, 10.0])
    template_symbols = ["Ga", "As", "H", "H"]
    template = np.array([[0.0, 0.0, 0.0], [0.25, 0.25, 0.25], [0.5, 0.5, 0.5], [0.8, 0.8, 0.8]])
    # reordered by VASP (species grouped differently) and displaced by relaxation, one atom across the cell edge
    relaxed_symbols = ["H", "Ga", "H", "As"]
    relaxed = np.array([[0.801, 0.79, 0.8], [0.995, 0.0, 0.0], [0.5, 0.505, 0.5], [0.26, 0.25, 0.24]])
    assert assign_kinds(template_symbols, template, relaxed_symbols, relaxed, cell) == [3, 0, 2, 1]


def test_an_atom_that_moved_too_far_is_an_error_naming_the_site():
    cell = np.diag([10.0, 10.0, 10.0])
    with pytest.raises(ValueError, match=r"H site 0 .*limit"):
        assign_kinds(["H"], [[0.0, 0.0, 0.0]], ["H"], [[0.5, 0.5, 0.5]], cell)


def test_different_atoms_are_an_error():
    cell = np.eye(3) * 10
    with pytest.raises(ValueError, match="different atoms"):
        assign_kinds(["H", "H"], [[0, 0, 0], [0.5, 0, 0]], ["H", "O"], [[0, 0, 0], [0.5, 0, 0]], cell)


def _structure(names_and_positions, cell=10.0):
    from aiida import orm

    structure = orm.StructureData(cell=np.diag([cell] * 3))
    for symbol, name, position in names_and_positions:
        structure.append_atom(position=position, symbols=symbol, name=name)
    return structure


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_calcfunction_restores_the_pseudo_hydrogen_kinds():
    from psteros.backends.vasp_kinds import restore_vasp_kinds

    template = _structure([("Ga", "Ga", (0, 0, 0)), ("H", "H0p75", (2, 2, 2)), ("H", "H1p25", (6, 6, 6))])
    relaxed = _structure([("H", "H", (6.1, 6.0, 5.9)), ("Ga", "Ga", (0.05, 0, 0)), ("H", "H", (2.0, 2.1, 2.0))])
    assert has_custom_kinds(template) and not has_custom_kinds(relaxed)
    restored = restore_vasp_kinds(relaxed, template)
    assert [site.kind_name for site in restored.sites] == ["H1p25", "Ga", "H0p75"]
    assert np.allclose([site.position for site in restored.sites], [site.position for site in relaxed.sites])


def _code(tmp_path):
    from aiida import orm

    computer = orm.Computer(
        label=f"local-{tmp_path.name}", hostname="localhost", transport_type="core.local",
        scheduler_type="core.direct", workdir=str(tmp_path),
    ).store()
    orm.InstalledCode(computer=computer, filepath_executable="/bin/true", label="vasp",
                      default_calc_job_plugin="vasp.vasp").store()
    return computer.label


def _recipes(computer, mapping):
    execution = psteros.ExecutionPolicy(computer=computer, queue="debug", max_concurrent_jobs=1)

    def recipe(name, incar):
        return psteros.SurfaceWorkflowConfig(
            backend="vasp", name=name, execution=execution,
            calculation=psteros.VaspCalculationConfig(
                code_label=f"vasp@{computer}", incar=incar, potential_mapping=mapping),
        )

    return recipe("kinds", {"ENCUT": 400, "IBRION": 2, "NSW": 50}), recipe("kinds_static", {"ENCUT": 400, "NSW": 0, "IBRION": -1})


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_the_static_step_gets_restored_kinds_only_when_there_are_custom_kinds(tmp_path):
    pytest.importorskip("aiida_vasp")
    computer = _code(tmp_path)
    relax, static = _recipes(computer, {"Ga": "Ga_d", "H0p75": "H.75"})
    with_h = _structure([("Ga", "Ga", (0, 0, 0)), ("H", "H0p75", (2, 2, 2))])
    plain = _structure([("Ga", "Ga", (0, 0, 0)), ("Ga", "Ga", (3, 3, 3))])
    graph = psteros.build_relax_static_workgraph({"pseudo": with_h, "plain": plain}, relax, static)
    names = {task.name for task in graph.tasks}
    assert "pseudo_restore_kinds" in names
    assert "plain_restore_kinds" not in names  # graphs of structures without custom kinds are unchanged
    links = {(link.from_task.name, link.to_task.name) for link in graph.links}
    assert ("pseudo_relax_vasp", "pseudo_restore_kinds") in links
    assert ("pseudo_restore_kinds", "pseudo_static_vasp") in links
    assert ("plain_relax_vasp", "plain_static_vasp") in links


def test_a_restart_finds_the_potcar_of_an_unambiguous_custom_kind():
    from psteros.backends.vasp import restart_potential_aliases

    mapping = {"Ga": "Ga_d", "As": "As", "H1p25": "H1.25"}
    kinds = {"Ga": "Ga", "As": "As", "H1p25": "H"}
    assert restart_potential_aliases(mapping, kinds) == {**mapping, "H": "H1.25"}
    # nothing to alias without custom kinds: the mapping is unchanged
    assert restart_potential_aliases({"Ga": "Ga_d"}, {"Ga": "Ga"}) == {"Ga": "Ga_d"}
    # two custom kinds of one element, or a plain kind of that element: ambiguous, no alias
    two = {**mapping, "H0p75": "H.75"}
    assert "H" not in restart_potential_aliases(two, {**kinds, "H0p75": "H"})
    assert restart_potential_aliases({**mapping, "H": "H"}, {**kinds, "H": "H"}) == {**mapping, "H": "H"}
    # an explicit entry of the user is kept
    assert restart_potential_aliases({**mapping, "H": "H.5"}, kinds)["H"] == "H.5"


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_the_graph_input_carries_the_restart_alias(tmp_path):
    pytest.importorskip("aiida_vasp")
    computer = _code(tmp_path)
    relax, static = _recipes(computer, {"Ga": "Ga_d", "H0p75": "H.75"})
    graph = psteros.build_relax_static_workgraph(
        {"pseudo": _structure([("Ga", "Ga", (0, 0, 0)), ("H", "H0p75", (2, 2, 2))]),
         "plain": _structure([("Ga", "Ga", (0, 0, 0)), ("Ga", "Ga", (3, 3, 3))])}, relax, static)
    assert graph.tasks["pseudo_relax_vasp"].inputs.potential_mapping.value.get_dict() == {
        "Ga": "Ga_d", "H0p75": "H.75", "H": "H.75"}
    # no custom kind in the structure: the mapping is the recipe's, with no alias
    assert graph.tasks["plain_relax_vasp"].inputs.potential_mapping.value.get_dict() == {
        "Ga": "Ga_d", "H0p75": "H.75"}

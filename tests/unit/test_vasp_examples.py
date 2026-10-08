"""Tier-1/2 tests for the VASP example campaigns."""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import pytest

import psteros

EXAMPLES = Path(__file__).resolve().parents[2] / "examples"


def load(relative: str, name: str):
    spec = importlib.util.spec_from_file_location(name, EXAMPLES / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class FakeStructure:
    def __init__(self, structure):
        self.structure = structure

    def get_pymatgen_structure(self):
        return self.structure


def test_sno2_phase_diagram_from_vasp_outputs() -> None:
    analysis = load("vasp_surface_phase_diagram/phase_diagram.py", "vasp_sno2_phase_diagram")
    bulk, metal = psteros.rutile_sno2_bulk(), psteros.alpha_sn_bulk()
    misc = lambda energy: {"total_energies": {"energy_extrapolated": energy}}  # noqa: E731
    refs = {
        "sno2_bulk_static_misc": misc(-2 * 20.0), "sno2_bulk_relaxed_structure": FakeStructure(bulk),
        "alpha_sn_static_misc": misc(-8 * 4.0), "alpha_sn_relaxed_structure": FakeStructure(metal),
        "o2_static_misc": misc(-9.8), "o2_relaxed_structure": FakeStructure(psteros.triplet_o2_cell()),
    }
    slabs = {}
    for termination, label in zip(("o", "sno", "sn2o"), analysis.SLABS):
        slab, _ = psteros.sno2_110_slab(termination=termination, triple_layers=3)
        n_sn = slab.composition["Sn"]
        n_o = slab.composition["O"]
        slabs[f"{label}_static_misc"] = misc(-20.0 * n_sn / 2 - 4.9 * (n_o - 2 * n_sn) + 3.0)
        slabs[f"{label}_relaxed_structure"] = slab

    class Graph:
        def __init__(self, outputs):
            self.outputs = outputs

    diagram = analysis.build_diagram(Graph(refs), Graph(slabs))
    assert set(diagram.curves) == set(analysis.SLABS)
    assert diagram.references.oxygen_poor_limit_ev < 0


def _arguments(code):
    return argparse.Namespace(code=code, potential_family="PBE", computer="local", queue="debug", ranks=1,
                              walltime=600, refs_pk=None, a=4.737, c=3.186, submit=False)


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_vasp_examples_build_their_graphs(tmp_path) -> None:
    pytest.importorskip("aiida_vasp")
    from aiida import orm

    computer = orm.Computer(
        label=f"local-{tmp_path.name}", hostname="localhost", transport_type="core.local",
        scheduler_type="core.direct", workdir=str(tmp_path),
    ).store()
    orm.InstalledCode(computer=computer, filepath_executable="/bin/true", label="vasp",
                      default_calc_job_plugin="vasp.vasp").store()
    code = f"vasp@{computer.label}"

    campaign = load("vasp_surface_phase_diagram/campaign.py", "vasp_sno2_campaign")
    args = _arguments(code)
    refs = campaign.build_refs(args)
    assert {"sno2_bulk_relax_vasp", "sno2_bulk_static_vasp", "o2_static_vasp"} <= {t.name for t in refs.tasks}
    o2 = refs.tasks["o2_static_vasp"].inputs.parameters.value.get_dict()["incar"]
    assert o2["ISPIN"] == 2 and o2["NSW"] == 0
    slabs = campaign.build_slabs(args)
    dynamics = slabs.tasks["slab_o_relax_vasp"].inputs.parameters.value.get_dict()["dynamics"]["positions_dof"]
    assert [False, False, False] in dynamics and [True, True, True] in dynamics

    zno = load("charge_neutral_terminations/zno_vasp_study.py", "zno_vasp_study")
    zno.CODE, zno.COMPUTER = code, computer.label
    study = zno.study()
    relax, static = zno.recipes(study)
    graph = psteros.build_relax_static_workgraph(study.structures, relax, static)
    assert graph.tasks["ref_Zn_relax_vasp"].inputs.parameters.value.get_dict()["incar"]["ISIF"] == 3
    assert "ISIF" not in graph.tasks["ref_Zn_static_vasp"].inputs.parameters.value.get_dict()["incar"]

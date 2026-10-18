"""Part D check of docs/plans/MACE_SURFACE_PLAN.md: one end-to-end run of ``psteros.build_mlip_surface_workgraph``.

A wiring check, not a result: two terminations (o_bridge, o_cus) and the three references, the CP2K recipe of the
labels on Phileas (at most 2 jobs at a time), a 2-pass replay fine-tune with one seed on Lovelace ``testegpu``, one
active-learning round and the DFT check. Checks after the run (``check``):

* the builder's first statics give the energies of the C2 graphs for the same structures (same recipe);
* every task finished; ``read_mlip_results`` and both sources of ``read_mlip_energies`` read the graph;
* the second fine-tune trained on the five first statics plus the five active-learning statics.

Usage (PYTHONPATH=<worktree>):
    python c8_builder.py submit [--dry]
    python c8_builder.py check
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import psteros

import c2_cp2k
import c4_gpu
import cluster

HERE = Path(__file__).resolve().parent
RECORD = HERE / "c8_builder.json"
SLABS = ("o_bridge", "o_cus")
QUEUE = c4_gpu.QUEUES["testegpu"]


def structures() -> dict:
    refs, labels = c4_gpu.graph("c2_refs.json"), c4_gpu.graph("c2_labels.json")
    result = {label: c4_gpu.input_structure(refs, label) for label in c4_gpu.REFERENCES}
    result.update({label: c4_gpu.input_structure(labels, label) for label in SLABS})
    return result


def recipes():
    static = c2_cp2k.config("ENERGY_FORCE", "ruo2_c8", cluster.MAX_CP2K_JOBS,
                            {"o2": psteros.CalculationOverride(parameters=c2_cp2k.O2_PARAMETERS)})
    computer = c4_gpu.COMPUTERS["lovelace"]
    fine_tune = psteros.MaceFineTuneConfig(
        foundation_model=computer["foundation"], foundation_model_sha256=c4_gpu.FOUNDATION_SHA256,
        replay_file=c4_gpu.REPLAY_FILE, epochs=2, seeds=(20260920,))
    job = psteros.MlipJobConfig(computer["code"], psteros.ExecutionPolicy(
        computer="lovelace", queue=QUEUE["queue_name"], with_mpi=False,
        max_wallclock_seconds=QUEUE["max_wallclock_seconds"], resources=QUEUE["resources"],
        custom_scheduler_commands=QUEUE["custom_scheduler_commands"],
        extra_options={"redirect_stderr": True, "import_sys_environment": False,
                       "environment_variables": {"OMP_NUM_THREADS": "8", "MKL_NUM_THREADS": "1",
                                                 "PYTHONUNBUFFERED": "1", "PYTHONNOUSERSITE": "1",
                                                 "PYTORCH_CUDA_ALLOC_CONF": "backend:cudaMallocAsync"}}))
    return static, fine_tune, job


def submit(dry: bool) -> None:
    items = structures()
    static, fine_tune, job = recipes()
    fixed = {label: psteros.central_sites(items[label]) for label in SLABS}
    graph = psteros.build_mlip_surface_workgraph(items, static, fine_tune, job=job, metallic_bulk="bulk",
                                                 fixed_sites=fixed, submit=not dry)
    print(len(graph.tasks), "tasks")
    if not dry:
        RECORD.write_text(json.dumps({"workgraph_pk": graph.pk, "labels": list(items)}, indent=2) + "\n")
        print("WorkGraph", graph.pk)


def check() -> None:
    from aiida import orm

    record = json.loads(RECORD.read_text())
    graph, labels = orm.load_node(record["workgraph_pk"]), record["labels"]
    failed = [link.link_label for link in graph.base.links.get_outgoing().all()
              if hasattr(link.node, "is_finished_ok") and not link.node.is_finished_ok]
    first = psteros.read_cp2k_results(graph, labels)[0]  # the <label>_static_parameters outputs
    refs, slabs = c4_gpu.graph("c2_refs.json"), c4_gpu.graph("c2_labels.json")
    before = {**psteros.read_cp2k_results(refs, list(c4_gpu.REFERENCES))[0], **psteros.read_cp2k_results(slabs, list(SLABS))[0]}
    relaxations = psteros.read_mlip_results(graph, labels)
    model = psteros.read_mlip_energies(graph, labels, energy_source="mlip")
    dft = psteros.read_mlip_energies(graph, labels, energy_source="mlip_dft")
    training = {link.link_label: link.node for link in graph.base.links.get_outgoing().all()}
    frames = {name: training[name].inputs.settings.get_dict()["frames"]
              for name in ("mlip_round0_training", "mlip_round1_training")}
    report = {
        "workgraph_pk": graph.pk, "state": graph.process_state.value, "exit_status": graph.exit_status,
        "not_finished_ok": failed,
        "first_static_minus_c2_ev": {label: first[label] - before[label] for label in labels},
        "relaxations": {label: {"energy_ev": r.energy_ev, "max_force": r.max_force_ev_per_angstrom, "steps": r.steps,
                                "converged": r.converged} for label, r in relaxations.items()},
        "mlip_ev": model, "mlip_dft_ev": dft, "frames": frames,
    }
    (HERE / "c8_check.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    from aiida import load_profile

    load_profile(c2_cp2k.PROFILE)
    if sys.argv[1] == "submit":
        cluster.with_retry(submit, "--dry" in sys.argv)
    else:
        cluster.with_retry(check)

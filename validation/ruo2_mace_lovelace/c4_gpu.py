"""Steps C4 and C5 of docs/plans/MACE_SURFACE_PLAN.md: MACE fine-tunes and relaxations on Lovelace GPUs.

Training frames come from the finished CP2K graphs of ``c2_cp2k.py``: the static calculations of bulk RuO2, hcp
Ru and O2 (on their CP2K-relaxed structures), of the six terminations as built and of bulk110; with ``--rattled``
also the 30 rattled terminations. Each structure is the input structure of its CP2K job (from the provenance);
energies in eV and forces in eV/A from the psteros readers. Bulk RuO2 must be metallic (``require_metallic``).

Studies (one aiida-shell job each, ``gpu_driver.py`` on the node):
    check   replay fine-tune, one seed, 2 passes; no relaxation (Lovelace ``testegpu``): set-up and speed
    main    M1 (ordinary), M2 (replay), M3 (replay + rattled), 3 seeds each, then every model and M0 (MH-1
            zero-shot) relax the six terminations, bulk RuO2 and hcp Ru (cell fixed) and O2 (Bohr ``gpu_a100``,
            one A100, the families one after the other; user, 2026-10-17)

Usage (PYTHONPATH=<worktree>):
    python c4_gpu.py prepare check|main      # writes local/jobs/<study>/ (git-ignored)
    python c4_gpu.py submit check|main [--dry]
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import psteros
from psteros import mlip

import cluster

HERE = Path(__file__).resolve().parent
JOBS = HERE / "local" / "jobs"
PROFILE = "psteros_mace_ruo2"
FOUNDATION_SHA256 = "a522eb7f59c7879963d41586528f4980baf33e086c94aa92e3eafdeccad3be47"
REPLAY_FILE = ("/home/thiagotdorini/projects/research/ag3po4-h2o-o2-surfaces/projects/archived/"
               "ag3po4_h2o_o2_mace_active_learning/goals/03_mace_finetuning/replay/replay-data-mh-1-omat-pbe.xyz")
REFERENCES = ("bulk", "ru", "o2")
# the MACE venv and the MH-1 copy on each GPU computer (same mace-torch 0.3.16, torch 2.5.1+cu118, ASE 3.29)
COMPUTERS = {
    "lovelace": {"code": "mace-python-cu118@lovelace",
                 "runtime": "/home/lovelace/proj/proj962/dorinitt/envs/ag3po4_mace_cu118/bin/python",
                 "foundation": "/work/dorinitt/goal03_inputs/models/mace-mh-1.model"},
    "bohr": {"code": "mace-python-cu118@bohr",
             "runtime": "/home/des01/zacarias/dorinitt/envs/ag3po4_g001_mace_cu118/bin/python",
             "foundation": "/home/des01/zacarias/dorinitt/models/mace-mh-1/mace-mh-1.model"},
}
QUEUES = {
    "gpu_a100": {"resources": {"num_machines": 1, "num_cores_per_machine": 56, "num_mpiprocs_per_machine": 1},
                 "custom_scheduler_commands": "#PBS -j oe\n#PBS -l place=excl", "queue_name": "gpu_a100", "max_wallclock_seconds": 48 * 3600},
    "testegpu": {"resources": {"num_machines": 1, "num_cores_per_machine": 8, "num_mpiprocs_per_machine": 8,
                               "ngpus_per_machine": 1},
                 "custom_scheduler_commands": "#PBS -j oe", "queue_name": "testegpu", "max_wallclock_seconds": 1800},
    "duasgpus": {"resources": {"num_machines": 1, "num_cores_per_machine": 32, "num_mpiprocs_per_machine": 32,
                               "ngpus_per_machine": 2},
                 "custom_scheduler_commands": "#PBS -j oe", "queue_name": "duasgpus", "max_wallclock_seconds": 24 * 3600},
}
STUDIES = {"check": {"queue": "testegpu", "gpus": 1, "computer": "lovelace"},
           "main": {"queue": "gpu_a100", "gpus": 1, "computer": "bohr"}}


def graph(name: str):
    from aiida import orm

    return orm.load_node(json.loads((HERE / name).read_text())["workgraph_pk"])


def input_structure(node, label: str):
    """The input structure of the CP2K job of ``label`` (its static stage when there is one)."""

    from psteros.backends.cp2k import _stage_output

    retrieved = _stage_output(node, label, "retrieved")
    if retrieved is None:
        raise LookupError(f"{label} has no finished CP2K calculation yet")
    return retrieved.creator.inputs.structure.get_pymatgen()


def labelled(node, label: str) -> mlip.LabelledStructure:
    structure = input_structure(node, label)
    energy = psteros.read_cp2k_results(node, [label])[0][label]
    forces = psteros.read_cp2k_forces(node, [label])[label]
    return mlip.LabelledStructure(label, structure, energy, forces)


def frames(rattled: bool, skip_missing: bool = False) -> list[mlip.LabelledStructure]:
    """Training frames; ``skip_missing`` (the check study only) leaves out labels that are not finished."""

    refs, labels = graph("c2_refs.json"), graph("c2_labels.json")
    dos = psteros.read_cp2k_dos(refs, ["bulk"])["bulk"]
    psteros.require_metallic("bulk", partially_occupied=dos.partially_occupied_bins())
    result = [labelled(refs, label) for label in REFERENCES]
    names = json.loads((HERE / "c2_labels.json").read_text())["labels"]
    for label in names:
        if "_r" in label and label.rsplit("_r", 1)[1].isdigit() and not rattled:
            continue
        try:
            result.append(labelled(labels, label))
        except LookupError:
            if not skip_missing:
                raise
            print("check: left out (not finished):", label)
    return result


def relaxation_structures() -> dict:
    refs, labels = graph("c2_refs.json"), graph("c2_labels.json")
    structures = {label: input_structure(refs, label) for label in REFERENCES}
    for label in psteros.RUTILE_110_TERMINATIONS:
        structures[label] = input_structure(labels, label)
    return structures


def prepare(study: str) -> None:
    folder = JOBS / study
    if folder.exists():
        raise FileExistsError(f"{folder} exists; prepare each study once")
    # read every AiiDA input first (retrying SQLite locks), then write the folders
    plain = cluster.with_retry(frames, False, skip_missing=(study == "check"))
    rattled = cluster.with_retry(frames, True) if study == "main" else None
    structures = cluster.with_retry(relaxation_structures) if study == "main" else None
    train, valid, record = mlip.select_replay_structures(REPLAY_FILE, elements=["Ru", "O"])
    FOUNDATION = COMPUTERS[STUDIES[study]["computer"]]["foundation"]
    base = dict(foundation_model=FOUNDATION, foundation_model_sha256=FOUNDATION_SHA256, replay_file=REPLAY_FILE)
    finetune = []
    if study == "check":
        config = mlip.MaceFineTuneConfig(**base, epochs=2, seeds=(20260920,))
        mlip.prepare_finetune_job(folder / "ft_check", plain, config, name="ruo2_check",
                                  replay_train=train, replay_valid=valid)
        finetune = ["ft_check"]
        relax = None
    else:
        families = {
            "m1": (mlip.MaceFineTuneConfig(**{**base, "replay_file": None}, replay=False), plain),
            "m2": (mlip.MaceFineTuneConfig(**base), plain),
            "m3": (mlip.MaceFineTuneConfig(**base), rattled),
        }
        models = {"m0": (FOUNDATION, "omat_pbe")}
        for family, (config, items) in families.items():
            replay = dict(replay_train=train, replay_valid=valid) if config.replay else {}
            mlip.prepare_finetune_job(folder / f"ft_{family}", items, config, name=f"ruo2_{family}", **replay)
            finetune.append(f"ft_{family}")
            for seed in config.seeds:
                models[f"{family}_seed{seed}"] = (f"../ft_{family}/models/ruo2_{family}_seed{seed}.model", "Default")
        fixed = {label: psteros.central_sites(structures[label]) for label in psteros.RUTILE_110_TERMINATIONS}
        mlip.prepare_relax_job(folder / "relax", structures, models, fixed_sites=fixed)
        relax = "relax"
    # M2 (the proposed method) first; with two GPUs, M2 then M1 on one and M3 on the other
    order = [name for name in ("ft_m2", "ft_m3", "ft_m1") if name in finetune] or finetune
    (folder / "driver.json").write_text(json.dumps(
        {"gpus": STUDIES[study]["gpus"], "finetune": order, "relax": relax}, indent=2) + "\n")
    shutil.copy(HERE / "gpu_driver.py", folder / "gpu_driver.py")
    (folder / "replay_selection.json").write_text(json.dumps(record, indent=2) + "\n")
    print(sorted(path.name for path in folder.iterdir()), record["train"], "replay train,", record["valid"], "valid")


def submit(study: str, dry: bool) -> None:
    from aiida import engine, orm
    from aiida_shell import ShellJob

    folder = JOBS / study
    plan = json.loads((folder / "driver.json").read_text())
    folders = list(plan["finetune"]) + ([plan["relax"]] if plan["relax"] else [])
    nodes = {"driver": orm.SinglefileData(file=str(folder / "gpu_driver.py")),
             "driver_json": orm.SinglefileData(file=str(folder / "driver.json"))}
    filenames = {"driver": "gpu_driver.py", "driver_json": "driver.json"}
    for name in folders:
        nodes[name] = orm.FolderData(tree=str(folder / name))
        filenames[name] = name
    outputs = ["report"]  # gpu_driver.py collects summaries, logs, evaluations and relaxations there
    computer = COMPUTERS[STUDIES[study]["computer"]]
    code = orm.load_code(computer["code"])
    assert str(code.filepath_executable) == computer["runtime"] and code.with_mpi is False
    options = dict(QUEUES[STUDIES[study]["queue"]], withmpi=False, redirect_stderr=True, import_sys_environment=False,
                   environment_variables={"OMP_NUM_THREADS": "8", "MKL_NUM_THREADS": "1", "PYTHONUNBUFFERED": "1",
                                          "PYTHONNOUSERSITE": "1",
                                          "PYTORCH_CUDA_ALLOC_CONF": "backend:cudaMallocAsync"})
    inputs = dict(code=code, nodes=nodes, filenames=filenames, arguments=["-u", "{driver}"], outputs=outputs,
                  metadata={"label": f"ruo2_{study}", "description": f"RuO2 MACE {study} (C4/C5)", "options": options})
    if dry:
        inputs["metadata"].update(dry_run=True, store_provenance=False)
        _, node = engine.run_get_node(ShellJob, **inputs)
        script = (Path(node.dry_run_info["folder"]) / node.dry_run_info["script_filename"]).read_text()
        print(script)
        return
    node = engine.submit(ShellJob, **inputs)
    (HERE / f"c4_{study}.json").write_text(json.dumps({"shelljob_pk": node.pk, "queue": STUDIES[study]["queue"],
                                                        "driver": plan}, indent=2) + "\n")
    print("ShellJob", node.pk)


if __name__ == "__main__":
    from aiida import load_profile

    load_profile(PROFILE)
    action, study = sys.argv[1], sys.argv[2]
    assert study in STUDIES
    if action == "prepare":
        prepare(study)
    else:
        submit(study, "--dry" in sys.argv)

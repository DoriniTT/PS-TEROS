"""Where the CP2K calculations of the validation run (one CP2K build for every CP2K energy of the study).

2026-10-17: moved from Lovelace (CP2K 2024.3 module, ``par128``: 3 of 42 C1 jobs done in 80 min) to Phileas
(CP2K 2026.2 module, ``standard`` partition, account m26257, idle nodes), at most 2 CP2K jobs at a time (user).
The MACE jobs stay on the Lovelace GPUs.
"""

from __future__ import annotations

import psteros

TAG = "phileas"
CP2K_CODE = "CP2K-2026.2@phileas"
MAX_CP2K_JOBS = 2
_SHAPES = {
    # 8 MPI ranks on a shared node, short checks
    "debug": {"resources": {"num_machines": 1, "num_mpiprocs_per_machine": 8, "num_cores_per_mpiproc": 1}},
    # one standard node: 96 MPI ranks, 230 GiB (Slurm keeps 16 GB of the 257 GB)
    "node": {"resources": {"num_machines": 1, "num_mpiprocs_per_machine": 96, "num_cores_per_mpiproc": 1},
             "max_memory_kb": 241172480},
}


def execution(shape: str, concurrent: int, walltime_seconds: int) -> psteros.ExecutionPolicy:
    if concurrent > MAX_CP2K_JOBS:
        raise ValueError(f"at most {MAX_CP2K_JOBS} CP2K jobs at a time")
    extra = {key: value for key, value in _SHAPES[shape].items() if key != "resources"}
    return psteros.ExecutionPolicy(
        computer="phileas", queue="standard", account="m26257", max_concurrent_jobs=concurrent,
        max_wallclock_seconds=walltime_seconds, resources=_SHAPES[shape]["resources"], extra_options=extra)

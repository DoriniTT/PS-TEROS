"""examples/mlip_routes/node_local_python.sh: unpack an environment once per node and run its Python."""

from __future__ import annotations

import os
import shutil
import subprocess
import tarfile
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "examples" / "mlip_routes" / "node_local_python.sh"


def _environment(root: Path, message: str) -> Path:
    python = root / "my_env" / "bin" / "python"
    python.parent.mkdir(parents=True)
    python.write_text(f'#!/bin/bash\necho "{message} $*"\n')
    python.chmod(0o755)
    tar = root / "my_env.tar"
    with tarfile.open(tar, "w") as archive:
        archive.add(root / "my_env", arcname="my_env")
    shutil.rmtree(root / "my_env")
    return tar


@pytest.mark.skipif(shutil.which("flock") is None, reason="needs flock")
def test_node_local_python_unpacks_once_and_runs(tmp_path) -> None:
    tar = _environment(tmp_path / "shared", "first")
    cache = tmp_path / "node"
    environment = {**os.environ, "PSTEROS_ENV_TAR": str(tar), "PSTEROS_ENV_CACHE": str(cache)}
    run = lambda: subprocess.run(["bash", str(SCRIPT), "-u", "finetune.py"], env=environment, capture_output=True,
                                 text=True, check=True).stdout.strip()
    assert run() == "first -u finetune.py"
    stamp = (cache / "my_env" / ".psteros_stamp").read_text()
    (cache / "my_env" / "bin" / "python").write_text('#!/bin/bash\necho "reused $*"\n')
    assert run() == "reused -u finetune.py"  # same tar: the unpacked copy is kept
    os.utime(tar, (1_000_000, 1_000_000))
    assert run() == "first -u finetune.py"  # the tar changed: unpacked again
    assert (cache / "my_env" / ".psteros_stamp").read_text() != stamp
    missing = subprocess.run(["bash", str(SCRIPT)], env={**os.environ, "PSTEROS_ENV_TAR": ""}, capture_output=True,
                             text=True)
    assert missing.returncode != 0 and "set PSTEROS_ENV_TAR" in missing.stderr

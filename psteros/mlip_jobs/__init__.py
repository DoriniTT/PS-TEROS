"""Self-contained scripts that psteros ships to the GPU computer (they need only ase, torch and mace-torch).

``finetune.py`` fine-tunes a MACE foundation model (optionally with a replay head); ``relax.py`` relaxes structures
with one or more models. :func:`script_path` gives the file of each, to copy into a job.
"""

from pathlib import Path

SCRIPTS = ("finetune.py", "relax.py")


def script_path(name: str) -> Path:
    """Path of the job script ``name`` (``"finetune.py"`` or ``"relax.py"``)."""

    if name not in SCRIPTS:
        raise ValueError(f"unknown job script {name!r}; choose from {SCRIPTS}")
    return Path(__file__).with_name(name)

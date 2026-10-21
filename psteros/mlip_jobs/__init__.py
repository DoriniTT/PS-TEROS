"""Self-contained scripts that psteros ships to the GPU computer (they need only ase, torch and mace-torch).

``finetune.py`` fine-tunes a MACE foundation model (optionally with a replay head), ``uma_finetune.py`` a UMA
model of fairchem and ``nequip_finetune.py`` a NequIP or Allegro package; ``relax.py`` relaxes structures with one or more models; ``vibrations.py`` gives the forces of
displaced structures for harmonic modes. The UMA and NequIP scripts need fairchem-core or nequip instead of mace-torch.
:func:`script_path` gives the file of each, to copy into a job.
"""

from pathlib import Path

SCRIPTS = ("finetune.py", "relax.py", "vibrations.py", "uma_finetune.py", "nequip_finetune.py")


def script_path(name: str) -> Path:
    """Path of the job script ``name`` (``"finetune.py"``, ``"relax.py"`` or ``"vibrations.py"``)."""

    if name not in SCRIPTS:
        raise ValueError(f"unknown job script {name!r}; choose from {SCRIPTS}")
    return Path(__file__).with_name(name)

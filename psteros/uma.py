"""UMA (fairchem) foundation models on the MLIP route (no AiiDA, no torch).

The MLIP route fine-tunes a foundation model on the DFT statics of a study. With
a :class:`UmaFineTuneConfig` instead of a :class:`psteros.MaceFineTuneConfig`,
that model is a UMA model of Meta's fairchem (``fairchem-core`` on the GPU
computer); everything else of the route is the same.

The fine-tune starts from the pretrained model, as the MACE one does. fairchem's
own fine-tune starts from the raw weights, builds a new energy head from random
weights and fits new element references; the job here starts from the averaged
weights that fairchem's inference uses, merges UMA's per-task head on ``task``
and copies it into the new head, and sets its element references to the
foundation's for ``task`` plus one shift per element, fitted in least squares to
E_DFT − E_UMA over the training frames, so that the untrained model reproduces
the foundation model (shifted), which the job checks before training. The replay structures are labelled by the
foundation model with the same shifts, so the DFT labels and the replay labels
share one energy reference.

UMA's weights are gated: accept the licence at https://huggingface.co/facebook/UMA
and download the checkpoint (``checkpoints/uma-s-1p2p1.pt``, ...) to the GPU
computer.
"""

from __future__ import annotations

import json
from collections.abc import Iterable, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any

from psteros.mace_models import MaceFoundationHead

UMA_REPOSITORY = "https://huggingface.co/facebook/UMA"


@dataclass(frozen=True)
class UmaFoundationModel:
    """A released UMA checkpoint: ``filename`` in the Hugging Face repository ``facebook/UMA``, its licence, its
    ``tasks`` (described like the heads of :class:`psteros.MaceFoundationModel`) and the default task."""

    name: str
    filename: str
    licence: str
    tasks: tuple[MaceFoundationHead, ...]
    default_task: str
    use: str

    def __post_init__(self) -> None:
        names = [task.name for task in self.tasks]
        if not names or len(set(names)) != len(names) or self.default_task not in names:
            raise ValueError(f"{self.name}: tasks must be distinct and include default_task")

    @property
    def url(self) -> str:
        return f"{UMA_REPOSITORY}/resolve/main/{self.filename}"

    @property
    def task_names(self) -> tuple[str, ...]:
        return tuple(task.name for task in self.tasks)

    def task(self, name: str | None = None) -> MaceFoundationHead:
        """The task ``name`` (the default task when ``None``); ``ValueError`` for a task the model does not have."""

        wanted = self.default_task if name is None else name
        for task in self.tasks:
            if task.name == wanted:
                return task
        raise ValueError(f"{self.name} has no task {wanted!r}; its tasks are {list(self.task_names)}")


_TASKS = (
    MaceFoundationHead("omat", "PBE", "mp", "OMat24 (inorganic materials)"),
    MaceFoundationHead("oc20", "RPBE", "none", "OC20 (catalysis: slabs and adsorbates)"),
    MaceFoundationHead("omol", "wB97M-V", "none", "OMol25 (molecules)", for_slabs=False),
)
_OC22 = MaceFoundationHead("oc22", "PBE", None, "OC22 (oxide catalysis: slabs and adsorbates)")  # UMA 1.2 and later

_MODELS = (
    UmaFoundationModel("uma-s-1p2p1", "checkpoints/uma-s-1p2p1.pt", "FAIR Chemistry License v1",
                       (*_TASKS[:2], _OC22, _TASKS[2]), "omat",
                       "The small UMA model (2.3 GB file)."),
    UmaFoundationModel("uma-m-1p1", "checkpoints/uma-m-1p1.pt", "FAIR Chemistry License v1", _TASKS, "omat",
                       "The medium UMA model (11 GB file), slower and heavier to fine-tune."),
)

#: UMA checkpoints of Meta's release (``facebook/UMA``), by name; the job script was written for fairchem-core 2.23.
#: The checkpoints also have the tasks ``oc25``, ``odac`` and ``omc``, left out here (their DFT is not described); a
#: :class:`UmaFineTuneConfig` without ``foundation_name`` can still name them.
UMA_FOUNDATION_MODELS = MappingProxyType({model.name: model for model in _MODELS})


def uma_foundation_model(name: str) -> UmaFoundationModel:
    """The model ``name`` of :data:`UMA_FOUNDATION_MODELS`; ``ValueError`` listing the known names otherwise."""

    try:
        return UMA_FOUNDATION_MODELS[name]
    except (KeyError, TypeError):
        raise ValueError(f"unknown UMA foundation model {name!r}; known: {sorted(UMA_FOUNDATION_MODELS)}") from None


def uma_level_warnings(model: str, task: str | None, calculation: Any, elements: Iterable[str]) -> list[str]:
    """:func:`psteros.foundation_level_warnings` for the task ``task`` of the UMA model ``model``."""

    from psteros.mace_models import _level_warnings

    entry = uma_foundation_model(model)
    info = entry.task(task)
    alternatives = [(f"{other.name} (task {item.name!r})", item) for other in _MODELS for item in other.tasks]
    return _level_warnings(f"{entry.name} (task {info.name!r})", info, calculation, elements, alternatives)


@dataclass(frozen=True)
class UmaFineTuneConfig:
    """Recipe of a UMA fine-tune on the MLIP route (the counterpart of :class:`psteros.MaceFineTuneConfig`).

    ``foundation_model`` is the path of the UMA checkpoint on the GPU computer (checked against
    ``foundation_model_sha256`` when given) and ``task`` the UMA task whose embedding and element references the
    fine-tune starts from (``"omat"`` for inorganic materials). With ``replay``, :func:`psteros.select_replay_structures`
    picks replay geometries from the local ``replay_file`` (any extended-XYZ file, e.g. the replay data of MACE-MH-1);
    the job labels them with the foundation model and trains on them together with the DFT frames, every
    ``replay_validation_every``-th one kept for validation.

    The training runs ``epochs`` passes with AdamW (``learning_rate``, ``weight_decay``), a cosine schedule,
    ``batch_size`` structures per step and fairchem's losses (per-atom energy MAE weighted ``energy_weight``,
    force L2 weighted ``forces_weight``); the kept model is the exponential moving average (``ema_decay``) of the
    last step, which fairchem's inference loads. The graph settings of the backbone (neighbours, cutoff) are the
    foundation's. One model per seed.
    ``foundation_name`` (optional) names a model of :data:`UMA_FOUNDATION_MODELS`.
    """

    foundation_model: str
    foundation_model_sha256: str | None = None
    task: str = "omat"
    replay: bool = True
    replay_file: str | None = None
    replay_min_atoms: int = 2
    replay_max_atoms: int = 128
    replay_validation_every: int = 5
    epochs: int = 20
    learning_rate: float = 1.0e-4
    weight_decay: float = 1.0e-3
    batch_size: int = 2
    energy_weight: float = 20.0
    forces_weight: float = 2.0
    ema_decay: float = 0.999
    seeds: tuple[int, ...] = (20260920,)
    foundation_name: str | None = None

    @classmethod
    def for_foundation(cls, name: str, foundation_model: str, *, task: str | None = None,
                       **options: Any) -> "UmaFineTuneConfig":
        """The recipe for the UMA model ``name`` of :data:`UMA_FOUNDATION_MODELS` (its default task unless given)."""

        model = uma_foundation_model(name)
        return cls(foundation_model, task=model.task(task).name, foundation_name=name, **options)

    def __post_init__(self) -> None:
        if not self.foundation_model or not self.task:
            raise ValueError("UmaFineTuneConfig needs foundation_model and task")
        if self.foundation_name is not None:
            uma_foundation_model(self.foundation_name).task(self.task)
        if self.replay and not self.replay_file:
            raise ValueError("UmaFineTuneConfig: replay=True needs replay_file")
        if not 1 <= self.replay_min_atoms <= self.replay_max_atoms:
            raise ValueError("UmaFineTuneConfig: need 1 <= replay_min_atoms <= replay_max_atoms")
        if self.replay_validation_every < 2:
            raise ValueError("UmaFineTuneConfig: replay_validation_every must be at least 2")
        for name in ("epochs", "batch_size"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"UmaFineTuneConfig: {name} must be a positive integer")
        for name in ("learning_rate", "energy_weight", "forces_weight"):
            if getattr(self, name) <= 0:
                raise ValueError(f"UmaFineTuneConfig: {name} must be positive")
        if self.weight_decay < 0:
            raise ValueError("UmaFineTuneConfig: weight_decay must not be negative")
        if not 0 < self.ema_decay < 1:
            raise ValueError("UmaFineTuneConfig: ema_decay must be between 0 and 1")
        seeds = tuple(self.seeds)
        if not seeds or len(set(seeds)) != len(seeds) or any(isinstance(s, bool) or not isinstance(s, int) for s in seeds):
            raise ValueError("UmaFineTuneConfig: seeds must be distinct integers, at least one")
        object.__setattr__(self, "seeds", seeds)


def prepare_uma_finetune_job(
    directory: str | Path,
    frames: Sequence[Any],
    config: UmaFineTuneConfig,
    *,
    name: str,
    replay_train: Sequence[Any] = (),
    replay_valid: Sequence[Any] = (),
    device: str = "cuda",
    seeds: Sequence[int] | None = None,
) -> Path:
    """Write the folder of one UMA fine-tune job (``job.json``, ``train.extxyz``, replay geometries, ``finetune.py``).

    As :func:`psteros.prepare_finetune_job`: one run per seed (``<name>_seed<seed>``), the ``"train"`` frames of
    ``frames`` (:class:`psteros.LabelledStructure`) and, with ``config.replay``, the replay geometries. The job
    writes each model to ``models/<name>_seed<seed>.model`` (a fairchem inference checkpoint). Returns the folder,
    which must not exist yet.
    """

    from psteros.mlip import write_extxyz
    from psteros.mlip_jobs import script_path

    train = [frame for frame in frames if frame.role == "train"]
    if not train:
        raise ValueError("no training frames (role 'train')")
    if config.replay and (not replay_train or not replay_valid):
        raise ValueError("a replay fine-tune needs replay_train and replay_valid structures")
    if not config.replay and (replay_train or replay_valid):
        raise ValueError("replay structures were given to a fine-tune without replay")
    if device not in ("cuda", "cpu"):
        raise ValueError("device must be 'cuda' or 'cpu'")
    chosen = tuple(config.seeds if seeds is None else seeds)
    unknown = sorted(set(chosen).difference(config.seeds))
    if unknown:
        raise ValueError(f"seeds {unknown} are not seeds of the recipe {config.seeds}")
    folder = Path(directory)
    folder.mkdir(parents=True, exist_ok=False)
    write_extxyz(train, folder / "train.extxyz")
    if config.replay:
        from ase.io import write

        write(folder / "replay_train_geometries.extxyz", list(replay_train), format="extxyz")
        write(folder / "replay_valid_geometries.extxyz", list(replay_valid), format="extxyz")
    recipe = asdict(config)
    job = {
        "family": "uma", "foundation_model": config.foundation_model,
        "foundation_model_sha256": config.foundation_model_sha256, "task": config.task, "replay": config.replay,
        "device": device,
        "training": {key: recipe[key] for key in ("epochs", "learning_rate", "weight_decay", "batch_size",
                                                  "energy_weight", "forces_weight", "ema_decay")},
        "runs": [{"name": f"{name}_seed{seed}", "seed": seed} for seed in chosen],
    }
    (folder / "job.json").write_text(json.dumps(job, indent=2) + "\n")
    (folder / "finetune.py").write_text(script_path("uma_finetune.py").read_text())
    return folder

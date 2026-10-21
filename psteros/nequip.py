"""NequIP and Allegro foundation models (the NequIP framework) on the MLIP route (no AiiDA, no torch).

The MLIP route fine-tunes a foundation model on the DFT statics of a study. With
a :class:`NequipFineTuneConfig` instead of a :class:`psteros.MaceFineTuneConfig`,
that model is a packaged NequIP or Allegro model (``nequip`` and, for Allegro,
``nequip-allegro`` on the GPU computer); everything else of the route is the same.

The fine-tune follows the NequIP fine-tuning guide (the packaged model, its
cutoff and type names, a lower learning rate) and starts from the foundation
model, as the MACE one does: the per-type energy shifts of the model are the
foundation's plus one shift per element, fitted in least squares to
E_DFT − E_foundation over the training frames, so that the untrained model
reproduces the foundation model (shifted) exactly, which the job checks before
training. The replay structures are labelled by the foundation model with the
same shifts.

The packages are on https://www.nequip.net (CC BY 4.0); download one with
``nequip-package info nequip.net:mir-group/NequIP-OAM-L:0.1`` on a computer with
internet (it lands in ``~/.nequip/model_cache``) and give its path.
"""

from __future__ import annotations

import json
from collections.abc import Iterable, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any

from psteros.mace_models import MaceFoundationHead

NEQUIP_SITE = "https://www.nequip.net/models"


@dataclass(frozen=True)
class NequipFoundationModel:
    """A NequIP-framework foundation model of nequip.net: its ``model_id`` (``mir-group/NequIP-OAM-L:0.1``), the
    architecture (``"NequIP"`` or ``"Allegro"``, which needs ``nequip-allegro``), licence and the DFT of its training
    data (described like a head of :class:`psteros.MaceFoundationModel`)."""

    name: str
    model_id: str
    architecture: str
    licence: str
    level: MaceFoundationHead
    use: str

    def __post_init__(self) -> None:
        if self.architecture not in ("NequIP", "Allegro"):
            raise ValueError(f"{self.name}: architecture must be 'NequIP' or 'Allegro'")

    @property
    def url(self) -> str:
        return f"{NEQUIP_SITE}/{self.model_id}"

    @property
    def download(self) -> str:
        """The ``nequip.net:`` identifier that nequip's tools download."""

        return f"nequip.net:{self.model_id}"


_OAM = MaceFoundationHead("default", "PBE", "mp", "OMat24, then sAlex and MPtrj")

_MODELS = tuple(
    NequipFoundationModel(f"{architecture.lower()}-oam-{size.lower()}", f"mir-group/{architecture}-OAM-{size}:0.1",
                          architecture, "CC BY 4.0", _OAM, use)
    for architecture, size, use in (
        ("NequIP", "S", "The smallest and fastest NequIP OAM model."),
        ("NequIP", "M", "A medium NequIP OAM model."),
        ("NequIP", "L", "The NequIP OAM model recommended by its authors for most inorganic solids."),
        ("NequIP", "XL", "The largest NequIP OAM model."),
        ("Allegro", "L", "Allegro (strictly local, the fastest for large slabs); needs nequip-allegro."),
    )
)

#: The NequIP-framework foundation models of nequip.net, by name; the job script was written for nequip 0.19.
NEQUIP_FOUNDATION_MODELS = MappingProxyType({model.name: model for model in _MODELS})


def nequip_foundation_model(name: str) -> NequipFoundationModel:
    """The model ``name`` of :data:`NEQUIP_FOUNDATION_MODELS`; ``ValueError`` listing the known names otherwise."""

    try:
        return NEQUIP_FOUNDATION_MODELS[name]
    except (KeyError, TypeError):
        raise ValueError(f"unknown NequIP foundation model {name!r}; known: {sorted(NEQUIP_FOUNDATION_MODELS)}") from None


def nequip_level_warnings(model: str, calculation: Any, elements: Iterable[str]) -> list[str]:
    """:func:`psteros.foundation_level_warnings` for a model of :data:`NEQUIP_FOUNDATION_MODELS`."""

    from psteros.mace_models import _level_warnings

    entry = nequip_foundation_model(model)
    alternatives = [(other.name, other.level) for other in _MODELS]
    return _level_warnings(entry.name, entry.level, calculation, elements, alternatives)


@dataclass(frozen=True)
class NequipFineTuneConfig:
    """Recipe of a NequIP or Allegro fine-tune on the MLIP route (the counterpart of :class:`psteros.MaceFineTuneConfig`).

    ``foundation_model`` is the path of the packaged model (``.nequip.zip``) on the GPU computer, checked against
    ``foundation_model_sha256`` when given. With ``replay``, :func:`psteros.select_replay_structures` picks replay
    geometries from the local ``replay_file`` (any extended-XYZ file, e.g. the replay data of MACE-MH-1); the job
    labels them with the foundation model and trains on them together with the DFT frames, every
    ``replay_validation_every``-th one kept for validation.

    The training runs ``epochs`` passes of ``nequip-train`` with Adam at the constant ``learning_rate``,
    ``batch_size`` structures per step and nequip's energy-force loss (per-atom energy weighted ``energy_weight``,
    forces ``forces_weight``); the kept model is the exponential moving average of the weights (``ema_decay``) after
    the last pass, packaged as ``.nequip.zip``. The cutoff and the atom types are the foundation's. One model per
    seed. ``foundation_name`` (optional) names a model of :data:`NEQUIP_FOUNDATION_MODELS`.
    """

    foundation_model: str
    foundation_model_sha256: str | None = None
    replay: bool = True
    replay_file: str | None = None
    replay_min_atoms: int = 2
    replay_max_atoms: int = 128
    replay_validation_every: int = 5
    epochs: int = 30
    learning_rate: float = 1.0e-4
    batch_size: int = 2
    energy_weight: float = 1.0
    forces_weight: float = 1.0
    ema_decay: float = 0.999
    seeds: tuple[int, ...] = (20260920,)
    foundation_name: str | None = None

    @classmethod
    def for_foundation(cls, name: str, foundation_model: str, **options: Any) -> "NequipFineTuneConfig":
        """The recipe for the model ``name`` of :data:`NEQUIP_FOUNDATION_MODELS`, packaged at ``foundation_model``."""

        nequip_foundation_model(name)
        return cls(foundation_model, foundation_name=name, **options)

    def __post_init__(self) -> None:
        if not self.foundation_model:
            raise ValueError("NequipFineTuneConfig needs foundation_model")
        if not str(self.foundation_model).endswith(".nequip.zip"):
            raise ValueError("NequipFineTuneConfig: foundation_model must be a packaged model (.nequip.zip)")
        if self.foundation_name is not None:
            nequip_foundation_model(self.foundation_name)
        if self.replay and not self.replay_file:
            raise ValueError("NequipFineTuneConfig: replay=True needs replay_file")
        if not 1 <= self.replay_min_atoms <= self.replay_max_atoms:
            raise ValueError("NequipFineTuneConfig: need 1 <= replay_min_atoms <= replay_max_atoms")
        if self.replay_validation_every < 2:
            raise ValueError("NequipFineTuneConfig: replay_validation_every must be at least 2")
        for name in ("epochs", "batch_size"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"NequipFineTuneConfig: {name} must be a positive integer")
        for name in ("learning_rate", "energy_weight", "forces_weight"):
            if getattr(self, name) <= 0:
                raise ValueError(f"NequipFineTuneConfig: {name} must be positive")
        if not 0 < self.ema_decay < 1:
            raise ValueError("NequipFineTuneConfig: ema_decay must be between 0 and 1")
        seeds = tuple(self.seeds)
        if not seeds or len(set(seeds)) != len(seeds) or any(isinstance(s, bool) or not isinstance(s, int) for s in seeds):
            raise ValueError("NequipFineTuneConfig: seeds must be distinct integers, at least one")
        object.__setattr__(self, "seeds", seeds)


def prepare_nequip_finetune_job(
    directory: str | Path,
    frames: Sequence[Any],
    config: NequipFineTuneConfig,
    *,
    name: str,
    replay_train: Sequence[Any] = (),
    replay_valid: Sequence[Any] = (),
    device: str = "cuda",
    seeds: Sequence[int] | None = None,
) -> Path:
    """Write the folder of one NequIP/Allegro fine-tune job (``job.json``, ``train.extxyz``, replay geometries,
    ``finetune.py``), as :func:`psteros.prepare_finetune_job`. The job writes each model to
    ``models/<name>_seed<seed>.nequip.zip`` (a nequip package). Returns the folder, which must not exist yet.
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
        "family": "nequip", "foundation_model": config.foundation_model,
        "foundation_model_sha256": config.foundation_model_sha256, "replay": config.replay, "device": device,
        "training": {key: recipe[key] for key in ("epochs", "learning_rate", "batch_size", "energy_weight",
                                                  "forces_weight", "ema_decay")},
        "runs": [{"name": f"{name}_seed{seed}", "seed": seed} for seed in chosen],
    }
    (folder / "job.json").write_text(json.dumps(job, indent=2) + "\n")
    (folder / "finetune.py").write_text(script_path("nequip_finetune.py").read_text())
    return folder

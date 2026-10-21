"""The MACE foundation models of the MLIP route (no AiiDA, no torch).

The MLIP route fine-tunes a MACE foundation model on the DFT statics of a
study (:class:`psteros.MaceFineTuneConfig`). Any MACE model file works; this
module lists the released models that have been run through the route's
fine-tune job, with the DFT each of their heads was trained on, so that a
recipe can name one (``MaceFineTuneConfig.for_foundation("mace-omat-0-medium",
path)``) and be checked against the DFT of the study.

Why the level of theory matters: the job fits one energy shift per element
(the E0s of the new head) to DFT − foundation, which absorbs any constant
offset between the two, such as different pseudopotentials or a different
reference energy. It does not absorb a different functional or a different
+U scheme, which change the shape of the energy surface; the fine-tune then
has to learn that difference from the few DFT labels of the study. A head
trained on the functional of the study starts closer.

The Materials Project +U scheme, used by MPtrj, sAlex and OMat24 (and so by
MACE-MP-0, MPA-0, OMAT-0 and the ``omat_pbe`` head of MACE-MH-1), puts +U on
the d states of Co, Cr, Fe, Mn, Mo, Ni, V and W in oxides and fluorides only;
other compounds are plain PBE there.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from types import MappingProxyType
from typing import Any

#: Metals that carry +U in oxides and fluorides in the Materials Project scheme.
MP_HUBBARD_METALS = ("Co", "Cr", "Fe", "Mn", "Mo", "Ni", "V", "W")
#: Anions that switch the Materials Project +U on.
MP_HUBBARD_ANIONS = ("O", "F")
HUBBARD_SCHEMES = ("mp", "none")

_RELEASES = "https://github.com/ACEsuit/mace-foundations/releases/download"


@dataclass(frozen=True)
class MaceFoundationHead:
    """One head of a MACE foundation model and the DFT of its training data.

    ``name`` is the head as the model file names it (``"omat_pbe"``; most
    single-head models call theirs ``"default"``, and MACE calls the head of
    a file without heads, such as MACE-MP-0 of 2023, ``"Default"``). ``functional`` is the exchange-correlation functional
    of the training data (``"PBE"``, ``"r2SCAN"``, ``"wB97M-VV10"``, ...).
    ``hubbard`` is how +U was used: ``"mp"`` (the Materials Project scheme,
    see :data:`MP_HUBBARD_METALS`), ``"none"``, or ``None`` when the release
    does not say. ``for_slabs`` is false for heads trained on molecules only.
    """

    name: str
    functional: str
    hubbard: str | None
    training_data: str
    for_slabs: bool = True

    def __post_init__(self) -> None:
        if not self.name or not self.functional or not self.training_data:
            raise ValueError("MaceFoundationHead needs name, functional and training_data")
        if self.hubbard is not None and self.hubbard not in HUBBARD_SCHEMES:
            raise ValueError(f"MaceFoundationHead {self.name}: hubbard must be one of {HUBBARD_SCHEMES} or None")


@dataclass(frozen=True)
class MaceFoundationModel:
    """A released MACE foundation model.

    ``url`` is the model file of the release and ``sha256`` its checksum
    (the fine-tune job checks the file against it). ``licence`` is ``"MIT"``
    or ``"ASL"`` (Academic Software License). ``min_mace_torch`` is the
    oldest mace-torch that loads the file. ``default_head`` is the head the
    route uses unless told otherwise. ``replay_url`` is the replay data of
    the release (extended XYZ); the fine-tune job relabels the geometries it
    picks with the foundation model, so any of these files works with any
    model, but the model's own file holds the structures it was trained on.
    """

    name: str
    url: str
    sha256: str
    licence: str
    min_mace_torch: str
    heads: tuple[MaceFoundationHead, ...]
    default_head: str
    replay_url: str
    use: str

    def __post_init__(self) -> None:
        names = [head.name for head in self.heads]
        if not names or len(set(names)) != len(names):
            raise ValueError(f"{self.name}: heads must be distinct and at least one")
        if self.default_head not in names:
            raise ValueError(f"{self.name}: default_head {self.default_head!r} is not one of {names}")
        if len(self.sha256) != 64:
            raise ValueError(f"{self.name}: sha256 must be 64 hexadecimal characters")

    @property
    def filename(self) -> str:
        """File name of the model in its release."""

        return self.url.rsplit("/", 1)[-1]

    @property
    def head_names(self) -> tuple[str, ...]:
        return tuple(head.name for head in self.heads)

    def head(self, name: str | None = None) -> MaceFoundationHead:
        """The head ``name`` (the default head when ``None``); ``ValueError`` for a head the model does not have."""

        wanted = self.default_head if name is None else name
        for head in self.heads:
            if head.name == wanted:
                return head
        raise ValueError(f"{self.name} has no head {wanted!r}; its heads are {list(self.head_names)}")


def _single(name: str, functional: str, hubbard: str | None, training_data: str) -> tuple[MaceFoundationHead, ...]:
    return (MaceFoundationHead(name, functional, hubbard, training_data),)


_MP_REPLAY = f"{_RELEASES}/mace_mp_0b/mp_traj_combined.xyz"

_MODELS = (
    MaceFoundationModel(
        name="mace-mh-1",
        url=f"{_RELEASES}/mace_mh_1/mace-mh-1.model",
        sha256="a522eb7f59c7879963d41586528f4980baf33e086c94aa92e3eafdeccad3be47",
        licence="ASL",
        min_mace_torch="0.3.14",
        heads=(
            MaceFoundationHead("omat_pbe", "PBE", "mp", "OMat24 (replay)"),
            MaceFoundationHead("mp_pbe_refit_add", "PBE", "mp", "MPtrj"),
            MaceFoundationHead("oc20_usemppbe", "PBE", None, "OC20: metal slabs and adsorbates"),
            MaceFoundationHead("matpes_r2scan", "r2SCAN", "none", "MatPES r2SCAN"),
            MaceFoundationHead("omol", "wB97M-VV10", "none", "1 % of OMol25 (molecules)", for_slabs=False),
            MaceFoundationHead("spice_wB97M", "wB97M-D3(BJ)", "none", "SPICE (molecules)", for_slabs=False),
        ),
        default_head="omat_pbe",
        replay_url=f"{_RELEASES}/mace_mh_1/replay-data-mh-1-omat-pbe.xyz",
        use="The default of the route: one model for bulk, surfaces and molecules; heads for r2SCAN and adsorbates.",
    ),
    MaceFoundationModel(
        name="mace-omat-0-medium",
        url=f"{_RELEASES}/mace_omat_0/mace-omat-0-medium.model",
        sha256="d4b14be9afa294eebdbe31a0280b26a0fa29715771e978cbfd3ca24e0d90307a",
        licence="ASL",
        min_mace_torch="0.3.10",
        heads=_single("default", "PBE", "mp", "OMat24"),
        default_head="default",
        replay_url=f"{_RELEASES}/mace_omat_0/mp_traj_combined_omat.xyz",
        use="Single-head model trained on far-from-equilibrium OMat24 structures (PBE/PBE+U).",
    ),
    MaceFoundationModel(
        name="mace-omat-0-small",
        url=f"{_RELEASES}/mace_omat_0/mace-omat-0-small.model",
        sha256="0abfde07862cf1e93b8b4d03cb702f29ce9c344ff2fc4de2ec0d7166d6c113a5",
        licence="ASL",
        min_mace_torch="0.3.10",
        heads=_single("omat_pbe", "PBE", "mp", "OMat24"),
        default_head="omat_pbe",
        replay_url=f"{_RELEASES}/mace_omat_0/mp_traj_combined_omat.xyz",
        use="Smaller, faster OMAT-0 for large slabs.",
    ),
    MaceFoundationModel(
        name="mace-matpes-pbe-0",
        url=f"{_RELEASES}/mace_matpes_0/MACE-matpes-pbe-omat-ft.model",
        sha256="e618ad582b84239905b9c3b77ce6e9ce111b0ecd1533223a1a6aac7a696b8aa0",
        licence="ASL",
        min_mace_torch="0.3.10",
        heads=_single("default", "PBE", "none", "MatPES PBE (OMAT-0 fine-tuned)"),
        default_head="default",
        replay_url=f"{_RELEASES}/mace_matpes_0/matpes-pbe-replay-data.xyz",
        use="Plain PBE without +U: for a DFT recipe without +U on transition-metal oxides.",
    ),
    MaceFoundationModel(
        name="mace-matpes-r2scan-0",
        url=f"{_RELEASES}/mace_matpes_0/MACE-matpes-r2scan-omat-ft.model",
        sha256="8f147ecffa1d06a696e5648b13b81abda9c96fda9d0fb3faba3cd3fb27d9bba9",
        licence="ASL",
        min_mace_torch="0.3.10",
        heads=_single("default", "r2SCAN", "none", "MatPES r2SCAN (OMAT-0 fine-tuned)"),
        default_head="default",
        replay_url=f"{_RELEASES}/mace_matpes_0/matpes-r2scan-replay-data.extxyz",
        use="For r2SCAN DFT statics.",
    ),
    MaceFoundationModel(
        name="mace-mpa-0-medium",
        url=f"{_RELEASES}/mace_mpa_0/mace-mpa-0-medium.model",
        sha256="75428afe3a1d7d8062e19bcaabd5c433623cabf308242ec9fb493e38604fb638",
        licence="MIT",
        min_mace_torch="0.3.10",
        heads=_single("default", "PBE", "mp", "MPtrj and sAlex"),
        default_head="default",
        replay_url=_MP_REPLAY,
        use="The best MIT-licensed model (PBE/PBE+U).",
    ),
    MaceFoundationModel(
        name="mace-mp-0b3-medium",
        url=f"{_RELEASES}/mace_mp_0b3/mace-mp-0b3-medium.model",
        sha256="2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012",
        licence="MIT",
        min_mace_torch="0.3.9",
        heads=_single("default", "PBE", "mp", "MPtrj"),
        default_head="default",
        replay_url=_MP_REPLAY,
        use="MIT-licensed, recommended by its authors for fine-tuning; weak on surfaces without it.",
    ),
    MaceFoundationModel(
        name="mace-mp-0-medium",
        url=f"{_RELEASES}/mace_mp_0/2023-12-03-mace-128-L1_epoch-199.model",
        sha256="01bfe22100139f424713cf921144e5509cbe353d67aa9fa1be9c6e1e0ed35845",
        licence="MIT",
        min_mace_torch="0.3.6",
        heads=_single("Default", "PBE", "mp", "MPtrj"),
        default_head="Default",
        replay_url=_MP_REPLAY,
        use="The original MACE-MP-0 (2023), to compare with earlier work.",
    ),
)

#: The MACE foundation models the route has been run with (fine-tune and relaxation jobs, mace-torch 0.3.16), by name.
MACE_FOUNDATION_MODELS: Mapping[str, MaceFoundationModel] = MappingProxyType({model.name: model for model in _MODELS})


def mace_foundation_model(name: str) -> MaceFoundationModel:
    """The model ``name`` of :data:`MACE_FOUNDATION_MODELS`; ``ValueError`` listing the known names otherwise."""

    try:
        return MACE_FOUNDATION_MODELS[name]
    except (KeyError, TypeError):
        raise ValueError(f"unknown MACE foundation model {name!r}; known: {sorted(MACE_FOUNDATION_MODELS)}") from None


def _true(value: Any) -> bool:
    if isinstance(value, str):
        return value.strip().strip(".").upper() in ("TRUE", "T")
    return bool(value)


_FUNCTIONAL_NAMES = {"PBE": "PBE", "PE": "PBE", "PBESOL": "PBEsol", "PS": "PBEsol", "RPBE": "RPBE", "RP": "RPBE",
                     "R2SCAN": "r2SCAN", "SCAN": "SCAN"}


def _has_key(tree: Any, fragment: str) -> bool:
    if isinstance(tree, Mapping):
        return any(fragment in str(key).upper() or _has_key(value, fragment) for key, value in tree.items())
    if isinstance(tree, (list, tuple)):
        return any(_has_key(item, fragment) for item in tree)
    return False


def recipe_level_of_theory(calculation: Any) -> tuple[str | None, bool | None]:
    """``(functional, uses_hubbard_u)`` of a QE, VASP or CP2K recipe; ``None`` where the recipe does not tell.

    VASP: ``LHFCALC`` (hybrid), ``METAGGA`` (``R2SCAN``, ``SCAN``), ``GGA``
    (``PE``, ``PS``, ``RP``), else the functional of the POTCAR family;
    ``LDAU`` for +U. QE: ``SYSTEM/input_dft``, else the pseudopotential
    family's name; ``SYSTEM/lda_plus_u``. CP2K: ``FORCE_EVAL/DFT/XC/XC_FUNCTIONAL``
    (with the ``PARAMETRIZATION`` of ``PBE``, or ``R2SCAN`` in a LIBXC
    functional); any ``*PLUS_U*`` section or keyword for +U.
    """

    from psteros.config import Cp2kCalculationConfig, QeCalculationConfig, VaspCalculationConfig

    if isinstance(calculation, VaspCalculationConfig):
        incar = dict(calculation.incar or {})
        nested = next((value for key, value in incar.items() if str(key).lower() == "incar"), None)
        tags = {str(key).upper(): value for key, value in (nested if isinstance(nested, Mapping) else incar).items()}
        uses_u = _true(tags.get("LDAU", False))
        if _true(tags.get("LHFCALC", False)):
            return "hybrid", uses_u
        metagga = str(tags.get("METAGGA", "")).strip().upper()
        if metagga and metagga not in ("--", "NONE"):
            return _FUNCTIONAL_NAMES.get(metagga), uses_u
        gga = str(tags.get("GGA", "")).strip().upper()
        if gga:
            return _FUNCTIONAL_NAMES.get(gga), uses_u
        return ("PBE" if "PBE" in calculation.potential_family.upper() else None), uses_u
    if isinstance(calculation, QeCalculationConfig):
        system = {str(key).lower(): value for key, value in calculation.parameters.get("SYSTEM", {}).items()}
        uses_u = _true(system.get("lda_plus_u", False))
        dft = str(system.get("input_dft", "")).strip().upper().replace("-", "")
        if dft:
            return _FUNCTIONAL_NAMES.get(dft), uses_u
        family = calculation.pseudo_family.upper()
        if "PBESOL" in family:
            return "PBEsol", uses_u
        return ("PBE" if "PBE" in family else None), uses_u
    if isinstance(calculation, Cp2kCalculationConfig):
        uses_u = _has_key(calculation.parameters, "PLUS_U")
        xc = calculation.parameters["FORCE_EVAL"]["DFT"].get("XC", {})
        section = xc.get("XC_FUNCTIONAL", {}) if isinstance(xc, Mapping) else {}
        if isinstance(section, str):
            return _FUNCTIONAL_NAMES.get(section.strip().upper()), uses_u
        if "_" in section:
            return _FUNCTIONAL_NAMES.get(str(section["_"]).strip().upper()), uses_u
        named = [key for key in section if not str(key).startswith("_")]
        if len(named) != 1:
            return None, uses_u
        key = str(named[0]).upper()
        if key == "PBE":
            body = section[named[0]]
            parametrization = str(body.get("PARAMETRIZATION", "ORIG") if isinstance(body, Mapping) else "ORIG").upper()
            return {"ORIG": "PBE", "PBESOL": "PBEsol"}.get(parametrization), uses_u
        if key == "LIBXC" and "R2SCAN" in str(section[named[0]]).upper():
            return "r2SCAN", uses_u
        return _FUNCTIONAL_NAMES.get(key), uses_u
    raise TypeError("calculation must be a QE, VASP or CP2K recipe")


def foundation_level_warnings(model: str, head: str | None, calculation: Any, elements: Iterable[str]) -> list[str]:
    """What separates the DFT of a study from the training data of a foundation head (empty when nothing does).

    ``model`` and ``head`` name a head of :data:`MACE_FOUNDATION_MODELS`
    (the default head when ``head`` is ``None``), ``calculation`` is the DFT
    recipe of the statics and ``elements`` the elements of the study. Checks,
    in order: a head trained on molecules only; a different functional
    (nothing more is checked then); the Materials Project +U of the head on
    oxides or fluorides of :data:`MP_HUBBARD_METALS` in the study while the
    recipe has no +U; a head trained without +U while the recipe has +U. A
    recipe whose functional cannot be read gives no warning about it.
    """

    entry = mace_foundation_model(model)
    info = entry.head(head)
    alternatives = [(f"{other.name} (head {item.name!r})", item) for other in _MODELS for item in other.heads]
    return _level_warnings(f"{entry.name} (head {info.name!r})", info, calculation, elements, alternatives)


def _level_warnings(label: str, info: MaceFoundationHead, calculation: Any, elements: Iterable[str],
                    alternatives: Sequence[tuple[str, MaceFoundationHead]]) -> list[str]:
    """The checks of :func:`foundation_level_warnings` for the head ``info`` (``label`` in the messages), suggesting
    the matching heads among ``alternatives`` (label, head)."""

    messages = []
    if not info.for_slabs:
        messages.append(f"{label} was trained on molecules ({info.training_data}), not on solids or slabs")
    functional, uses_u = recipe_level_of_theory(calculation)
    if functional is None:
        return messages
    if functional != info.functional:
        closer = [name for name, item in alternatives if item.functional == functional and item.for_slabs]
        messages.append(
            f"the DFT recipe uses {functional} but {label} was trained on {info.functional} data; the fine-tune has "
            f"to learn the difference between the functionals from your labels"
            + (f" (trained on {functional}: {', '.join(closer)})" if closer else "")
        )
        return messages
    present = set(elements)
    hubbard_metals = sorted(present.intersection(MP_HUBBARD_METALS)) if present.intersection(MP_HUBBARD_ANIONS) else []
    if info.hubbard == "mp" and hubbard_metals and not uses_u:
        closer = [name for name, item in alternatives
                  if item.functional == functional and item.hubbard == "none" and item.for_slabs]
        messages.append(
            f"{label} was trained with the Materials Project +U on {', '.join(hubbard_metals)} in oxides and "
            f"fluorides, but the DFT recipe has no +U"
            + (f" (trained on {functional} without +U: {', '.join(closer)})" if closer else "")
        )
    if info.hubbard == "none" and uses_u:
        messages.append(f"{label} was trained without +U, but the DFT recipe uses +U")
    return messages

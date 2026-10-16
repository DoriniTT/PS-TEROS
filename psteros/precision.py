"""Precision levels of a VASP recipe: ``fast``, ``moderate`` and ``precise``.

A level fixes the numerical settings that decide the cost and the precision of a calculation, so that a recipe does not have to
spell them out: the plane-wave cutoff, the k-point spacing, the electronic convergence and, for a relaxation, the force
convergence.  Everything else a recipe needs (``ISIF``, ``NSW``, the code, the POTCAR family) is still given by the user, and
every value the user gives explicitly wins over the level.

aiida-vasp has a single VASP protocol, ``balanced`` (``moderate`` is an alias of it); it has no ``fast`` or ``stringent`` one.
PS-TEROS therefore defines the three levels itself, and ``moderate`` is the ``balanced`` protocol of aiida-vasp:

============  =========  ==========================  ==============  ============================
level         ``ENCUT``  ``kpoints_spacing`` (A^-1)  ``EDIFF`` (eV)  relaxation ``EDIFFG`` (eV/A)
============  =========  ==========================  ==============  ============================
``fast``      400        0.44  (0.07 x 2 pi)         1e-5            -0.05
``moderate``  550        0.31  (0.05 x 2 pi)         1e-6            -0.03
``precise``   700        0.19  (0.03 x 2 pi)         1e-7            -0.01
============  =========  ==========================  ==============  ============================

``moderate`` takes ``ENCUT``, the k-point spacing, the smearing and the electronic settings of aiida-vasp's ``balanced``;
``fast`` and ``precise`` follow the ends of the scan of aiida-vasp's convergence work chain (cutoff 300 to 700 eV, k-point spacing
0.07 to 0.02 in its units), starting ``fast`` at 400 eV, the cutoff of the oxygen POTCAR, and ``precise`` at 0.03.  The
force convergence of ``fast`` and ``moderate`` is that of aiida-vasp's relaxation protocols.  They are starting points, not a
substitute for a convergence test of the property you compute.

Differences from the ``balanced`` protocol, on purpose: the functional is PBE (``GGA = PE``; ``balanced`` sets PBEsol), the
calculation is not spin-polarised unless the structure holds an element that aiida-vasp gives a magnetic moment (then ``ISPIN = 2``
with its ``MAGMOM``), and ``EDIFF`` is the same for every structure instead of growing with the number of atoms, so that the
energies of cells of different sizes are comparable.  The POTCAR of an element that the recipe does not name is the one of the
aiida-vasp protocol (``Sn_d``, ``Sr_sv``, ``Ti_pv``, ``O``, ...).

Nothing here needs an AiiDA profile; only the POTCAR and magnetic-moment tables are read, from the protocol file of the
installed aiida-vasp.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Mapping, Sequence

TWO_PI = 2.0 * math.pi


@dataclass(frozen=True)
class PrecisionLevel:
    """One precision level.

    ``kpoints_spacing_aiida_vasp`` is the spacing in the unit of the aiida-vasp protocols (which multiply it by 2 pi);
    :attr:`kpoints_spacing` is the same value in the unit of PS-TEROS (A^-1 with the 2 pi of VASP's ``KSPACING``).
    """

    name: str
    encut_ev: float
    kpoints_spacing_aiida_vasp: float
    ediff_ev: float
    force_cutoff_ev_per_angstrom: float
    description: str = ""

    @property
    def kpoints_spacing(self) -> float:
        return self.kpoints_spacing_aiida_vasp * TWO_PI

    @property
    def ediffg_ev_per_angstrom(self) -> float:
        """``EDIFFG`` of a relaxation: negative, a force criterion."""

        return -self.force_cutoff_ev_per_angstrom


#: The three levels, by name.
PRECISION_LEVELS: Mapping[str, PrecisionLevel] = {
    "fast": PrecisionLevel("fast", 400.0, 0.07, 1.0e-5, 0.05, "quick checks and smoke tests"),
    "moderate": PrecisionLevel("moderate", 550.0, 0.05, 1.0e-6, 0.03, "the balanced protocol of aiida-vasp"),
    "precise": PrecisionLevel("precise", 700.0, 0.03, 1.0e-7, 0.01, "final numbers, at the cost of time"),
}

DEFAULT_PRECISION = "moderate"

#: Settings shared by the three levels (from the ``balanced`` protocol of aiida-vasp, with PBE).
COMMON_INCAR: Mapping[str, Any] = {
    "PREC": "Accurate",
    "ISMEAR": 0,
    "SIGMA": 0.05,
    "LREAL": False,
    "ALGO": "Normal",
    "NELM": 200,
    "LASPH": True,
    "LWAVE": False,
    "LCHARG": False,
    "GGA": "PE",
}


def precision_level(name: str) -> PrecisionLevel:
    """The level called ``name``; a ``ValueError`` that lists the known names otherwise."""

    try:
        return PRECISION_LEVELS[name]
    except (KeyError, TypeError):
        raise ValueError(f"precision must be one of {sorted(PRECISION_LEVELS)}, got {name!r}") from None


@lru_cache(maxsize=1)
def protocol_tables() -> dict[str, Any]:
    """The POTCAR and magnetic-moment tables of the ``balanced`` protocol of the installed aiida-vasp."""

    import importlib.util
    from pathlib import Path

    import yaml

    spec = importlib.util.find_spec("aiida_vasp")
    if spec is None or not spec.submodule_search_locations:
        raise ImportError("the precision levels read the POTCAR table of aiida-vasp, which is not installed")
    path = Path(next(iter(spec.submodule_search_locations))) / "protocols" / "vasp.yaml"
    with path.open() as handle:
        defaults = yaml.safe_load(handle)["default_inputs"]
    return {
        "potential_mapping": dict(defaults["potential_mapping"]),
        "magmom_mapping": dict(defaults["magmom_mapping"]),
        "protocol_encut": float(defaults["parameters"]["incar"]["encut"]),
        "protocol_kpoints_spacing": float(defaults["kpoints_spacing"]),
        "protocol_ediff_per_atom": float(defaults["meta_parameters"]["ediff_per_atom"]),
    }


def protocol_potential(symbol: str) -> str | None:
    """The POTCAR that the aiida-vasp protocol uses for an element, or ``None``."""

    return protocol_tables()["potential_mapping"].get(symbol)


def spin_incar(site_symbols: Sequence[str]) -> dict[str, Any]:
    """``ISPIN`` (and ``MAGMOM``): off, unless an element has a magnetic moment in the aiida-vasp protocol."""

    moments = protocol_tables()["magmom_mapping"]
    default = float(moments.get("default", 1.0))
    if not any(symbol in moments and symbol != "default" for symbol in site_symbols):
        return {"ISPIN": 1}
    return {"ISPIN": 2, "MAGMOM": [float(moments.get(symbol, default)) for symbol in site_symbols]}


def structure_symbols(structure: Any) -> list[str]:
    """Element symbol of every site of a pymatgen ``Structure`` or an AiiDA ``StructureData``."""

    if hasattr(structure, "get_kind"):
        return [structure.get_kind(site.kind_name).symbol for site in structure.sites]
    return [site.specie.symbol for site in structure]


def is_relaxation(incar: Mapping[str, Any]) -> bool:
    """Whether an INCAR relaxes the ions (``IBRION`` 1, 2 or 3 with ``NSW > 0``)."""

    values = {str(key).upper(): value for key, value in incar.items()}
    return int(values.get("IBRION", -1)) in (1, 2, 3) and int(values.get("NSW", 0)) > 0


def level_incar(
    name: str, user_incar: Mapping[str, Any], *, site_symbols: Sequence[str] | None = None
) -> dict[str, Any]:
    """The INCAR tags that level ``name`` adds to ``user_incar``: everything it sets that the user did not.

    ``EDIFFG`` is added only when ``user_incar`` relaxes the ions, and ``ISPIN`` (with ``MAGMOM`` if magnetic) only when
    ``site_symbols`` is given.  The keys are upper case; a tag the user wrote in any case is left alone.
    """

    level = precision_level(name)
    given = {str(key).upper() for key in user_incar}
    tags: dict[str, Any] = {**COMMON_INCAR, "ENCUT": level.encut_ev, "EDIFF": level.ediff_ev}
    if is_relaxation(user_incar):
        tags["EDIFFG"] = level.ediffg_ev_per_angstrom
    if site_symbols is not None and not {"ISPIN", "MAGMOM"} & given:
        tags.update(spin_incar(site_symbols))
    return {key: value for key, value in tags.items() if key not in given}

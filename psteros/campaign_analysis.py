"""Energies of a campaign as plain Python, for any material.

A campaign graph (:func:`psteros.campaign.build_vasp_campaign_workgraph`) has
the same layout for every material: references and slabs, by label and block.
This module turns its results into analysis inputs without assuming what the
material is:

* :class:`CampaignEntry` is one structure of a campaign: its group
  (``references`` or ``slabs``), phase, composition, energy (eV) and, for a
  slab, the area of one face (A^2);
* :func:`campaign_chemical_potentials` gives the energy per atom of every
  element that has a single-element reference (its elemental limit);
* :func:`campaign_references` builds the reference objects of the binary and
  ternary oxide phase diagrams from the compositions of the references;
* :func:`campaign_surface_energies` gives gamma of every slab at given
  chemical potentials, which is the whole answer for a unary system.

Every function takes either the PK of a campaign graph or entries typed by
hand, so nothing here needs AiiDA unless a PK is given.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Iterable, Mapping

from psteros.phase_diagram import BinaryOxideReferences, _integer_composition
from psteros.phase_diagram_ternary import CompetingPhase, TernaryOxideReferences
from psteros.thermodynamics import surface_energy_elemental

GROUPS = ("references", "slabs")
PHASES = ("gas", "solid")


@dataclass(frozen=True)
class CampaignEntry:
    """One structure of a campaign and its energy.

    ``energy_ev`` is the energy of the whole cell (eV) from ``block``, or
    ``None`` while that block has not finished (``state`` says why).
    ``surface_area_angstrom2`` is the area of one face of a slab (A^2), the
    plane of its first two lattice vectors; references have none.
    """

    label: str
    group: str
    phase: str
    composition: Mapping[str, int]
    energy_ev: float | None
    block: str = "static"
    state: str = "finished"
    surface_area_angstrom2: float | None = None

    def __post_init__(self) -> None:
        if not self.label:
            raise ValueError("entry label must not be empty")
        if self.group not in GROUPS:
            raise ValueError(f"{self.label}: group must be one of {GROUPS}, got {self.group!r}")
        if self.phase not in PHASES:
            raise ValueError(f"{self.label}: phase must be one of {PHASES}, got {self.phase!r}")
        object.__setattr__(self, "composition", _integer_composition(self.composition, f"{self.label}: composition"))
        if self.energy_ev is None:
            if self.state == "finished":
                raise ValueError(f"{self.label}: an entry without energy_ev cannot be in state 'finished'")
        else:
            object.__setattr__(self, "energy_ev", _finite(self.label, "energy_ev", self.energy_ev))
        if self.group == "slabs" and self.phase != "solid":
            raise ValueError(f"{self.label}: a slab is a solid, got phase {self.phase!r}")
        if self.surface_area_angstrom2 is not None:
            if self.group != "slabs":
                raise ValueError(f"{self.label}: only slabs have a surface area")
            area = _finite(self.label, "surface_area_angstrom2", self.surface_area_angstrom2)
            object.__setattr__(self, "surface_area_angstrom2", area)
            if area <= 0:
                raise ValueError(f"{self.label}: surface_area_angstrom2 must be a positive number")

    @property
    def atoms(self) -> int:
        return sum(self.composition.values())

    @property
    def energy_per_atom_ev(self) -> float | None:
        return None if self.energy_ev is None else self.energy_ev / self.atoms


def _finite(label: str, name: str, value: Any) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        number = math.nan
    if not math.isfinite(number):
        raise ValueError(f"{label}: {name} must be a finite number, got {value!r}")
    return number


def _entries(source: Any) -> tuple[CampaignEntry, ...]:
    """Entries of a campaign PK, or the given entries."""

    if isinstance(source, int) and not isinstance(source, bool):
        from psteros.campaign import campaign_entries

        return campaign_entries(source)
    entries = tuple(source)
    for entry in entries:
        if not isinstance(entry, CampaignEntry):
            raise TypeError(f"source must be a campaign PK or CampaignEntry objects, got {type(entry).__name__}")
    labels = [entry.label for entry in entries]
    repeated = sorted({label for label in labels if labels.count(label) > 1})
    if repeated:
        raise ValueError(f"labels {repeated} appear more than once")
    return entries


def _energy(entry: CampaignEntry) -> float:
    if entry.energy_ev is None:
        raise ValueError(f"{entry.label}: block {entry.block!r} has no energy yet ({entry.state})")
    return entry.energy_ev


def _is_o2(entry: CampaignEntry) -> bool:
    return entry.group == "references" and entry.phase == "gas" and dict(entry.composition) == {"O": 2}


def _check_labels(entries: Iterable[CampaignEntry], labels: Iterable[str], option: str) -> None:
    known = {entry.label for entry in entries}
    unknown = sorted(set(labels).difference(known))
    if unknown:
        raise ValueError(f"{option} names unknown labels {unknown}")


def _reservoirs(
    entries: tuple[CampaignEntry, ...],
    reservoirs: Mapping[str, str] | None,
    elements: Iterable[str] | None = None,
) -> dict[str, CampaignEntry]:
    """``{element: reference entry}`` of the elemental reservoirs, for ``elements`` (default: all).

    The candidates of an element are its single-element references; for
    oxygen, the O2 gas references when there are any (an O atom does not
    compete with O2).  Several candidates are an error unless ``reservoirs``
    chooses one.
    """

    if isinstance(elements, str):
        raise TypeError(f"elements must be a collection of element symbols, e.g. ({elements!r},), not a string")
    chosen = dict(reservoirs or {})
    _check_labels(entries, chosen.values(), "reservoirs")
    by_label = {entry.label: entry for entry in entries}
    for element, label in chosen.items():
        entry = by_label[label]
        if entry.group != "references" or set(entry.composition) != {element}:
            raise ValueError(f"reservoirs: {label!r} is not a single-element reference of {element}")
    candidates: dict[str, list[CampaignEntry]] = {}
    for entry in entries:
        if entry.group == "references" and len(entry.composition) == 1:
            candidates.setdefault(next(iter(entry.composition)), []).append(entry)
    if any(_is_o2(entry) for entry in candidates.get("O", ())):
        candidates["O"] = [entry for entry in candidates["O"] if _is_o2(entry)]
    wanted = set(candidates) if elements is None else set(elements)
    result = {}
    for element in sorted(wanted.intersection(candidates)):
        found = candidates[element]
        if element in chosen:
            result[element] = by_label[chosen[element]]
        elif len(found) > 1:
            raise ValueError(
                f"several references of {element}: {sorted(entry.label for entry in found)}; "
                f"choose one with reservoirs={{{element!r}: label}}"
            )
        else:
            result[element] = found[0]
    return result


def campaign_chemical_potentials(
    source: Any, *, reservoirs: Mapping[str, str] | None = None, elements: Iterable[str] | None = None
) -> dict[str, float]:
    """``{element: eV per atom}`` of the elements that have a single-element reference.

    mu_X = E(reference of X) / (atoms in it): a solid gives its energy per
    atom (a bulk metal), a gas its energy per atom (E(O2)/2 for oxygen).
    These are the elemental limits of the chemical potentials.  The
    reference of oxygen is an O2 gas when there is one.  When an element has
    several references, ``reservoirs={"Sn": "sn_beta"}`` chooses one.
    ``elements`` limits the result (and the checks) to those elements.
    Multi-element references are not used here.
    """

    entries = _entries(source)
    found = _reservoirs(entries, reservoirs, elements)
    missing = sorted(set(elements or ()).difference(found))
    if missing:
        raise ValueError(f"no single-element reference of {missing}")
    return {element: _energy(entry) / entry.atoms for element, entry in found.items()}


def campaign_references(
    source: Any,
    *,
    host: str,
    reservoirs: Mapping[str, str] | None = None,
    exclude: Iterable[str] = (),
    independent: str | None = None,
) -> BinaryOxideReferences | TernaryOxideReferences:
    """Reference object of the surface phase diagram of the bulk ``host``.

    The host's elements decide the model: ``{M, O}`` gives
    :class:`~psteros.phase_diagram.BinaryOxideReferences` (the metal
    reference is optional), ``{A, B, O}`` gives
    :class:`~psteros.phase_diagram_ternary.TernaryOxideReferences`, whose
    competing phases are the other solid references made of the host's
    elements (compounds such as SrO and TiO2, and elemental phases other than
    the chosen reservoirs).  Oxygen comes from an O2 gas reference.
    References that do not fit the model (other elements, compound gases,
    any other solid for a binary oxide) are an error unless listed in
    ``exclude``.  Other hosts have no phase-diagram model yet: use
    :func:`campaign_chemical_potentials` and :func:`campaign_surface_energies`.
    """

    try:
        return _campaign_references(_entries(source), host, reservoirs, tuple(exclude), independent)
    except ValueError as error:
        message = str(error)
        if message.startswith((f"{host!r}", f"{host}:", f"host {host!r}")):
            raise
        raise ValueError(f"{host!r}: {message}") from error


def _campaign_references(
    entries: tuple[CampaignEntry, ...],
    host: str,
    reservoirs: Mapping[str, str] | None,
    exclude: tuple[str, ...],
    independent: str | None,
) -> BinaryOxideReferences | TernaryOxideReferences:
    _check_labels(entries, exclude, "exclude")
    _check_labels(entries, [host], "host")
    if host in exclude:
        raise ValueError(f"host {host!r} cannot be excluded")
    excluded_reservoirs = sorted(set((reservoirs or {}).values()).intersection(exclude))
    if excluded_reservoirs:
        raise ValueError(f"reservoirs {excluded_reservoirs} are also in exclude")
    entries = tuple(entry for entry in entries if entry.label not in exclude)
    bulk = next(entry for entry in entries if entry.label == host)
    if bulk.group != "references" or bulk.phase != "solid":
        raise ValueError(f"host {host!r} must be a solid reference, got a {bulk.phase} in {bulk.group}")
    elements = set(bulk.composition)
    metals = sorted(elements - {"O"})
    if "O" not in elements or len(metals) not in (1, 2):
        formula = "".join(f"{element}{count}" for element, count in sorted(bulk.composition.items()))
        raise ValueError(
            f"no phase-diagram model for {host!r} ({formula}) yet: only binary and ternary oxides have one; "
            "use campaign_chemical_potentials and campaign_surface_energies"
        )
    if len(metals) == 1 and independent is not None:
        raise ValueError(f"independent applies to ternary oxides; {host!r} is a binary oxide")

    references = [entry for entry in entries if entry.group == "references" and entry.label != host]
    foreign = sorted(entry.label for entry in references if not set(entry.composition) <= elements)
    if foreign:
        raise ValueError(f"references {foreign} contain elements outside {host!r}; list them in exclude")
    compound_gases = sorted(entry.label for entry in references if entry.phase == "gas" and len(entry.composition) > 1)
    if compound_gases:
        raise ValueError(
            f"compound gas references {compound_gases} are not reservoirs of this model; list them in exclude"
        )

    found = _reservoirs(entries, reservoirs, elements)
    oxygen = found.get("O")
    if oxygen is None or not _is_o2(oxygen):
        gases = sorted(entry.label for entry in references if entry.phase == "gas")
        raise ValueError(f"{host!r}: the oxygen reference must be an O2 gas reference; gas references: {gases}")
    reservoir_labels = {entry.label for entry in found.values()}
    others = [entry for entry in references if entry.phase == "solid" and entry.label not in reservoir_labels]
    if len(metals) == 1:
        if others:
            raise ValueError(
                f"references {sorted(entry.label for entry in others)} are other solids, which a binary "
                "oxide diagram does not use; list them in exclude"
            )
        metal = found.get(metals[0])
        return BinaryOxideReferences(
            bulk_energy_ev=_energy(bulk),
            bulk_composition=dict(bulk.composition),
            oxygen_molecule_energy_ev=_energy(oxygen),
            metal_energy_per_atom_ev=None if metal is None else _energy(metal) / metal.atoms,
        )
    missing = [element for element in metals if element not in found]
    if missing:
        raise ValueError(f"a ternary oxide needs single-element references of {missing}")
    return TernaryOxideReferences(
        bulk_energy_ev=_energy(bulk),
        bulk_composition=dict(bulk.composition),
        oxygen_molecule_energy_ev=_energy(oxygen),
        element_energies_per_atom_ev={
            element: _energy(found[element]) / found[element].atoms for element in metals
        },
        competing_phases=tuple(
            CompetingPhase(label=entry.label, energy_ev=_energy(entry), composition=dict(entry.composition))
            for entry in others
        ),
        independent=independent,
    )


def campaign_surface_energies(
    source: Any,
    *,
    chemical_potentials_ev: Mapping[str, float] | None = None,
    reservoirs: Mapping[str, str] | None = None,
    surfaces: int = 2,
) -> dict[str, float]:
    """``{slab label: gamma in eV/A^2}`` at the chemical potentials ``chemical_potentials_ev`` (eV per atom).

    gamma = (E_slab - sum_i N_i mu_i) / (surfaces * A), with A the area of one
    face.  Without ``chemical_potentials_ev``, mu is the energy per atom of
    the element's reference (:func:`campaign_chemical_potentials`), which is
    exact for a unary system; slabs of several elements then need explicit
    chemical potentials, because the elemental limits are not a state in
    equilibrium with their bulk.  Multiply by
    :data:`psteros.EV_PER_ANGSTROM2_TO_J_PER_M2` for J/m^2.
    """

    entries = _entries(source)
    slabs = [entry for entry in entries if entry.group == "slabs"]
    if chemical_potentials_ev is None:
        compounds = sorted(entry.label for entry in slabs if len(entry.composition) > 1)
        if compounds:
            raise ValueError(
                f"slabs {compounds} contain several elements: pass chemical_potentials_ev (eV per atom) at the "
                "conditions of interest; the elemental limits are not in equilibrium with their bulk"
            )
        needed = {element for entry in slabs for element in entry.composition}
        chemical_potentials_ev = {
            element: _energy(entry) / entry.atoms for element, entry in _reservoirs(entries, reservoirs, needed).items()
        }
    elif reservoirs is not None:
        raise ValueError("reservoirs choose the default chemical potentials; give it or chemical_potentials_ev, not both")
    gammas = {}
    for entry in slabs:
        missing = sorted(set(entry.composition).difference(chemical_potentials_ev))
        if missing:
            raise ValueError(f"{entry.label}: no chemical potential for {missing}")
        if entry.surface_area_angstrom2 is None:
            raise ValueError(f"{entry.label}: block {entry.block!r} has no surface area yet ({entry.state})")
        gammas[entry.label] = surface_energy_elemental(
            slab_energy_ev=_energy(entry),
            stoichiometry=dict(entry.composition),
            chemical_potentials_ev=chemical_potentials_ev,
            surface_area_angstrom2=entry.surface_area_angstrom2,
            surfaces=surfaces,
        )
    return gammas

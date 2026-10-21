"""Phonons on a q-point mesh with phonopy, for the vibrational free energies (no AiiDA).

The Gamma-point finite differences of :mod:`psteros.vibrations` see a crystal
only at the q-points commensurate with the cell of the calculation. Here the
force constants of a supercell (phonopy, symmetry-reduced displacements) are
Fourier-interpolated onto a q-point mesh of ``N_q`` points, and the
vibrational free energy per input cell is the mesh average

``F_vib(T) = 1/N_q sum_q sum_nu [h nu_qnu / 2 + k_B T ln(1 - exp(-h nu_qnu / k_B T))]``.

The mesh frequencies go into a :class:`psteros.HarmonicVibrations` with
``supercell_size = N_q``, so the free energies, :func:`psteros.solid_free_energy_ev`
and :func:`psteros.free_energies_with_vibrations` are those of the Gamma-point
modes, unchanged.

**Slabs vibrate as a whole** (full-slab phonons): every atom of the slab, in
an in-plane supercell ``(na, nb, 1)`` with a 2D mesh ``(ma, mb, 1)``, the
vacuum along ``c``. The frozen centre of a relaxation is only a constraint of
that relaxation; no region is counted as bulk, and for a stoichiometric slab
``gamma_vib = [F_vib(slab) - n F_vib(bulk)] / 2A``. Gas molecules are not
handled here: they keep their Gamma-point modes.

**Displacements and forces.** Each displacement moves one atom of the
supercell (phonopy's symmetry-reduced set; sites of one element with
different kind names are kept apart). The forces of the undisplaced supercell
are subtracted from those of every displacement, so that the residual forces
of a geometry that is not the minimum of the force field (DFT forces at an
MLIP geometry) do not enter the force constants.

**The mesh** is shifted by half a step (Monkhorst-Pack) on every axis with
more than one point: Gamma, whose three acoustic modes have zero frequency,
is not sampled, and F_vib converges much faster with the mesh than on a
Gamma-centred one (``shifted_mesh=False``, which leaves the 3 acoustic modes
at Gamma out).

**Stability.** The frequencies on the mesh and at the q-points commensurate
with the supercell (where the force constants are exact, Gamma included) are
searched for imaginary modes: :class:`PhononStability`.

phonopy is an optional dependency (``pip install "psteros[phonopy]"``).
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from itertools import product
from typing import Any, Mapping, Sequence

from psteros.vibrations import HarmonicVibrations, ImaginaryModes

#: One terahertz in cm^-1 (phonopy's frequencies are in THz).
THZ_TO_CM1 = 1.0e12 / 2.997_924_58e10

#: Smallest gap (A) along ``c`` that counts as the vacuum of a slab.
SLAB_VACUUM_ANGSTROM = 5.0


def _triple(value: Any, name: str) -> tuple[int, int, int]:
    value = tuple(value)
    if len(value) != 3 or any(not isinstance(n, int) or isinstance(n, bool) or n <= 0 for n in value):
        raise ValueError(f"{name} must be three positive integers, got {value!r}")
    return value


def _lattice(structure: Any) -> Any:
    import numpy as np

    return np.asarray(structure.lattice.matrix, dtype=float)


def _widths(lattice: Any) -> Any:
    """Distance between opposite faces of the cell, per axis (A): ``1 / |b_i|``, ``b`` without the 2 pi."""

    import numpy as np

    return 1.0 / np.linalg.norm(np.linalg.inv(lattice).T, axis=1)


def vacuum_along_c_angstrom(structure: Any) -> float:
    """The widest gap between the atoms of ``structure`` along ``c`` (A, perpendicular to the ``ab`` plane)."""

    import numpy as np

    heights = np.sort(np.mod(np.asarray(structure.frac_coords, dtype=float)[:, 2], 1.0))
    gaps = np.append(np.diff(heights), 1.0 + heights[0] - heights[-1])
    return float(gaps.max() * _widths(_lattice(structure))[2])


@dataclass(frozen=True)
class PhonopyConfig:
    """Phonons on a q-point mesh with phonopy: the opt-in mode of :class:`psteros.VibrationsConfig`.

    Every structure that is not a molecule is calculated with phonopy (the
    molecules keep their Gamma-point modes).

    ``supercell_length_angstrom``
        The supercell of each structure is the smallest diagonal one
        ``(n1, n2, n3)`` whose faces are at least this far apart on every
        periodic axis (``n_i = ceil(L / w_i)``, ``w_i`` the width of the cell
        along axis ``i``); a slab has ``n3 = 1``.
    ``mesh_length_angstrom``
        The q-point mesh has ``m_i = max(1, round(L |b_i|))`` points along
        axis ``i`` (``b_i`` the reciprocal vector without 2 pi, the "length"
        of phonopy and VASP); a slab has ``m3 = 1``.
    ``slabs``
        Labels of the slabs (vacuum along ``c``): in-plane supercell and 2D
        mesh. A slab needs a gap of at least ``SLAB_VACUUM_ANGSTROM`` along
        ``c``; a structure with such a gap that is not listed gets a warning.
    ``supercells``, ``meshes``
        Per label, an explicit supercell or mesh in place of the lengths (a
        slab's with 1 along ``c``).
    ``shifted_mesh``
        Shift the mesh by half a step on every axis with more than one point
        (default), or use a Gamma-centred mesh (``False``).
    ``plus_minus``
        Displace every atom by ``+d`` and ``-d`` (central differences); by
        default phonopy adds ``-d`` only where symmetry does not give it.
    ``symprec``
        Symmetry tolerance of phonopy (A). Raising it (1e-3) gives fewer
        displacements for a structure that is only nearly symmetric.
    """

    supercell_length_angstrom: float = 10.0
    mesh_length_angstrom: float = 40.0
    slabs: tuple[str, ...] = ()
    supercells: Mapping[str, tuple[int, int, int]] = field(default_factory=dict)
    meshes: Mapping[str, tuple[int, int, int]] = field(default_factory=dict)
    shifted_mesh: bool = True
    plus_minus: bool = False
    symprec: float = 1e-5

    def __post_init__(self) -> None:
        for name in ("supercell_length_angstrom", "mesh_length_angstrom", "symprec"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
                raise ValueError(f"PhonopyConfig.{name} must be a positive number, got {value!r}")
        slabs = tuple(self.slabs)
        if any(not isinstance(label, str) for label in slabs) or len(set(slabs)) != len(slabs):
            raise ValueError("PhonopyConfig.slabs must be distinct labels")
        object.__setattr__(self, "slabs", slabs)
        for name in ("supercells", "meshes"):
            values = {label: _triple(size, f"PhonopyConfig.{name}[{label!r}]") for label, size in dict(getattr(self, name)).items()}
            wrong = sorted(label for label, size in values.items() if label in slabs and size[2] != 1)
            if wrong:
                raise ValueError(f"PhonopyConfig.{name}: a slab has 1 along c (vacuum), got {', '.join(wrong)}")
            object.__setattr__(self, name, values)
        for name in ("shifted_mesh", "plus_minus"):
            if not isinstance(getattr(self, name), bool):
                raise ValueError(f"PhonopyConfig.{name} must be True or False")

    @property
    def labels(self) -> set[str]:
        """The labels that this config names."""

        return set(self.slabs) | set(self.supercells) | set(self.meshes)

    def check(self, label: str, structure: Any) -> None:
        """A listed slab needs its vacuum along ``c``; a structure with vacuum that is not listed gets a warning."""

        vacuum = vacuum_along_c_angstrom(structure)
        if label in self.slabs and vacuum < SLAB_VACUUM_ANGSTROM:
            raise ValueError(
                f"{label}: PhonopyConfig.slabs lists it, but it has no vacuum along c (widest gap {vacuum:.1f} A)"
            )
        if label not in self.slabs and vacuum >= SLAB_VACUUM_ANGSTROM:
            warnings.warn(
                f"{label}: a gap of {vacuum:.1f} A along c, but PhonopyConfig.slabs does not list it; it is sampled "
                "as a bulk (supercell and q-mesh along c too)", UserWarning, stacklevel=3,
            )

    def supercell(self, label: str, structure: Any) -> tuple[int, int, int]:
        """The phonopy supercell of ``label`` (explicit, or from ``supercell_length_angstrom``)."""

        if label in self.supercells:
            return self.supercells[label]
        sizes = [max(1, math.ceil(self.supercell_length_angstrom / width - 1e-9)) for width in _widths(_lattice(structure))]
        if label in self.slabs:
            sizes[2] = 1
        return tuple(int(n) for n in sizes)

    def mesh(self, label: str, structure: Any) -> tuple[int, int, int]:
        """The q-point mesh of ``label`` (explicit, or from ``mesh_length_angstrom``)."""

        if label in self.meshes:
            return self.meshes[label]
        sizes = [max(1, int(round(self.mesh_length_angstrom / width))) for width in _widths(_lattice(structure))]
        if label in self.slabs:
            sizes[2] = 1
        return tuple(int(n) for n in sizes)


def _kind_moments(structure: Any) -> list[float] | None:
    """Fake collinear moments (1, 2, ...) that keep sites of one element with different kind names apart."""

    kinds = structure.site_properties.get("kind_name")
    if not kinds:
        return None
    names = sorted(set(kinds))
    if len(names) == len({site.specie.symbol for site in structure}):
        return None
    return [1.0 + names.index(kind) for kind in kinds]


def _phonopy(structure: Any, supercell: Sequence[int], *, symprec: float, masses_amu: Sequence[float] | None) -> Any:
    try:
        from phonopy import Phonopy
        from phonopy.structure.atoms import PhonopyAtoms
    except ImportError as error:
        raise ImportError('phonons on a q-point mesh need phonopy: pip install "psteros[phonopy]"') from error
    import numpy as np

    if masses_amu is None:
        masses_amu = [float(site.specie.atomic_mass) for site in structure]
    if len(masses_amu) != len(structure) or any(mass <= 0 for mass in masses_amu):
        raise ValueError("masses_amu must be one positive mass per site")
    unit = PhonopyAtoms(
        symbols=[site.specie.symbol for site in structure], cell=_lattice(structure),
        scaled_positions=np.asarray(structure.frac_coords, dtype=float), masses=[float(m) for m in masses_amu],
        magnetic_moments=_kind_moments(structure),
    )
    with warnings.catch_warnings():
        # an anisotropic supercell has fewer point operations than the cell: expected, phonopy uses the supercell's
        warnings.filterwarnings("ignore", message=".*Point group symmetries of supercell and primitive", category=UserWarning)
        return Phonopy(unit, supercell_matrix=np.diag(_triple(supercell, "supercell")), primitive_matrix=np.eye(3),
                       symprec=symprec, log_level=0)


@dataclass(frozen=True)
class PhonopyDisplacements:
    """The displacements of one structure in its phonopy supercell (see :func:`phonopy_displacements`).

    ``lattice_angstrom`` and ``positions_angstrom`` are the cell (rows) and the
    Cartesian positions of the undisplaced supercell in phonopy's order,
    ``unit_sites`` the site of the input structure that each supercell atom
    repeats, and each of ``displacements`` moves one supercell atom by a
    Cartesian vector (A).
    """

    supercell: tuple[int, int, int]
    lattice_angstrom: tuple[tuple[float, float, float], ...]
    positions_angstrom: tuple[tuple[float, float, float], ...]
    unit_sites: tuple[int, ...]
    displacements: tuple[tuple[int, tuple[float, float, float]], ...]
    displacement_angstrom: float

    def __len__(self) -> int:
        return len(self.displacements)

    def displaced_positions(self, index: int) -> Any:
        """Cartesian positions (A) of the supercell of displacement ``index``."""

        import numpy as np

        positions = np.array(self.positions_angstrom, dtype=float)
        atom, vector = self.displacements[index]
        positions[atom] += vector
        return positions

    def structures(self, structure: Any) -> tuple[Any, list[Any]]:
        """The undisplaced supercell and the displaced ones as pymatgen structures, with the species and site
        properties (kind names) of the sites of ``structure`` they repeat."""

        from pymatgen.core import Structure

        species = [structure[site].species for site in self.unit_sites]
        properties = {key: [values[site] for site in self.unit_sites] for key, values in structure.site_properties.items()}

        def build(positions):
            return Structure(self.lattice_angstrom, species, positions, coords_are_cartesian=True,
                             site_properties=properties or None)

        return build(self.positions_angstrom), [build(self.displaced_positions(i)) for i in range(len(self))]

    def to_dict(self) -> dict[str, Any]:
        return {
            "supercell": list(self.supercell),
            "lattice_angstrom": [list(row) for row in self.lattice_angstrom],
            "positions_angstrom": [list(row) for row in self.positions_angstrom],
            "unit_sites": list(self.unit_sites),
            "displacements": [[atom, *vector] for atom, vector in self.displacements],
            "displacement_angstrom": self.displacement_angstrom,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "PhonopyDisplacements":
        if hasattr(data, "get_dict"):
            data = data.get_dict()
        return cls(
            supercell=tuple(int(n) for n in data["supercell"]),
            lattice_angstrom=tuple(tuple(float(x) for x in row) for row in data["lattice_angstrom"]),
            positions_angstrom=tuple(tuple(float(x) for x in row) for row in data["positions_angstrom"]),
            unit_sites=tuple(int(site) for site in data["unit_sites"]),
            displacements=tuple((int(row[0]), tuple(float(x) for x in row[1:4])) for row in data["displacements"]),
            displacement_angstrom=float(data["displacement_angstrom"]),
        )


def phonopy_displacements(
    structure: Any,
    supercell: Sequence[int],
    *,
    displacement_angstrom: float = 0.01,
    plus_minus: bool = False,
    symprec: float = 1e-5,
) -> PhonopyDisplacements:
    """phonopy's symmetry-reduced displacements of ``structure`` (pymatgen) in the supercell ``(n1, n2, n3)``.

    Each displacement moves one atom of the supercell by ``displacement_angstrom``;
    ``plus_minus`` adds the opposite displacement of each (otherwise phonopy
    adds it only where symmetry does not give it). Sites of one element with
    different ``kind_name`` site properties are not treated as equivalent.
    """

    if not 0 < displacement_angstrom <= 0.1:
        raise ValueError("displacement_angstrom must be in (0, 0.1] A")
    phonon = _phonopy(structure, supercell, symprec=symprec, masses_amu=None)
    phonon.generate_displacements(distance=float(displacement_angstrom), is_plusminus=True if plus_minus else "auto")
    cell = phonon.supercell
    unit_sites = [int(cell.u2u_map[cell.s2u_map[index]]) for index in range(len(cell))]
    return PhonopyDisplacements(
        supercell=_triple(supercell, "supercell"),
        lattice_angstrom=tuple(tuple(float(x) for x in row) for row in cell.cell),
        positions_angstrom=tuple(tuple(float(x) for x in row) for row in cell.positions),
        unit_sites=tuple(unit_sites),
        displacements=tuple((int(item["number"]), tuple(float(x) for x in item["displacement"]))
                            for item in phonon.dataset["first_atoms"]),
        displacement_angstrom=float(displacement_angstrom),
    )


def _grid(mesh: Sequence[int], shifted: bool) -> Any:
    """Fractional q-points of a mesh, shifted by half a step on the axes with more than one point."""

    import numpy as np

    axes = [[(i + (0.5 if shifted and n > 1 else 0.0)) / n for i in range(n)] for n in mesh]
    return np.array(list(product(*axes)), dtype=float)


@dataclass(frozen=True)
class PhononStability:
    """Imaginary modes of a phonopy calculation: a dynamical stability check.

    ``q_points`` is the number of q-points searched (the mesh and the points
    commensurate with the supercell), ``imaginary`` the ``(q, lowest
    frequency)`` of each with an imaginary mode (fractional q, cm^-1 < 0) and
    ``lowest_frequency_cm1`` the lowest frequency found (the 3 acoustic modes
    at Gamma left out). An imaginary mode means the structure is not a minimum
    (a saddle point of the relaxation, or a distortion the cell cannot follow).
    """

    q_points: int
    imaginary: tuple[tuple[tuple[float, float, float], float], ...]
    lowest_frequency_cm1: float

    @property
    def stable(self) -> bool:
        return not self.imaginary

    def summary(self) -> str:
        if self.stable:
            return f"stable: no imaginary mode at {self.q_points} q-points (lowest {self.lowest_frequency_cm1:.1f} cm^-1)"
        worst = min(self.imaginary, key=lambda item: item[1])
        q = ", ".join(f"{x:.3f}" for x in worst[0])
        return (f"NOT STABLE: imaginary modes at {len(self.imaginary)} of {self.q_points} q-points, "
                f"down to {-worst[1]:.1f}i cm^-1 at q = ({q})")

    def to_dict(self) -> dict[str, Any]:
        return {"q_points": self.q_points, "imaginary": [[list(q), value] for q, value in self.imaginary],
                "lowest_frequency_cm1": self.lowest_frequency_cm1}

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "PhononStability":
        if hasattr(data, "get_dict"):
            data = data.get_dict()
        if "stability" in data:
            data = data["stability"]
        return cls(int(data["q_points"]), tuple((tuple(float(x) for x in q), float(value)) for q, value in data["imaginary"]),
                   float(data["lowest_frequency_cm1"]))


@dataclass(frozen=True)
class PhononModes:
    """The result of :func:`phonopy_modes`: the mesh modes and the stability check.

    ``vibrations`` holds the frequencies of every q-point of the mesh
    (``supercell_size`` is the number of q-points), so its free energies are
    per input cell; :meth:`to_dict` is the ``<label>_vibrations`` output that
    :func:`psteros.read_vibrations` reads.
    """

    vibrations: HarmonicVibrations
    stability: PhononStability
    supercell: tuple[int, int, int]
    mesh: tuple[int, int, int]
    shifted_mesh: bool
    zero_modes: int
    displacements: int
    displacement_angstrom: float
    symprec: float

    def to_dict(self) -> dict[str, Any]:
        return {
            **self.vibrations.to_dict(),
            "method": "phonopy",
            "supercell": list(self.supercell),
            "mesh": list(self.mesh),
            "shifted_mesh": self.shifted_mesh,
            "zero_modes": self.zero_modes,
            "displacements": self.displacements,
            "displacement_angstrom": self.displacement_angstrom,
            "symprec": self.symprec,
            "stability": self.stability.to_dict(),
        }


def phonopy_modes(
    structure: Any,
    displacements: PhonopyDisplacements,
    forces_ev_per_angstrom: Sequence[Any],
    *,
    reference_forces_ev_per_angstrom: Sequence[Any] | None,
    mesh: Sequence[int],
    shifted_mesh: bool = True,
    symprec: float = 1e-5,
    masses_amu: Sequence[float] | None = None,
    imaginary_modes: ImaginaryModes = "raise",
    low_frequency_cutoff_cm1: float | None = None,
) -> PhononModes:
    """Force constants from the forces of the displacements, and the modes on a q-point mesh.

    ``forces_ev_per_angstrom[k]`` are the forces (eV/A, one row per supercell
    atom) of displacement ``k`` of ``displacements``, and
    ``reference_forces_ev_per_angstrom`` those of the undisplaced supercell,
    subtracted from each (``None`` only if they are zero by construction).
    ``mesh`` is the q-point mesh, shifted by half a step unless
    ``shifted_mesh`` is false (a Gamma-centred mesh loses the 3 acoustic
    modes at Gamma). ``masses_amu`` are those of the sites of ``structure``
    (by default pymatgen's). The vibrations raise on imaginary modes unless
    ``imaginary_modes="drop"``; the stability check lists them in any case.

    The force constants are phonopy's (symmetrised), the frequencies those of
    its dynamical matrix at each q of the mesh, and the free energy of the
    result is per input cell,
    ``F_vib(T) = 1/N_q sum_q sum_nu [h nu / 2 + k_B T ln(1 - exp(-h nu / k_B T))]``.
    """

    import numpy as np

    mesh = _triple(mesh, "mesh")
    phonon = _phonopy(structure, displacements.supercell, symprec=symprec, masses_amu=masses_amu)
    cell = phonon.supercell
    if (len(cell) != len(displacements.positions_angstrom)
            or not np.allclose(cell.positions, displacements.positions_angstrom, atol=1e-6)
            or not np.allclose(cell.cell, displacements.lattice_angstrom, atol=1e-6)):
        raise ValueError("the phonopy supercell of the structure is not the one of the displacements "
                         "(another structure, or another phonopy version)")
    forces = np.asarray(forces_ev_per_angstrom, dtype=float)
    expected = (len(displacements), len(cell), 3)
    if forces.shape != expected:
        raise ValueError(f"forces_ev_per_angstrom must have shape {expected}, got {forces.shape}")
    if reference_forces_ev_per_angstrom is not None:
        reference = np.asarray(reference_forces_ev_per_angstrom, dtype=float)
        if reference.shape != expected[1:]:
            raise ValueError(f"reference_forces_ev_per_angstrom must have shape {expected[1:]}, got {reference.shape}")
        forces = forces - reference[None]
    phonon.dataset = {
        "natom": len(cell),
        "first_atoms": [{"number": atom, "displacement": np.array(vector), "forces": forces[index]}
                        for index, (atom, vector) in enumerate(displacements.displacements)],
    }
    phonon.produce_force_constants()
    phonon.symmetrize_force_constants()

    def frequencies(points):
        phonon.run_qpoints(points)
        return np.asarray(phonon.qpoints.frequencies, dtype=float) * THZ_TO_CM1

    def acoustic(points, values):
        """``(q index, band)`` of the 3 acoustic modes at Gamma, if Gamma is among ``points``."""

        gamma = np.flatnonzero(np.all(np.abs(points - np.round(points)) < 1e-12, axis=1))
        if not len(gamma):
            return set()
        return {(int(gamma[0]), int(band)) for band in np.argsort(np.abs(values[gamma[0]]))[:3]}

    points = _grid(mesh, shifted_mesh)
    values = frequencies(points)
    zero = acoustic(points, values)
    kept = [float(values[q, band]) for q in range(len(points)) for band in range(values.shape[1]) if (q, band) not in zero]

    commensurate = _grid(displacements.supercell, False)
    lowest_per_q = {}
    for grid, grid_values in ((points, values), (commensurate, frequencies(commensurate))):
        skip = acoustic(grid, grid_values)
        for q in range(len(grid)):
            row = [grid_values[q, band] for band in range(grid_values.shape[1]) if (q, band) not in skip]
            if not row:  # Gamma of a one-atom cell: acoustic modes only
                continue
            key = tuple(round(float(x) % 1.0, 9) % 1.0 for x in grid[q])
            lowest_per_q[key] = min(float(min(row)), lowest_per_q.get(key, math.inf))
    imaginary = tuple(sorted(((q, value) for q, value in lowest_per_q.items() if value < 0), key=lambda item: item[1]))
    stability = PhononStability(len(lowest_per_q), imaginary, min(lowest_per_q.values()))

    composition: dict[str, int] = {}
    for site in structure:
        composition[site.specie.symbol] = composition.get(site.specie.symbol, 0) + 1
    vibrations = HarmonicVibrations(
        frequencies_cm1=tuple(kept), composition=composition, supercell_size=len(points),
        imaginary_modes=imaginary_modes, low_frequency_cutoff_cm1=low_frequency_cutoff_cm1,
    )
    return PhononModes(vibrations, stability, displacements.supercell, mesh, bool(shifted_mesh), len(zero),
                       len(displacements), displacements.displacement_angstrom, float(symprec))

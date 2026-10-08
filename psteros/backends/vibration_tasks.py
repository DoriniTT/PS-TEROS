"""Calcfunctions of :func:`psteros.build_vibrations_workgraph`.

This module imports AiiDA at import time, so psteros loads it only while
building a graph; the daemon imports it when it runs one, which is why
psteros must be installed in the environment of the AiiDA daemon.
"""

from __future__ import annotations

from itertools import product

from aiida import orm
from aiida_workgraph import task


def _copy_structure(structure: orm.StructureData, cell, positions) -> orm.StructureData:
    """A structure with the kinds of ``structure`` (kind names kept), a new cell and new site positions."""

    copy = orm.StructureData(cell=cell, pbc=structure.pbc)
    for kind in structure.kinds:
        copy.append_kind(kind)
    for kind_name, position in positions:
        copy.append_site(orm.Site(kind_name=kind_name, position=position))
    return copy


@task.calcfunction
def repeat_structure(structure: orm.StructureData, size: orm.List) -> orm.StructureData:
    """Supercell ``(na, nb, nc)`` of ``structure``, keeping the kind names."""

    na, nb, nc = size.get_list()
    a, b, c = (list(vector) for vector in structure.cell)
    positions = []
    for i, j, k in product(range(na), range(nb), range(nc)):
        shift = [i * a[x] + j * b[x] + k * c[x] for x in range(3)]
        positions.extend(
            (site.kind_name, [site.position[x] + shift[x] for x in range(3)]) for site in structure.sites
        )
    cell = [[na * value for value in a], [nb * value for value in b], [nc * value for value in c]]
    return _copy_structure(structure, cell, positions)


@task.calcfunction
def displace_site(structure: orm.StructureData, site: orm.Int, vector: orm.List) -> orm.StructureData:
    """``structure`` with site ``site`` moved by ``vector`` (Cartesian, A)."""

    shift = vector.get_list()
    positions = [
        (current.kind_name, [current.position[x] + (shift[x] if index == site.value else 0.0) for x in range(3)])
        for index, current in enumerate(structure.sites)
    ]
    return _copy_structure(structure, structure.cell, positions)


def _forces(folder: orm.FolderData, backend: str, symbols: list[str]) -> list[list[float]]:
    from psteros.vibrations_workflow import qe_output_forces, vasprun_forces

    if backend == "vasp":
        found, forces = vasprun_forces(folder.get_object_content("vasprun.xml"))
        if found != symbols:
            raise ValueError("the atoms of vasprun.xml are not in the order of the structure")
        return forces
    return qe_output_forces(folder.get_object_content("aiida.out"), len(symbols))


def _is_linear(structure: orm.StructureData) -> bool:
    import numpy as np

    if len(structure.sites) <= 2:
        return True
    cell = np.array(structure.cell)
    positions = np.array([site.position for site in structure.sites])
    fractional = (positions - positions[0]) @ np.linalg.inv(cell)
    centered = (fractional - np.round(fractional)) @ cell  # minimum image around the first atom
    centered -= centered.mean(axis=0)
    singular = np.linalg.svd(centered, compute_uv=False)
    return bool(singular[1] < 1e-3 * singular[0])


@task.calcfunction
def harmonic_modes(structure: orm.StructureData, settings: orm.Dict, **retrieved) -> orm.Dict:
    """Harmonic modes from the forces of the displaced static calculations.

    ``retrieved`` holds the retrieved folder of every displacement under the
    key of :func:`psteros.vibrations_workflow.displacement_key`.
    """

    from psteros.vibrations import harmonic_vibrations_from_forces
    from psteros.vibrations_workflow import displacement_key

    options = settings.get_dict()
    sites = list(options["displaced_sites"])
    backend = options["backend"]
    supercell_size = int(options["supercell_size"])
    symbols = [structure.get_kind(site.kind_name).symbol for site in structure.sites]
    masses = [structure.get_kind(site.kind_name).mass for site in structure.sites]

    forces = {"plus": [], "minus": []}
    for site in sites:
        for axis in range(3):
            for sign in forces:
                forces[sign].append(_forces(retrieved[displacement_key(site, axis, sign)], backend, symbols))

    if options.get("molecule"):
        zero_modes = 5 if _is_linear(structure) else 6
    elif len(sites) == len(symbols):
        zero_modes = 3
    else:
        zero_modes = 0

    composition: dict[str, int] = {}
    frozen: dict[str, int] = {}
    displaced = set(sites)
    for index, symbol in enumerate(symbols):
        composition[symbol] = composition.get(symbol, 0) + 1
        if index not in displaced:
            frozen[symbol] = frozen.get(symbol, 0) + 1
    if any(count % supercell_size for count in composition.values()):
        raise ValueError("the supercell composition is not a multiple of the supercell size")

    vibrations = harmonic_vibrations_from_forces(
        masses_amu=masses,
        displaced_sites=sites,
        displacement_angstrom=float(options["displacement_angstrom"]),
        forces_plus=forces["plus"],
        forces_minus=forces["minus"],
        composition={symbol: count // supercell_size for symbol, count in composition.items()},
        frozen_composition=frozen,
        zero_modes=zero_modes,
        supercell_size=supercell_size,
        imaginary_modes="drop",  # stored as negative numbers; read_vibrations applies the policy
    )
    return orm.Dict({
        **vibrations.to_dict(),
        "displaced_sites": sites,
        "displacement_angstrom": float(options["displacement_angstrom"]),
        "zero_modes": zero_modes,
    })

"""AiiDA tasks of the VASP reference blocks (aiida-vasp ``vasp.v2.vasp``).

This module imports AiiDA at import time, so psteros loads it only while
building a graph.  The daemon imports the calcfunctions by module path,
which is why psteros must be installed in the daemon's Python environment.
"""

from __future__ import annotations

import re
from typing import Any, Mapping

from aiida import orm
from aiida.engine import calcfunction

from psteros.thermochemistry import parse_vasp_frequencies_cm1

# Files kept in the retrieved folder of every block (aiida-vasp otherwise
# retrieves OUTCAR and vasprun.xml only temporarily, for parsing).
RETRIEVE = ["OUTCAR", "vasprun.xml", "CONTCAR", "OSZICAR"]


@calcfunction
def vasp_energy(misc: orm.Dict, retrieved: orm.FolderData) -> orm.Float:
    """Total energy (eV, sigma -> 0) from the VASP misc output, else the last TOTEN in OUTCAR."""

    energies = misc.get_dict()
    if isinstance(energies.get("total_energies"), dict):
        energies = energies["total_energies"]
    for key in ("energy_extrapolated", "energy_no_entropy", "energy"):
        if key in energies:
            return orm.Float(energies[key])
    if "OUTCAR" in retrieved.base.repository.list_object_names():
        matches = re.findall(
            r"free\s+energy\s+TOTEN\s+=\s+([-\d.]+)", retrieved.base.repository.get_object_content("OUTCAR")
        )
        if matches:
            return orm.Float(float(matches[-1]))
    raise ValueError(f"no total energy in the VASP outputs; misc keys: {sorted(energies)}")


@calcfunction
def vasp_frequencies(retrieved: orm.FolderData) -> orm.List:
    """Frequencies (cm^-1, imaginary negative) of the dynamical matrix in the retrieved OUTCAR."""

    return orm.List(list(parse_vasp_frequencies_cm1(retrieved.base.repository.get_object_content("OUTCAR"))))


@calcfunction
def make_supercell(structure: orm.StructureData, size: orm.List) -> orm.StructureData:
    """Repeat a structure ``size = [nx, ny, nz]`` times along its lattice vectors.

    The atoms are grouped by element, in the order the elements first appear in
    ``structure``.  ``ase`` repeats cell by cell (``Sn Sn O O O O Sn Sn ...``)
    and VASP reads each run of equal elements of the POSCAR as an ion type of
    its own, so an ungrouped supercell has no symmetry for VASP: it found only
    the identity, planned 216 instead of a handful of displacements for a
    72-atom SnO2 cell and asked for days of walltime.
    """

    atoms = structure.get_ase().repeat(tuple(size.get_list()))
    symbols = atoms.get_chemical_symbols()
    rank: dict[str, int] = {}
    for symbol in symbols:
        rank.setdefault(symbol, len(rank))
    # sorted() is stable: the order inside an element stays cell by cell.
    atoms = atoms[sorted(range(len(atoms)), key=lambda index: rank[symbols[index]])]
    return orm.StructureData(ase=atoms)


def add_vasp_block_task(
    workgraph: Any,
    *,
    name: str,
    structure: Any,
    code: Any,
    incar: Mapping[str, Any],
    namespaces: Mapping[str, Any],
    kpoints_spacing: float,
    potential_family: str,
    potential_mapping: Mapping[str, str],
    options: Mapping[str, Any],
    settings: Mapping[str, Any],
    clean_workdir: bool,
) -> Any:
    """Add one ``vasp.v2.vasp`` task with the INCAR in aiida-vasp's ``incar`` namespace."""

    from aiida.plugins import WorkflowFactory
    from aiida_workgraph import task

    parameters = {**dict(namespaces), "incar": dict(incar)}
    return workgraph.add_task(
        task(WorkflowFactory("vasp.v2.vasp")),
        name=name,
        structure=structure,
        code=code,
        parameters=orm.Dict(parameters),
        kpoints_spacing=orm.Float(kpoints_spacing),
        potential_family=orm.Str(potential_family),
        potential_mapping=orm.Dict(dict(potential_mapping)),
        options=orm.Dict(dict(options)),
        settings=orm.Dict(dict(settings)),
        clean_workdir=orm.Bool(clean_workdir),
    )

"""Restore the kind names that aiida-vasp drops from a relaxed structure.

The relaxed ``StructureData`` of a ``VaspWorkChain`` has kinds named after the
elements, so a pseudo-hydrogen (kind ``H0p75``) comes back as ``H`` and the
static calculation that follows cannot find its POTCAR (``KeyError: 'H'``).
:func:`restore_vasp_kinds` gives every relaxed atom the kind of the starting
structure atom it came from.
"""

from __future__ import annotations

import numpy as np
from aiida import orm
from aiida.engine import calcfunction

from psteros.backends.vasp import assign_kinds


def _fractional(structure: orm.StructureData) -> np.ndarray:
    positions = np.array([site.position for site in structure.sites], dtype=float)
    return np.linalg.solve(np.array(structure.cell, dtype=float).T, positions.T).T


@calcfunction
def restore_vasp_kinds(relaxed: orm.StructureData, template: orm.StructureData) -> orm.StructureData:
    """The relaxed structure with the kinds (name, element, mass) of ``template``."""

    symbols = [relaxed.get_kind(site.kind_name).symbol for site in relaxed.sites]
    template_symbols = [template.get_kind(site.kind_name).symbol for site in template.sites]
    order = assign_kinds(template_symbols, _fractional(template), symbols, _fractional(relaxed), relaxed.cell)
    result = orm.StructureData(cell=relaxed.cell)
    result.pbc = relaxed.pbc
    for site, source in zip(relaxed.sites, order):
        kind = template.get_kind(template.sites[source].kind_name)
        result.append_atom(
            position=site.position, symbols=kind.symbols, weights=kind.weights, name=kind.name, mass=kind.mass
        )
    return result

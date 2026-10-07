#!/usr/bin/env python
"""
Preview charge-neutral slab terminations before running DFT.

Pure Python (pymatgen only, no AiiDA). For three surfaces this script prints
the table of terminations, writes the slabs as POSCAR files and saves side
views:

1. Ag3PO4(110): PO4 kept whole; one Ag per face removed from a charged slab.
2. GaAs(100): needs a 2x1 cell; half of the outer Ga or As layer remains.
3. ZnO(0001): a polar direction, reported instead of returning charged slabs.

Usage:
    python preview_terminations.py
"""

import os

from pymatgen.core import Lattice, Structure

from psteros.core.terminations import (
    NoChargeNeutralTerminationError,
    find_charge_neutral_terminations,
)

HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT = os.path.join(HERE, 'output')


def show(name, terminations):
    print(f"\n{'=' * 78}\n{name}\n{'=' * 78}")
    print(terminations)
    folder = os.path.join(OUTPUT, name)
    terminations.write(folder)
    terminations.plot(os.path.join(OUTPUT, f"{name}.png"))
    print(f"\nSlabs in {folder}/, side views in {OUTPUT}/{name}.png")


def main():
    # 1. Ag3PO4(110). PO4 groups are kept whole with unit_bonds.
    ag3po4 = Structure.from_file(os.path.join(HERE, '..', 'vasp', 'structures', 'ag3po4.cif'))
    show('Ag3PO4_110', find_charge_neutral_terminations(
        ag3po4, (1, 1, 0), min_slab_thickness=10.0, min_vacuum_thickness=15.0,
        oxidation_states={'Ag': 1, 'P': 5, 'O': -2},
        unit_bonds={('P', 'O'): 1.9},
    ))

    # 2. GaAs(100). Oxidation states are guessed (Ga +3, As -3); the 1x1 cell
    #    cannot be neutralised, so a 2x1 cell is chosen automatically.
    gaas = Structure.from_spacegroup('F-43m', Lattice.cubic(5.653), ['Ga', 'As'],
                                     [[0, 0, 0], [0.25, 0.25, 0.25]])
    show('GaAs_100', find_charge_neutral_terminations(gaas, (1, 0, 0), 10.0))

    # 3. ZnO(0001) is polar: no slab of it has two equivalent faces.
    zno = Structure.from_spacegroup('P6_3mc', Lattice.hexagonal(3.25, 5.21), ['Zn', 'O'],
                                    [[1 / 3, 2 / 3, 0], [1 / 3, 2 / 3, 0.382]])
    print(f"\n{'=' * 78}\nZnO_0001\n{'=' * 78}")
    try:
        find_charge_neutral_terminations(zno, (0, 0, 1), 10.0)
    except NoChargeNeutralTerminationError as error:
        print(error)


if __name__ == '__main__':
    main()

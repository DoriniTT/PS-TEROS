#!/usr/bin/env python
"""
Absolute surface energies of GaAs (111) and (-1-1-1) with VASP.

1. Without arguments: build the slabs and every reference calculation, print
   what will run, and write the slabs and their side views to ./output.
2. With --submit: also submit the AiiDA WorkGraph (edit CODE, COMPUTER and
   QUEUE first; the POTCAR family must contain H.75 and H1.25).
3. With --analyse PK: read the finished graph and write the phase diagram.

The bulk and the elemental references must be relaxed with the same settings
beforehand; the bulk lattice is kept for the slabs.

Method: Zhang et al., Sci. Rep. 6, 20055 (2016); arXiv:1510.08961.
"""

import argparse
import os

from pymatgen.core import Lattice, Structure

import psteros

CODE = "vasp@cluster"
COMPUTER = "cluster"
QUEUE = "normal"
POTCARS = {"Ga": "Ga_d", "As": "As"}
HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT = os.path.join(HERE, "output")


def structures():
    """Replace these by your relaxed structures."""
    bulk = Structure.from_spacegroup("F-43m", Lattice.cubic(5.75), ["Ga", "As"], [[0, 0, 0], [0.25, 0.25, 0.25]])
    gallium = Structure.from_spacegroup("Cmce", Lattice.orthorhombic(4.52, 7.66, 4.53), ["Ga"], [[0, 0.1549, 0.081]])
    arsenic = Structure.from_spacegroup("R-3m", Lattice.hexagonal(3.76, 10.55), ["As"], [[0, 0, 0.227]])
    return bulk, gallium, arsenic


def study():
    bulk, gallium, arsenic = structures()
    return psteros.PolarSurfaceStudy(
        bulk,
        faces=[(1, 1, 1), (-1, -1, -1)],
        references={"Ga": gallium, "As": arsenic},
        bilayers=9,
        pseudo_hydrogen_method="molecules",
        nonpolar_check=(1, 1, 0),
    )


def config(calculations):
    return psteros.SurfaceWorkflowConfig(
        backend="vasp",
        name="gaas_polar",
        calculation=psteros.VaspCalculationConfig(
            code_label=CODE,
            incar={"ENCUT": 400, "PREC": "Accurate", "EDIFF": 1e-6, "ISMEAR": 0, "SIGMA": 0.05,
                   "IBRION": 2, "NSW": 200, "EDIFFG": -0.005, "LREAL": False},
            potential_family="PBE",
            potential_mapping=calculations.potential_mapping(POTCARS),
            kpoints_spacing=0.2,
        ),
        execution=psteros.ExecutionPolicy(computer=COMPUTER, queue=QUEUE, max_concurrent_jobs=1),
        role_overrides=calculations.vasp_overrides(),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--submit", action="store_true", help="submit the WorkGraph")
    parser.add_argument("--analyse", type=int, metavar="PK", help="analyse a finished WorkGraph")
    args = parser.parse_args()

    calculations = study()
    os.makedirs(OUTPUT, exist_ok=True)
    for face, slabs in calculations.face_sets.items():
        print(f"\n{slabs}")
        slabs.write(os.path.join(OUTPUT, face))
        slabs.plot(os.path.join(OUTPUT, f"{face}.png"))
    print("\nCalculations:")
    for label, role in calculations.roles.items():
        print(f"  {label:34s} {role:20s} {len(calculations.structures[label])} atoms")
    print(f"\nPOTCARs: {calculations.potential_mapping(POTCARS)}")

    if args.submit or args.analyse:
        from aiida import load_profile, orm

        load_profile()
        if args.submit:
            graph = psteros.build_surface_workgraph(calculations.structures, config(calculations), submit=True)
            print(f"\nSubmitted WorkGraph {graph.pk}; analyse it with --analyse {graph.pk}")
        else:
            energies, relaxed = psteros.read_vasp_results(orm.load_node(args.analyse), calculations.structures)
            result = calculations.analyse(energies, relaxed)
            print(f"\n{result}")
            for face, report in result.bottom_checks.items():
                print(f"\n{face}\n{report}")
            result.diagram.plot(os.path.join(OUTPUT, "gaas_polar_phase_diagram.png"), title="GaAs polar faces")
            result.diagram.to_csv(os.path.join(OUTPUT, "gaas_polar_phase_diagram.csv"))
            print(f"\nPhase diagram written to {OUTPUT}/")


if __name__ == "__main__":
    main()

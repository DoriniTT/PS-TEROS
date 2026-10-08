#!/usr/bin/env python
"""
Surface phase diagram of the non-polar ZnO faces with VASP.

1. Without arguments: cut the charge-neutral terminations of ZnO (10-10) and
   (11-20), print every calculation of the set and write the slabs and their
   side views to ./output/zno_study.
2. With --submit: also submit the relaxation -> static WorkGraph (edit CODE,
   POTCAR_FAMILY, COMPUTER and QUEUE first).
3. With --analyse PK: read the finished graph and write the phase diagram.

The bulk must be relaxed beforehand with the same settings; the slabs are cut
at its cell.
"""

import argparse
import os

from pymatgen.core import Lattice, Structure

import psteros

CODE = "vasp@my-cluster"
POTCAR_FAMILY = "PBE"
POTCARS = {"Zn": "Zn", "O": "O"}
COMPUTER = "my-cluster"
QUEUE = "my-queue"
HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT = os.path.join(HERE, "output", "zno_study")

ELECTRONIC = {"ENCUT": 520, "PREC": "Accurate", "EDIFF": 1.0e-6, "ISMEAR": 0, "SIGMA": 0.05, "LREAL": False}
RELAX = {"IBRION": 2, "NSW": 200, "EDIFFG": -0.01}
STATIC = {"IBRION": -1, "NSW": 0}


def structures():
    """Replace these by your relaxed structures."""
    bulk = Structure.from_spacegroup("P6_3mc", Lattice.hexagonal(3.25, 5.21), ["Zn", "O"],
                                     [[1 / 3, 2 / 3, 0], [1 / 3, 2 / 3, 0.382]])
    zinc = Structure.from_spacegroup("P6_3/mmc", Lattice.hexagonal(2.66, 4.95), ["Zn"], [[1 / 3, 2 / 3, 0.25]])
    return bulk, zinc


def study():
    bulk, zinc = structures()
    return psteros.ChargeNeutralSurfaceStudy(
        bulk,
        miller_indices=[(1, 0, 0), (1, 1, 0)],          # (10-10) and (11-20)
        references={"Zn": zinc, "O": psteros.triplet_o2_cell(cell_length=12.0)},
        min_slab_thickness=12.0,
        oxidation_states={"Zn": 2, "O": -2},
    )


def recipes(calculations):
    execution = psteros.ExecutionPolicy(computer=COMPUTER, queue=QUEUE, max_concurrent_jobs=1)

    def vasp(ionic):
        return psteros.VaspCalculationConfig(
            code_label=CODE, incar={**ELECTRONIC, **ionic}, potential_family=POTCAR_FAMILY,
            potential_mapping=calculations.potential_mapping(POTCARS), kpoints_spacing=0.25, max_iterations=3,
        )

    relax = psteros.SurfaceWorkflowConfig(backend="vasp", calculation=vasp(RELAX), execution=execution,
                                          name="zno_surfaces", role_overrides=calculations.vasp_overrides("relax"))
    static = psteros.SurfaceWorkflowConfig(backend="vasp", calculation=vasp(STATIC), execution=execution,
                                           name="zno_surfaces_static",
                                           role_overrides=calculations.vasp_overrides("static"))
    return relax, static


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--submit", action="store_true", help="submit the WorkGraph")
    parser.add_argument("--analyse", type=int, metavar="PK", help="analyse a finished WorkGraph")
    args = parser.parse_args()

    calculations = study()
    print(calculations)
    os.makedirs(OUTPUT, exist_ok=True)
    for prefix, terminations in calculations.termination_sets.items():
        terminations.write(os.path.join(OUTPUT, prefix))
        terminations.plot(os.path.join(OUTPUT, f"{prefix}.png"))
    print("\nCalculations:")
    for label, role in calculations.roles.items():
        print(f"  {label:20s} {role:16s} {len(calculations.structures[label])} atoms")

    if args.submit or args.analyse:
        from aiida import load_profile, orm

        load_profile()
        if args.submit:
            relax, static = recipes(calculations)
            graph = psteros.build_relax_static_workgraph(calculations.structures, relax, static, submit=True)
            print(f"\nSubmitted WorkGraph {graph.pk}; analyse it with --analyse {graph.pk}")
        else:
            energies, relaxed = psteros.read_vasp_results(orm.load_node(args.analyse), calculations.structures)
            result = calculations.analyse(energies, relaxed)
            print(f"\n{result}")
            result.diagram.plot(os.path.join(OUTPUT, "zno_phase_diagram.png"), title="ZnO non-polar faces")
            result.diagram.to_csv(os.path.join(OUTPUT, "zno_phase_diagram.csv"))
            print(f"\nPhase diagram written to {OUTPUT}/")


if __name__ == "__main__":
    main()

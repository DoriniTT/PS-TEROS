"""SnO2 reference systems with VASP: energies and vibrational free energies.

Runs relax -> static -> vibrations for triplet O2, rutile SnO2 and alpha-Sn
(see docs/source/reference-thermochemistry.rst).  The settings are a smoke
test; converge ENCUT, k-points and the supercells before using the numbers.

Usage
-----
    python references.py --profile P --code VASP@cluster --computer cluster --queue debug --submit
    python references.py --profile P --show <PK of the graph> --temperature 600 --pressure 0.2
"""

from __future__ import annotations

import argparse

import psteros

INCAR = {
    "encut": 520,
    "prec": "Accurate",
    "ediff": 1e-7,
    "ismear": 0,
    "sigma": 0.05,
    "lreal": False,
    "ibrion": 2,
    "nsw": 100,
    "ediffg": -0.005,
    "lwave": False,
    "lcharg": False,
}
CELL_RELAX = psteros.CalculationOverride(parameters={"INCAR": {"isif": 3}})
TRIPLET_O2 = psteros.CalculationOverride(
    parameters={"INCAR": {"ispin": 2, "nupdown": 2}},
    kpoints_distance=5.0,  # Gamma only in the 12 A box
)


def build(args):
    recipe = psteros.SurfaceWorkflowConfig(
        backend="vasp",
        calculation=psteros.VaspCalculationConfig(
            code_label=args.code,
            incar=INCAR,
            potential_family=args.potential_family,
            potential_mapping={"Sn": "Sn_d", "O": "O"},
            kpoints_spacing=0.03,  # aiida-vasp units of 2*pi/A: about 0.19 1/A
        ),
        execution=psteros.ExecutionPolicy(
            computer=args.computer,
            queue=args.queue,
            resources={"num_machines": 1, "num_mpiprocs_per_machine": args.ranks},
            max_wallclock_seconds=args.walltime,
        ),
        name="sno2",
    )
    references = {
        "o2": psteros.ReferenceSystem(
            psteros.triplet_o2_cell(cell_length=12.0), "gas", symmetry_number=2, spin=1.0, override=TRIPLET_O2
        ),
        "sno2": psteros.ReferenceSystem(
            psteros.rutile_sno2_bulk(), "solid", supercell=(2, 2, 3), block_overrides={"relax": CELL_RELAX}
        ),
        "sn": psteros.ReferenceSystem(
            psteros.alpha_sn_bulk(), "solid", supercell=(2, 2, 2), block_overrides={"relax": CELL_RELAX}
        ),
    }
    graph = psteros.build_vasp_reference_workgraph(references, recipe, submit=args.submit)
    print(f"submitted WorkGraph {graph.pk}" if args.submit else "built (not submitted); pass --submit")


def show(args):
    for label, blocks in psteros.reference_results(args.show).items():
        for block, result in blocks.items():
            value = result.get("energy") if result["kind"] != "vibrations" else result.get("frequencies")
            print(f"{label:>6} {block:<11} {result['state']:<14} {value}")
    systems = psteros.reference_thermochemistry(args.show)
    print(f"\nT = {args.temperature} K, p(O2) = {args.pressure} bar")
    for label, terms in psteros.free_energies(systems, args.temperature, args.pressure).items():
        print(
            f"{label:>6}  E = {terms.electronic_energy_ev:12.5f}  ZPE = {terms.zero_point_energy_ev:8.4f}"
            f"  H(T)-H(0) = {terms.thermal_enthalpy_ev:8.4f}  -TS = {terms.entropy_term_ev:8.4f}"
            f"  G = {terms.free_energy_ev:12.5f} eV"
        )
    delta_mu = psteros.delta_mu_oxygen_ev(systems["o2"], args.temperature, args.pressure)
    print(f"\nDelta mu_O(T, p) = {delta_mu:.4f} eV on the axis mu_O = E(O2)/2 + Delta mu_O")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--profile")
    parser.add_argument("--code")
    parser.add_argument("--potential-family", default="PBE")
    parser.add_argument("--computer")
    parser.add_argument("--queue")
    parser.add_argument("--ranks", type=int, default=32)
    parser.add_argument("--walltime", type=int, default=6 * 3600)
    parser.add_argument("--submit", action="store_true")
    parser.add_argument("--show", type=int, help="PK of a finished reference graph to analyse")
    parser.add_argument("--temperature", type=float, default=298.15)
    parser.add_argument("--pressure", type=float, default=1.0)
    args = parser.parse_args()

    from aiida import load_profile

    load_profile(args.profile)
    if args.show is not None:
        show(args)
    else:
        missing = [name for name in ("code", "computer", "queue") if getattr(args, name) is None]
        if missing:
            parser.error(f"building a graph needs --{', --'.join(missing)}")
        build(args)


if __name__ == "__main__":
    main()

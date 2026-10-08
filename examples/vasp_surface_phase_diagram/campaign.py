"""SnO2(110) surface phase diagram with VASP: the calculations.

A small, end-to-end campaign built only on the public psteros API, sized as a
smoke test: 3 triple-layer slabs, 10 A vacuum, 400 eV and coarse k-points.
Converge these before drawing physical conclusions.

Phase ``refs``   bulk rutile SnO2 and alpha-Sn (cell relaxation, ISIF=3) and
                 triplet O2 (Gamma only), each followed by a static calculation.
Phase ``slabs``  the o / sno / sn2o terminations of SnO2(110) built on the
                 *relaxed* bulk lattice, central triple layer fixed with
                 selective dynamics, relaxation -> static calculation.

Usage
-----
    python campaign.py refs  --profile P --code vasp@my-cluster --potential-family PBE --submit
    python campaign.py slabs --profile P --code vasp@my-cluster --potential-family PBE \\
        --refs-pk <PK of the finished refs graph> --submit

Then run ``phase_diagram.py`` on the two graph PKs.
"""

from __future__ import annotations

import argparse

import psteros

TERMINATIONS = ("o", "sno", "sn2o")
TRIPLE_LAYERS = 3
VACUUM = 10.0
POTCARS = {"Sn": "Sn_d", "O": "O"}

ELECTRONIC = {
    "ENCUT": 400,
    "PREC": "Accurate",
    "EDIFF": 1.0e-5,
    "ISMEAR": 0,
    "SIGMA": 0.05,
    "LREAL": False,
    "LWAVE": False,
    "LCHARG": False,
}
RELAX = {"IBRION": 2, "NSW": 100, "ISIF": 2, "EDIFFG": -0.02}
STATIC = {"IBRION": -1, "NSW": 0}

CELL_RELAXATION = psteros.CalculationOverride(parameters={"INCAR": {"ISIF": 3}})
TRIPLET_O2 = psteros.CalculationOverride(
    parameters={"INCAR": {"ISPIN": 2, "MAGMOM": [1.0, 1.0]}},
    kpoints_distance=10.0,  # Gamma only in the 12 A box
)


def execution(args) -> psteros.ExecutionPolicy:
    return psteros.ExecutionPolicy(
        computer=args.computer,
        queue=args.queue,
        resources={"num_machines": 1, "num_mpiprocs_per_machine": args.ranks},
        max_wallclock_seconds=args.walltime,
        with_mpi=True,
        max_concurrent_jobs=1,
    )


def recipes(args, name, relax_overrides=None, static_overrides=None):
    """Relaxation and static recipes sharing code, POTCARs, cutoff and k-points."""

    def vasp(ionic):
        return psteros.VaspCalculationConfig(
            code_label=args.code,
            incar={**ELECTRONIC, **ionic},
            potential_family=args.potential_family,
            potential_mapping=POTCARS,
            kpoints_spacing=0.3,
            # Restarts let the work chain continue a relaxation stopped by the walltime.
            max_iterations=3,
        )

    relax = psteros.SurfaceWorkflowConfig(
        backend="vasp", calculation=vasp(RELAX), execution=execution(args),
        name=name, role_overrides=relax_overrides or {},
    )
    static = psteros.SurfaceWorkflowConfig(
        backend="vasp", calculation=vasp(STATIC), execution=execution(args),
        name=f"{name}_static", role_overrides=static_overrides or {},
    )
    return relax, static


def build_refs(args):
    structures = {
        "sno2_bulk": psteros.rutile_sno2_bulk(),
        "alpha_sn": psteros.alpha_sn_bulk(),
        "o2": psteros.triplet_o2_cell(cell_length=12.0),
    }
    relax, static = recipes(
        args, "sno2_refs",
        relax_overrides={"sno2_bulk": CELL_RELAXATION, "alpha_sn": CELL_RELAXATION, "o2": TRIPLET_O2},
        static_overrides={"o2": TRIPLET_O2},
    )
    return psteros.build_relax_static_workgraph(structures, relax, static, submit=args.submit)


def relaxed_bulk_lattice(refs_pk: int) -> tuple[float, float]:
    """Return (a, c) of the relaxed rutile cell from a finished refs graph."""
    from aiida import orm

    relaxed = orm.load_node(refs_pk).outputs.sno2_bulk_relaxed_structure.get_pymatgen_structure()
    lengths = sorted(relaxed.lattice.abc)  # rutile: a = b > c
    return (lengths[1] + lengths[2]) / 2.0, lengths[0]


def build_slabs(args):
    a, c = relaxed_bulk_lattice(args.refs_pk) if args.refs_pk is not None else (args.a, args.c)
    print(f"building slabs on bulk lattice a={a:.4f} A, c={c:.4f} A")
    structures = {
        f"slab_{termination}": psteros.sno2_110_slab(
            termination=termination, triple_layers=TRIPLE_LAYERS, vacuum_angstrom=VACUUM, a=a, c=c
        )[0]
        for termination in TERMINATIONS
    }
    fixed = {
        label: psteros.CalculationOverride(fixed_sites=psteros.central_sites(slab, half_width=1.5))
        for label, slab in structures.items()
    }
    relax, static = recipes(args, "sno2_110_slabs", relax_overrides=fixed)
    return psteros.build_relax_static_workgraph(structures, relax, static, submit=args.submit)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("phase", choices=("refs", "slabs"))
    parser.add_argument("--profile", required=True)
    parser.add_argument("--code", required=True, help="AiiDA label of a vasp.vasp code")
    parser.add_argument("--potential-family", required=True, help="uploaded POTCAR family")
    parser.add_argument("--computer", help="name of your AiiDA computer (for the record)")
    parser.add_argument("--queue", help="queue or partition; the scheduler's default when left out")
    parser.add_argument("--ranks", type=int, default=4, help="MPI ranks per job")
    parser.add_argument("--walltime", type=int, default=3600, help="seconds per job")
    parser.add_argument("--refs-pk", type=int, help="slabs phase: PK of the finished refs graph")
    parser.add_argument("--a", type=float, default=4.737, help="bulk a when --refs-pk is not given")
    parser.add_argument("--c", type=float, default=3.186, help="bulk c when --refs-pk is not given")
    parser.add_argument("--submit", action="store_true", help="submit; otherwise only build and list the graph")
    args = parser.parse_args(argv)

    from aiida import load_profile

    load_profile(args.profile)
    graph = build_refs(args) if args.phase == "refs" else build_slabs(args)
    tasks = [task.name for task in graph.tasks if task.name not in ("graph_inputs", "graph_outputs", "graph_ctx")]
    print(f"graph {graph.name!r}: {len(tasks)} tasks, one active job at a time")
    print("  " + ", ".join(tasks))
    if args.submit:
        print(f"submitted: PK={graph.pk}")
    return graph


if __name__ == "__main__":
    main()

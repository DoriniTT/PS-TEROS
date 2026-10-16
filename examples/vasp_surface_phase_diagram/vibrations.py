"""SnO2(110) surface phase diagram with VASP: the vibrations (optional step).

Gamma-point finite-difference vibrations of the structures relaxed by
``campaign.py``, with the same POTCARs, cutoff and k-point spacing and a
tighter electronic convergence (EDIFF = 1e-7):

* the slabs: only the sites that were free in the relaxation are displaced;
  the fixed central triple layer (one Sn2O4 bulk cell) is counted as bulk;
* bulk rutile SnO2 and alpha-Sn: every site of a Gamma-only supercell;
* O2: the stretch only (its zero-point energy enters the O2 reference).

Each displaced site costs 6 static calculations. ``--parallel-jobs`` lets
several of them run at once (one by default).

    python vibrations.py --profile P --code vasp@my-cluster --potential-family PBE \\
        --refs-pk <REFS_PK> --slabs-pk <SLABS_PK> --submit

Then pass the PK of this graph to ``phase_diagram.py --vib-pk ... --temperature ...``.
"""

from __future__ import annotations

import argparse

import psteros
from campaign import (
    ELECTRONIC, POTCARS, STATIC, TERMINATIONS, TRIPLE_LAYERS, TRIPLET_O2, VACUUM, relaxed_bulk_lattice,
)

VIBRATIONS_EDIFF = 1.0e-7


def fixed_sites(refs_pk: int) -> dict[str, tuple[int, ...]]:
    """The sites campaign.py fixed in each slab, from the same slabs built on the same lattice."""

    a, c = relaxed_bulk_lattice(refs_pk)
    return {
        f"slab_{termination}": tuple(psteros.central_sites(
            psteros.sno2_110_slab(termination=termination, triple_layers=TRIPLE_LAYERS,
                                  vacuum_angstrom=VACUUM, a=a, c=c)[0],
            half_width=1.5,
        ))
        for termination in TERMINATIONS
    }


def static_recipe(args, fixed: dict[str, tuple[int, ...]]) -> psteros.SurfaceWorkflowConfig:
    """The static recipe of campaign.py with a tight EDIFF; fixed sites mark the frozen slab centres."""

    overrides = {label: psteros.CalculationOverride(fixed_sites=sites) for label, sites in fixed.items()}
    overrides["o2"] = TRIPLET_O2
    return psteros.SurfaceWorkflowConfig(
        backend="vasp",
        name="sno2_110_vibrations",
        calculation=psteros.VaspCalculationConfig(
            code_label=args.code,
            incar={**ELECTRONIC, **STATIC, "EDIFF": VIBRATIONS_EDIFF},
            potential_family=args.potential_family,
            potential_mapping=POTCARS,
            kpoints_spacing=0.3,
        ),
        execution=psteros.ExecutionPolicy(
            computer=args.computer,
            queue=args.queue,
            resources={"num_machines": 1, "num_mpiprocs_per_machine": args.ranks},
            max_wallclock_seconds=args.walltime,
            max_concurrent_jobs=args.parallel_jobs,
        ),
        role_overrides=overrides,
    )


def build_vibrations(args):
    from aiida import orm

    fixed = fixed_sites(args.refs_pk)
    _, references = psteros.read_vasp_results(orm.load_node(args.refs_pk), ("sno2_bulk", "alpha_sn", "o2"))
    _, slabs = psteros.read_vasp_results(orm.load_node(args.slabs_pk), tuple(fixed))
    return psteros.build_vibrations_workgraph(
        {**references, **slabs},
        static_recipe(args, fixed),
        psteros.VibrationsConfig(
            displacement_angstrom=0.01,
            supercells={"sno2_bulk": tuple(args.bulk_supercell)},
            molecules=("o2",),
        ),
        submit=args.submit,
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--code", required=True, help="AiiDA label of a vasp.vasp code")
    parser.add_argument("--potential-family", required=True, help="uploaded POTCAR family")
    parser.add_argument("--refs-pk", type=int, required=True, help="PK of the finished refs graph")
    parser.add_argument("--slabs-pk", type=int, required=True, help="PK of the finished slabs graph")
    parser.add_argument("--computer", help="name of your AiiDA computer (for the record)")
    parser.add_argument("--queue", help="queue or partition; the scheduler's default when left out")
    parser.add_argument("--ranks", type=int, default=4, help="MPI ranks per job")
    parser.add_argument("--walltime", type=int, default=3600, help="seconds per job")
    parser.add_argument("--bulk-supercell", type=int, nargs=3, default=(2, 2, 2),
                        help="Gamma-only supercell of the rutile cell (alpha-Sn uses its 8-atom cell)")
    parser.add_argument("--parallel-jobs", type=int, default=1, help="static calculations running at once")
    parser.add_argument("--submit", action="store_true", help="submit; otherwise only build and count the graph")
    args = parser.parse_args(argv)

    from aiida import load_profile

    load_profile(args.profile)
    graph = build_vibrations(args)
    jobs = sum(task.name.endswith("_vasp") for task in graph.tasks)
    print(f"graph {graph.name!r}: {jobs} static calculations, {args.parallel_jobs} at a time")
    if args.submit:
        print(f"submitted: PK={graph.pk}")
    return graph


if __name__ == "__main__":
    main()

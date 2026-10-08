"""SnO2(110) surface phase diagram with VASP and vibrational references: the calculations.

An end-to-end campaign built only on the public psteros API (the VASP version
of ``examples/qe_surface_phase_diagram``).  It runs on Lovelace (queue
``par128``, one 128-core node, one job at a time).

Phase ``o2``     triplet O2 alone (relax -> static -> vibrations): a cheap
                 smoke test of every code path.
Phase ``refs``   O2, bulk rutile SnO2 and alpha-Sn, each relax -> static ->
                 vibrations (``psteros.build_vasp_reference_workgraph``).
Phase ``slabs``  the o / sno / sn2o terminations of SnO2(110), built on the
                 *relaxed* SnO2 lattice of the refs graph, one relaxation each
                 (``psteros.build_surface_workgraph``).

Every graph is built and printed (INCAR, k-points, scheduler options,
settings and structure of each VASP task) before anything is submitted.

Usage
-----
    python campaign.py o2                       # build and print only
    python campaign.py o2 --submit
    python campaign.py refs --submit
    python campaign.py slabs --refs-pk <PK of the finished refs graph> --submit

Then run ``analysis.py`` on the graph PKs.
"""

from __future__ import annotations

import argparse

import psteros

PROFILE = "psteros_sno2_vibrations"
CODE = "VASP-6.5.1@lovelace"
POTENTIAL_FAMILY = "PBE"
POTENTIAL_MAPPING = {"Sn": "Sn_d", "O": "O"}

TERMINATIONS = ("o", "sno", "sn2o")
TRIPLE_LAYERS = 3
VACUUM = 15.0

# One INCAR for the references and the slabs: same cutoff, precision,
# smearing and POTCARs everywhere, so that their energies can be combined.
INCAR = {
    "encut": 520,
    "prec": "Accurate",
    "ediff": 1e-7,
    "ismear": 0,
    "sigma": 0.05,
    "lreal": False,
    "lasph": True,
    "ibrion": 2,
    "nsw": 100,
    "ncore": 16,
    "lwave": False,
    "lcharg": False,
}
KPOINTS_SPACING = 0.03  # aiida-vasp units of 2*pi/A: about 0.19 1/A
REFERENCE_EDIFFG = -0.005  # eV/A, tight: the frequencies are computed at this minimum
SLAB_EDIFFG = -0.02

CELL_RELAX = psteros.CalculationOverride(parameters={"INCAR": {"isif": 3}})
TRIPLET_O2 = psteros.CalculationOverride(
    parameters={"INCAR": {"ispin": 2, "nupdown": 2}},
    kpoints_distance=5.0,  # Gamma only in the 12 A box
)

DEFAULT_WALLTIME_HOURS = {"o2": 4, "refs": 12, "slabs": 12}


def execution(args) -> psteros.ExecutionPolicy:
    return psteros.ExecutionPolicy(
        computer=args.computer,
        queue=args.queue,
        resources={"num_machines": 1, "num_cores_per_machine": 128, "num_mpiprocs_per_machine": 128},
        max_wallclock_seconds=int(args.walltime_hours * 3600),
        max_concurrent_jobs=1,
        # Lmod shell functions exported by the login environment break the qstat parser.
        extra_options={"import_sys_environment": False},
    )


def recipe(args, name: str, ediffg: float) -> psteros.SurfaceWorkflowConfig:
    return psteros.SurfaceWorkflowConfig(
        backend="vasp",
        calculation=psteros.VaspCalculationConfig(
            code_label=args.code,
            incar={**INCAR, "ediffg": ediffg},
            potential_family=args.potential_family,
            potential_mapping=POTENTIAL_MAPPING,
            kpoints_spacing=KPOINTS_SPACING,
        ),
        execution=execution(args),
        name=name,
    )


def reference_systems(labels) -> dict[str, psteros.ReferenceSystem]:
    systems = {
        "o2": psteros.ReferenceSystem(
            psteros.triplet_o2_cell(cell_length=12.0), "gas", override=TRIPLET_O2, symmetry_number=2, spin=1.0
        ),
        "sno2": psteros.ReferenceSystem(
            psteros.rutile_sno2_bulk(), "solid", block_overrides={"relax": CELL_RELAX}, supercell=(2, 2, 3)
        ),
        "sn": psteros.ReferenceSystem(
            psteros.alpha_sn_bulk(), "solid", block_overrides={"relax": CELL_RELAX}, supercell=(2, 2, 2)
        ),
    }
    return {label: systems[label] for label in labels}


def reference_blocks():
    return (psteros.Relax(), psteros.Static(), psteros.Vibrations(incar={"ediff": 1e-8}))


def build_references(args, labels, submit: bool):
    return psteros.build_vasp_reference_workgraph(
        reference_systems(labels),
        recipe(args, "sno2_refs" if len(labels) > 1 else "o2_alone", REFERENCE_EDIFFG),
        blocks=reference_blocks(),
        submit=submit,
    )


def relaxed_bulk_lattice(refs_pk: int) -> tuple[float, float]:
    """(a, c) of the relaxed rutile SnO2 cell of a finished refs graph."""

    relaxed = psteros.reference_results(refs_pk)["sno2"]["relax"]["structure"]
    if relaxed is None:
        raise SystemExit(f"graph {refs_pk} has no relaxed SnO2 structure yet")
    lengths = sorted(relaxed.get_pymatgen_structure().lattice.abc)  # rutile: c < a = b
    return (lengths[1] + lengths[2]) / 2.0, lengths[0]


def build_slabs(args, submit: bool):
    a, c = relaxed_bulk_lattice(args.refs_pk) if args.refs_pk is not None else (args.a, args.c)
    print(f"building slabs on the SnO2 lattice a={a:.4f} A, c={c:.4f} A")
    structures = {
        f"slab_{termination}": psteros.sno2_110_slab(
            termination=termination, triple_layers=TRIPLE_LAYERS, vacuum_angstrom=VACUUM, a=a, c=c
        )[0]
        for termination in TERMINATIONS
    }
    config = recipe(args, "sno2_110_slabs", SLAB_EDIFFG)
    config = psteros.SurfaceWorkflowConfig(
        backend=config.backend,
        calculation=psteros.VaspCalculationConfig(
            code_label=config.calculation.code_label,
            incar={**config.calculation.incar, "isif": 2},  # fixed cell, all atoms free
            potential_family=config.calculation.potential_family,
            potential_mapping=config.calculation.potential_mapping,
            kpoints_spacing=config.calculation.kpoints_spacing,
        ),
        execution=config.execution,
        name=config.name,
    )
    return psteros.build_surface_workgraph(structures, config, submit=submit)


def describe(graph) -> None:
    """Print what every VASP task of an unsubmitted graph will run."""

    names = [task.name for task in graph.tasks if task.name not in ("graph_inputs", "graph_outputs", "graph_ctx")]
    print(f"graph {graph.name!r}: {len(names)} tasks, max_number_jobs={graph.max_number_jobs}")
    for task in graph.tasks:
        if task.name.endswith("_supercell"):
            print(f"- {task.name}: repeat {task.inputs['size'].value.get_list()}")
        if not task.name.endswith("_vasp"):
            continue
        inputs = {key: task.inputs[key].value for key in task.inputs._get_keys() if key in _SHOWN}
        parameters = inputs["parameters"].get_dict()
        structure = inputs["structure"]
        print(f"- {task.name}")
        if structure is not None:
            print(f"    structure: {structure.get_formula()}, {len(structure.sites)} atoms, cell {_cell(structure)}")
        else:
            print("    structure: from the previous block (or its supercell)")
        print(f"    INCAR: {parameters['incar']}")
        extra = {key: value for key, value in parameters.items() if key != "incar"}
        if extra:
            print(f"    other parameters: {extra}")
        print(f"    kpoints_spacing: {inputs['kpoints_spacing'].value}")
        print(f"    potentials: {inputs['potential_family'].value} {inputs['potential_mapping'].get_dict()}")
        print(f"    options: {inputs['options'].get_dict()}")
        print(f"    settings: {inputs['settings'].get_dict()}")


_SHOWN = ("parameters", "structure", "kpoints_spacing", "potential_family", "potential_mapping", "options", "settings")


def _cell(structure) -> str:
    return " x ".join(f"{length:.2f}" for length in structure.get_pymatgen_structure().lattice.abc) + " A"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("phase", choices=("o2", "refs", "slabs"))
    parser.add_argument("--profile", default=PROFILE)
    parser.add_argument("--code", default=CODE)
    parser.add_argument("--potential-family", default=POTENTIAL_FAMILY)
    parser.add_argument("--computer", default="lovelace")
    parser.add_argument("--queue", default="par128")
    parser.add_argument("--walltime-hours", type=float, help="per job; default 4 (o2) or 12")
    parser.add_argument("--refs-pk", type=int, help="slabs phase: PK of the refs graph with the relaxed SnO2")
    parser.add_argument("--a", type=float, default=4.737, help="SnO2 a when --refs-pk is not given")
    parser.add_argument("--c", type=float, default=3.186, help="SnO2 c when --refs-pk is not given")
    parser.add_argument("--submit", action="store_true", help="submit; otherwise only build and print the graph")
    args = parser.parse_args(argv)
    if args.walltime_hours is None:
        args.walltime_hours = DEFAULT_WALLTIME_HOURS[args.phase]

    from aiida import load_profile

    load_profile(args.profile)
    labels = ("o2",) if args.phase == "o2" else ("o2", "sno2", "sn")

    def build(submit: bool):
        return build_slabs(args, submit) if args.phase == "slabs" else build_references(args, labels, submit)

    describe(build(submit=False))
    if args.submit:
        graph = build(submit=True)
        print(f"submitted: PK={graph.pk}")
        return graph
    print("not submitted (use --submit)")


if __name__ == "__main__":
    main()

"""SnO2(110) surface phase diagram with VASP and vibrational references, in one WorkGraph.

The campaign of ``examples/vasp_sno2_phase_diagram_vibrations`` as a single call of
``psteros.build_vasp_campaign_workgraph``, on Lovelace (queue ``par128``, one
128-core node, one job at a time):

references  O2 (triplet, gas), rutile SnO2 and alpha-Sn (solid): relax -> static -> vibrations
slabs       the o / sno / sn2o terminations of SnO2(110): relax -> static

The slabs are finished structures, cut from the SnO2 lattice given by --a and --c
(by default the relaxed lattice of the earlier two-graph run). The graph does not
rebuild them from the bulk it relaxes.

Every VASP task is built and printed (INCAR, k-points, scheduler options, settings
and structure) before anything is submitted. The graph has one PK, and
``analysis.py`` reads everything back from it.

Usage
-----
    python campaign.py                          # build and print only
    python campaign.py --submit                 # the whole campaign
    python campaign.py --smoke --submit         # triplet O2 alone: a cheap test of every code path
    python analysis.py --pk <PK> --out results/sno2_110
"""

from __future__ import annotations

import argparse
import math

import psteros

PROFILE = "psteros_sno2_vibrations"
CODE = "VASP-6.5.1@lovelace"
POTENTIAL_FAMILY = "PBE"
POTENTIAL_MAPPING = {"Sn": "Sn_d", "O": "O"}

TERMINATIONS = ("o", "sno", "sn2o")
TRIPLE_LAYERS = 3
VACUUM = 15.0
# Relaxed rutile SnO2 lattice (A) of the earlier two-graph run (LOG.md of the example):
# the slabs are cut from it, because the graph does not cut them from the bulk it relaxes.
RELAXED_A = 4.8301
RELAXED_C = 3.2434

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
TWO_PI = 2.0 * math.pi
# A^-1 with the 2*pi, as VASP's KSPACING: 0.1885 1/A (7x7x11 for the relaxed SnO2 cell, 6x6x6 for alpha-Sn)
KPOINTS_SPACING = 0.03 * TWO_PI
REFERENCE_EDIFFG = -0.005  # eV/A, tight: the frequencies are computed at this minimum
SLAB_EDIFFG = -0.02

CELL_RELAX = psteros.CalculationOverride(parameters={"INCAR": {"isif": 3}})
# SnO2 is a wide-gap insulator and the supercell is 9.7 A wide, so a 2x2x2 mesh is enough for the forces
# of the displaced cells; the fine mesh of the shared recipe would make the ~12 displacements of the
# 72-atom cell need days.  The relaxation and the static energy keep the fine mesh.
SNO2_VIBRATIONS = psteros.CalculationOverride(
    kpoints_distance=0.06 * TWO_PI, metadata={"max_wallclock_seconds": 24 * 3600}
)
TRIPLET_O2 = psteros.CalculationOverride(
    parameters={"INCAR": {"ispin": 2, "nupdown": 2}},
    kpoints_distance=5.0 * TWO_PI,  # Gamma only in the 12 A box
)
# The slabs: fixed cell (ISIF = 2), all atoms free, a looser EDIFFG than the references.
SLAB = psteros.CalculationOverride(parameters={"INCAR": {"isif": 2, "ediffg": SLAB_EDIFFG}})


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


def recipe(args) -> psteros.SurfaceWorkflowConfig:
    return psteros.SurfaceWorkflowConfig(
        backend="vasp",
        calculation=psteros.VaspCalculationConfig(
            code_label=args.code,
            incar={**INCAR, "ediffg": REFERENCE_EDIFFG},
            potential_family=args.potential_family,
            potential_mapping=POTENTIAL_MAPPING,
            kpoints_spacing=KPOINTS_SPACING,
        ),
        execution=execution(args),
        name="sno2_110",
    )


def reference_systems(smoke: bool) -> dict[str, psteros.ReferenceSystem]:
    systems = {
        "o2": psteros.ReferenceSystem(
            psteros.triplet_o2_cell(cell_length=12.0), "gas", override=TRIPLET_O2, symmetry_number=2, spin=1.0
        ),
    }
    if not smoke:
        systems["sno2"] = psteros.ReferenceSystem(
            psteros.rutile_sno2_bulk(), "solid",
            block_overrides={"relax": CELL_RELAX, "vibrations": SNO2_VIBRATIONS}, supercell=(2, 2, 3),
        )
        systems["sn"] = psteros.ReferenceSystem(
            psteros.alpha_sn_bulk(), "solid", block_overrides={"relax": CELL_RELAX}, supercell=(2, 2, 2)
        )
    return systems


def slab_systems(a: float, c: float) -> dict[str, psteros.SlabSystem]:
    """The three terminations, cut from the SnO2 lattice (a, c) in A."""

    print(f"slabs cut from the SnO2 lattice a={a:.4f} A, c={c:.4f} A")
    return {
        f"slab_{termination}": psteros.SlabSystem(
            psteros.sno2_110_slab(
                termination=termination, triple_layers=TRIPLE_LAYERS, vacuum_angstrom=VACUUM, a=a, c=c
            )[0],
            override=SLAB,
        )
        for termination in TERMINATIONS
    }


def build(args, submit: bool):
    """The one graph of the campaign, or of the smoke test (``--smoke``: triplet O2 alone)."""

    return psteros.build_vasp_campaign_workgraph(
        reference_systems(args.smoke),
        {} if args.smoke else slab_systems(args.a, args.c),
        recipe(args),
        reference_blocks=(psteros.Relax(), psteros.Static(), psteros.Vibrations(incar={"ediff": 1e-8})),
        slab_blocks=(psteros.Relax(), psteros.Static()),
        submit=submit,
    )


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
    parser.add_argument("--profile", default=PROFILE)
    parser.add_argument("--code", default=CODE)
    parser.add_argument("--potential-family", default=POTENTIAL_FAMILY)
    parser.add_argument("--computer", default="lovelace")
    parser.add_argument("--queue", default="par128")
    parser.add_argument("--walltime-hours", type=float, default=12.0, help="per job (default 12)")
    parser.add_argument("--a", type=float, default=RELAXED_A, help="SnO2 a (A) the slabs are cut from")
    parser.add_argument("--c", type=float, default=RELAXED_C, help="SnO2 c (A) the slabs are cut from")
    parser.add_argument("--smoke", action="store_true", help="triplet O2 alone, no slabs: a cheap test of every code path")
    parser.add_argument("--submit", action="store_true", help="submit; otherwise only build and print the graph")
    args = parser.parse_args(argv)

    from aiida import load_profile

    load_profile(args.profile)
    describe(build(args, submit=False))
    if args.submit:
        graph = build(args, submit=True)
        print(f"submitted: PK={graph.pk}")
        return graph
    print("not submitted (use --submit)")


if __name__ == "__main__":
    main()

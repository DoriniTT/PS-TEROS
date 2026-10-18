"""RuO2(110) with a fine-tuned MACE model as one WorkGraph (``psteros.build_mlip_surface_workgraph``).

CP2K statics of bulk RuO2, hcp Ru, triplet O2 and the six terminations, a replay fine-tune of MACE-MH-1 and the
relaxations on a GPU computer, one active-learning round and a CP2K static at each relaxed geometry; then the
surface energies from those statics (``"mlip_dft"``, the route validated in validation/ruo2_mace_lovelace).
The references must be relaxed with the same CP2K recipe first (``build_cp2k_relax_static_workgraph``); here the
lattice constants are typed in. The CP2K settings of ``ruo2_mlip.py`` are a smoke test: converge them first.

Usage
-----
    python ruo2_mlip_graph.py build --profile P --code cp2k@cpu-cluster --mace-code mace-python@gpu-cluster \\
        --model /path/on/gpu/mace-mh-1.model --replay replay-data-mh-1-omat-pbe.xyz      # add --submit
    python ruo2_mlip_graph.py read --profile P --pk <PK>
"""

from __future__ import annotations

import argparse

import psteros
from ruo2_mlip import DFT, O2, structures


def build(args) -> None:
    static = psteros.SurfaceWorkflowConfig(
        backend="cp2k", name="ruo2_mlip", role_overrides={"o2": O2},
        execution=psteros.ExecutionPolicy(resources={"num_machines": 1, "num_mpiprocs_per_machine": 128},
                                          max_wallclock_seconds=2 * 3600, max_concurrent_jobs=2),
        calculation=psteros.Cp2kCalculationConfig(
            code_label=args.code, parameters={"FORCE_EVAL": {"METHOD": "Quickstep", "DFT": DFT}},
            basis_sets={"Ru": "DZVP-MOLOPT-PBE-GTH-q16", "O": "DZVP-MOLOPT-PBE-GTH-q6"},
            potentials={"Ru": "GTH-PBE-q16", "O": "GTH-PBE-q6"},
            basis_set_file="BASIS_MOLOPT_UZH", potential_file="POTENTIAL_UZH", kpoints_spacing=0.20),
    )
    fine_tune = psteros.MaceFineTuneConfig(args.model, replay_file=args.replay)
    job = psteros.MlipJobConfig(args.mace_code, psteros.ExecutionPolicy(
        queue=args.queue, with_mpi=False, max_wallclock_seconds=12 * 3600,
        resources={"num_machines": 1, "num_mpiprocs_per_machine": 8, "ngpus_per_machine": 1},
        extra_options={"environment_variables": {"OMP_NUM_THREADS": "8", "PYTHONUNBUFFERED": "1"}}))
    items = structures()
    slabs = psteros.RUTILE_110_TERMINATIONS
    graph = psteros.build_mlip_surface_workgraph(
        items, static, fine_tune, job=job, metallic_bulk="bulk",
        fixed_sites={label: psteros.central_sites(items[label]) for label in slabs}, submit=args.submit)
    print(graph.name, len(graph.tasks), "tasks", f"submitted as {graph.pk}" if args.submit else "not submitted")


def read(args) -> None:
    from aiida import orm

    graph = orm.load_node(args.pk)
    items = structures()
    energies = psteros.read_mlip_energies(graph, list(items), energy_source="mlip_dft")
    relaxations = psteros.read_mlip_results(graph, list(psteros.RUTILE_110_TERMINATIONS))
    for label, relaxation in relaxations.items():
        print(f"{label:16s} E_DFT {energies[label]:.4f} eV  model {relaxation.energy_ev:.4f} eV  "
              f"{relaxation.steps} steps, largest force {relaxation.max_force_ev_per_angstrom:.4f} eV/A")
    print({label: energies[label] for label in ("bulk", "ru", "o2")})
    # These energies (one source for slabs and references) go into psteros.SlabTermination and
    # psteros.BinaryOxideReferences, see docs/source/phase-diagram.rst.


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="phase", required=True)
    first = sub.add_parser("build")
    for name in ("--profile", "--code", "--mace-code", "--model", "--replay"):
        first.add_argument(name, required=True)
    first.add_argument("--queue", default="gpu")
    first.add_argument("--submit", action="store_true")
    second = sub.add_parser("read")
    second.add_argument("--profile", required=True)
    second.add_argument("--pk", type=int, required=True)
    args = parser.parse_args()
    from aiida import load_profile

    load_profile(args.profile)
    {"build": build, "read": read}[args.phase](args)


if __name__ == "__main__":
    main()

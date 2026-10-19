"""RuO2(110) with a fine-tuned MACE model: CP2K labels, then the fine-tune and relaxation job folders.

Phase ``labels``  CP2K static calculations (ENERGY_FORCE) of bulk RuO2, hcp Ru, triplet O2 and the six
                  RuO2(110) terminations of psteros.RUTILE_110_TERMINATIONS, as built (no relaxation).
Phase ``jobs``    reads the finished ``labels`` graph and writes the job folders of a replay fine-tune of
                  MACE-MH-1 and of the relaxations, to copy to a GPU node and run there
                  (``python finetune.py``, then ``python relax.py``).

The CP2K settings are placeholders sized as a smoke test: converge the cutoff, basis, k-points and vacuum first.
psteros routes the MLIP method to metallic systems only; ``jobs`` checks that bulk RuO2 is a metal.

Usage
-----
    python ruo2_mlip.py labels --profile P --code cp2k@my-cluster            # builds the graph; add --submit
    python ruo2_mlip.py jobs --profile P --labels-pk <PK> --model /path/mace-mh-1.model \\
        --replay replay-data-mh-1-omat-pbe.xyz --out jobs
"""

from __future__ import annotations

import argparse

import psteros

A, C = 4.492, 3.107  # use the lattice of your relaxed bulk
BULK_LABELS = ("bulk", "ru", "o2")

DFT = {
    "QS": {"METHOD": "GPW", "EPS_DEFAULT": 1.0e-12},
    "MGRID": {"NGRIDS": 4, "CUTOFF": 500, "REL_CUTOFF": 60},
    "XC": {"XC_FUNCTIONAL": {"_": "PBE"}},
    "POISSON": {"PERIODIC": "XYZ", "PSOLVER": "PERIODIC"},
    "SCF": {"EPS_SCF": 1.0e-7, "MAX_SCF": 300, "ADDED_MOS": 40,
            "DIAGONALIZATION": {"ALGORITHM": "STANDARD"},
            "MIXING": {"METHOD": "BROYDEN_MIXING", "ALPHA": 0.2, "NBUFFER": 8},
            "SMEAR": {"METHOD": "FERMI_DIRAC", "ELECTRONIC_TEMPERATURE": 300}},
}
O2 = psteros.CalculationOverride(kpoints_distance=100.0,
                                 parameters={"FORCE_EVAL": {"DFT": {"UKS": True, "MULTIPLICITY": 3}}})


def structures():
    from pymatgen.core import Lattice, Structure

    result = {
        "bulk": psteros.rutile_bulk("Ru", a=A, c=C),
        "ru": Structure.from_spacegroup("P6_3/mmc", Lattice.hexagonal(2.706, 4.282), ["Ru"], [[1 / 3, 2 / 3, 0.25]]),
        "o2": psteros.triplet_o2_cell(),
    }
    for termination in psteros.RUTILE_110_TERMINATIONS:
        slab, _ = psteros.rutile_110_slab(metal="Ru", termination=termination, a=A, c=C, vacuum_angstrom=22.0)
        result[termination] = Structure(slab.lattice, slab.species, slab.frac_coords)
    return result


def labels(args) -> None:
    config = psteros.SurfaceWorkflowConfig(
        backend="cp2k", name="ruo2_labels", role_overrides={"o2": O2},
        execution=psteros.ExecutionPolicy(resources={"num_machines": 1, "num_mpiprocs_per_machine": 128},
                                          max_wallclock_seconds=2 * 3600, max_concurrent_jobs=3),
        calculation=psteros.Cp2kCalculationConfig(
            code_label=args.code, parameters={"FORCE_EVAL": {"METHOD": "Quickstep", "DFT": DFT}},
            basis_sets={"Ru": "DZVP-MOLOPT-PBE-GTH-q16", "O": "DZVP-MOLOPT-PBE-GTH-q6"},
            potentials={"Ru": "GTH-PBE-q16", "O": "GTH-PBE-q6"},
            basis_set_file="BASIS_MOLOPT_UZH", potential_file="POTENTIAL_UZH", kpoints_spacing=0.20),
    )
    graph = psteros.build_cp2k_surface_workgraph(structures(), config, submit=args.submit)
    print(graph.name, sorted(task.name for task in graph.tasks), "submitted" if args.submit else "not submitted")


def jobs(args) -> None:
    from aiida import orm

    node = orm.load_node(args.labels_pk)
    built = structures()
    energies, _ = psteros.read_cp2k_results(node, list(built))
    forces = psteros.read_cp2k_forces(node, list(built))
    dos = psteros.read_cp2k_dos(node, ["bulk"])["bulk"]
    psteros.require_metallic("bulk", partially_occupied=dos.partially_occupied_bins())
    frames = [psteros.LabelledStructure(label, built[label], energies[label], forces[label]) for label in built]
    config = psteros.MaceFineTuneConfig(args.model, replay_file=args.replay)
    train, valid, record = psteros.select_replay_structures(args.replay, elements=["Ru", "O"])
    print("replay structures:", record["train"], "train,", record["valid"], "validation")
    folder = psteros.prepare_finetune_job(f"{args.out}/finetune", frames, config, name="ruo2_m2",
                                          replay_train=train, replay_valid=valid)
    slabs = {label: built[label] for label in psteros.RUTILE_110_TERMINATIONS}
    models = {f"m2_seed{seed}": (f"models/ruo2_m2_seed{seed}.model", "Default") for seed in config.seeds}
    models["m0"] = (args.model, config.foundation_head)
    psteros.prepare_relax_job(f"{args.out}/relax", {**slabs, **{k: built[k] for k in BULK_LABELS}}, models,
                              fixed_sites={label: psteros.central_sites(slab) for label, slab in slabs.items()})
    print("job folders:", folder, f"{args.out}/relax (copy the models/ folder of the fine-tune next to it)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="phase", required=True)
    first = sub.add_parser("labels")
    first.add_argument("--profile", required=True)
    first.add_argument("--code", required=True)
    first.add_argument("--submit", action="store_true")
    second = sub.add_parser("jobs")
    second.add_argument("--profile", required=True)
    second.add_argument("--labels-pk", type=int, required=True)
    second.add_argument("--model", required=True)
    second.add_argument("--replay", required=True)
    second.add_argument("--out", default="jobs")
    args = parser.parse_args()
    from aiida import load_profile

    load_profile(args.profile)
    {"labels": labels, "jobs": jobs}[args.phase](args)


if __name__ == "__main__":
    main()

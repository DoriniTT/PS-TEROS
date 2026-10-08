"""SnO2(110) surface phase diagram with VASP: the analysis.

Reads the static energies and relaxed structures produced by ``campaign.py``
and writes the phase diagram twice:

* ``<stem>.png`` (or ``.pdf``/``.svg``): the figure drawn by psteros;
* ``<stem>.csv``: gamma(Delta mu_O) of every termination, the stable
  termination and the stability-window flag, for plotting elsewhere.

    python phase_diagram.py --profile P --refs-pk 1234 --slabs-pk 5678 --out results/sno2_110

With the graph of ``vibrations.py`` (``--vib-pk``) and a temperature, the
energies become free energies E + F_vib(T) (the O2 reference gets its
zero-point energy) and the diagram is the one at that temperature.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import psteros

REFERENCES = ("sno2_bulk", "alpha_sn", "o2")
SLABS = ("slab_o", "slab_sno", "slab_sn2o")


def build_diagram(refs_graph, slabs_graph, vibrations_graph=None, temperature_k=None) -> psteros.SurfacePhaseDiagram:
    energies, relaxed = psteros.read_vasp_results(refs_graph, REFERENCES)
    slab_energies, slabs = psteros.read_vasp_results(slabs_graph, SLABS)
    if vibrations_graph is not None:
        vibrations = psteros.read_vibrations(vibrations_graph, REFERENCES + SLABS)
        bulk = vibrations["sno2_bulk"]
        energies = {
            "sno2_bulk": psteros.solid_free_energy_ev(energies["sno2_bulk"], bulk, temperature_k),
            "alpha_sn": psteros.solid_free_energy_ev(energies["alpha_sn"], vibrations["alpha_sn"], temperature_k),
            "o2": psteros.molecule_reference_energy_ev(energies["o2"], vibrations["o2"]),
        }
        slab_energies = {
            label: psteros.solid_free_energy_ev(slab_energies[label], vibrations[label], temperature_k, bulk=bulk)
            for label in SLABS
        }
    bulk = relaxed["sno2_bulk"].get_pymatgen_structure()
    metal = relaxed["alpha_sn"].get_pymatgen_structure()
    references = psteros.BinaryOxideReferences(
        bulk_energy_ev=energies["sno2_bulk"],
        bulk_composition=bulk.composition,
        oxygen_molecule_energy_ev=energies["o2"],
        metal_energy_per_atom_ev=energies["alpha_sn"] / len(metal),
    )
    terminations = [
        psteros.SlabTermination.from_structure(label, slab_energies[label], slabs[label]) for label in SLABS
    ]
    return psteros.surface_phase_diagram(terminations, references)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--refs-pk", type=int, required=True)
    parser.add_argument("--slabs-pk", type=int, required=True)
    parser.add_argument("--out", default="sno2_110_phase_diagram", help="output path without suffix")
    parser.add_argument("--format", default="png", choices=("png", "pdf", "svg"))
    parser.add_argument("--vib-pk", type=int, help="PK of the finished vibrations.py graph (optional)")
    parser.add_argument("--temperature", type=float, help="temperature in K of the free energies (with --vib-pk)")
    args = parser.parse_args(argv)
    if (args.vib_pk is None) != (args.temperature is None):
        parser.error("--vib-pk and --temperature go together")

    from aiida import load_profile, orm

    load_profile(args.profile)
    vibrations_graph = orm.load_node(args.vib_pk) if args.vib_pk is not None else None
    diagram = build_diagram(orm.load_node(args.refs_pk), orm.load_node(args.slabs_pk), vibrations_graph, args.temperature)
    references = diagram.references
    formation = "Delta H_f" if vibrations_graph is None else f"Delta G_f({args.temperature:g} K)"
    print(f"{formation}({references.formula}) = {references.formation_enthalpy_ev:.3f} eV per formula unit")
    print(f"stability window: {references.oxygen_poor_limit_ev:.3f} <= Delta mu_O <= 0 eV")
    for delta_mu, below, above in diagram.transitions:
        print(f"transition at Delta mu_O = {delta_mu:.3f} eV: {below} -> {above}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    title = "SnO$_2$(110) surface phase diagram"
    if vibrations_graph is not None:
        title += f" at {args.temperature:g} K"
    print("figure:", diagram.plot(out.with_suffix(f".{args.format}"), title=title))
    print("data:  ", diagram.to_csv(out.with_suffix(".csv")))


if __name__ == "__main__":
    main()

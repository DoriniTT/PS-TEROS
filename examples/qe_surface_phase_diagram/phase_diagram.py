"""SnO2(110) surface phase diagram with Quantum ESPRESSO: the analysis.

Reads the static energies and relaxed structures produced by ``campaign.py``
and writes the phase diagram twice:

* ``<stem>.png`` (or ``.pdf``/``.svg``): the figure drawn by psteros;
* ``<stem>.csv``: gamma(Delta mu_O) of every termination, the stable
  termination and the stability-window flag, for plotting elsewhere.

    python phase_diagram.py --profile P --refs-pk 1234 --slabs-pk 5678 --out results/sno2_110

Several PKs may be given per option, for example a refs graph plus a follow-up
graph that recomputed some static energies.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import psteros

SLABS = ("slab_o", "slab_sno", "slab_sn2o")


def _child(graph, task_name):
    from aiida.common.links import LinkType

    for link in graph.base.links.get_outgoing(link_type=LinkType.CALL_WORK).all():
        if link.link_label == task_name:
            return link.node
    return None


def _output(graphs, label, kind):
    """Static parameters or relaxed structure of ``label`` from the first graph that has it.

    aiida-workgraph attaches no graph-level outputs when any task failed, so the
    work chain of the task (the call link carries the task name) is the fallback.
    """
    names = {
        "parameters": [(f"{label}_static_parameters", f"{label}_static_qe", "output_parameters"),
                       (f"{label}_parameters", f"{label}_qe", "output_parameters")],
        "structure": [(f"{label}_relaxed_structure", f"{label}_relax_qe", "output_structure")],
    }[kind]
    # Exit 501: vc-relax converged, final SCF above thresholds; aiida-qe keeps the structure as final.
    accepted = {0, 501} if kind == "structure" else {0}
    for graph in graphs:
        for graph_output, task_name, port in names:
            if graph_output in graph.outputs:
                return graph.outputs[graph_output]
            child = _child(graph, task_name)
            if child is not None and child.exit_status in accepted and port in child.outputs:
                return child.outputs[port]
    raise LookupError(f"no finished {kind} for {label!r} in graphs {[graph.pk for graph in graphs]}")


def _energy(graphs, label) -> float:
    return float(_output(graphs, label, "parameters")["energy"])  # eV


def build_diagram(refs_graphs, slab_graphs) -> psteros.SurfacePhaseDiagram:
    bulk = _output(refs_graphs, "sno2_bulk", "structure").get_pymatgen_structure()
    metal = _output(refs_graphs, "alpha_sn", "structure").get_pymatgen_structure()
    references = psteros.BinaryOxideReferences(
        bulk_energy_ev=_energy(refs_graphs, "sno2_bulk"),
        bulk_composition=bulk.composition,
        oxygen_molecule_energy_ev=_energy(refs_graphs, "o2"),
        metal_energy_per_atom_ev=_energy(refs_graphs, "alpha_sn") / len(metal),
    )
    terminations = [
        psteros.SlabTermination.from_structure(
            label, _energy(slab_graphs, label), _output(slab_graphs, label, "structure")
        )
        for label in SLABS
    ]
    return psteros.surface_phase_diagram(terminations, references)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--refs-pk", type=int, nargs="+", required=True)
    parser.add_argument("--slabs-pk", type=int, nargs="+", required=True)
    parser.add_argument("--out", default="sno2_110_phase_diagram", help="output path without suffix")
    parser.add_argument("--format", default="png", choices=("png", "pdf", "svg"))
    args = parser.parse_args(argv)

    from aiida import load_profile, orm

    load_profile(args.profile)
    diagram = build_diagram(
        [orm.load_node(pk) for pk in args.refs_pk], [orm.load_node(pk) for pk in args.slabs_pk]
    )
    references = diagram.references
    print(f"Delta H_f({references.formula}) = {references.formation_enthalpy_ev:.3f} eV per formula unit")
    print(f"stability window: {references.oxygen_poor_limit_ev:.3f} <= Delta mu_O <= 0 eV")
    for delta_mu, below, above in diagram.transitions:
        print(f"transition at Delta mu_O = {delta_mu:.3f} eV: {below} -> {above}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    print("figure:", diagram.plot(out.with_suffix(f".{args.format}"), title="SnO$_2$(110) surface phase diagram"))
    print("data:  ", diagram.to_csv(out.with_suffix(".csv")))


if __name__ == "__main__":
    main()

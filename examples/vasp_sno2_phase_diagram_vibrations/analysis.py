"""SnO2(110) surface phase diagram with VASP and vibrational references: the analysis.

Reads the graphs of ``campaign.py`` and writes, in ``--out``'s folder:

* ``reference_free_energies.csv``  every term of the free energy of O2, SnO2
  and alpha-Sn at 0, 298.15, 600 and 1000 K (p = 1 bar for the gas);
* ``<stem>.png`` / ``<stem>.csv``  the surface phase diagram with **DFT
  energies** of the slabs and of the bulk (the slabs have no vibrational free
  energy, so a bulk free energy inside gamma would be inconsistent, see
  ``docs/source/reference-thermochemistry.rst``);
* ``<stem>_transitions.csv``  the transitions of the diagram as Delta mu_O and,
  through ``psteros.oxygen_pressure_bar``, as p(O2) at 600 and 1000 K;
* ``oxygen_conditions.csv``  Delta mu_O of O2 gas at a few (T, p).

It also prints the O-poor limit of the stability window at 0 K (DFT) and at
600 / 1000 K (free energies of SnO2 and Sn, bare E(O2) on the axis), and
compares the stable terminations with the QE example.

    python analysis.py --refs-pk 1234 --slabs-pk 5678 --out results/sno2_110
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import psteros

SLABS = ("slab_o", "slab_sno", "slab_sn2o")
TEMPERATURES_K = (0.0, 298.15, 600.0, 1000.0)
PRESSURE_BAR = 1.0
CONDITIONS = ((298.15, 1.0), (298.15, 1e-10), (600.0, 1.0), (600.0, 1e-10), (1000.0, 1.0), (1000.0, 1e-6))
PRESSURE_READING_K = (600.0, 1000.0)
QE_RESULT = Path(__file__).resolve().parents[1] / "qe_surface_phase_diagram" / "results" / "sno2_110_phase_diagram.csv"


def _child(graph, task_name):
    from aiida.common.links import LinkType

    children = [
        link.node
        for link in graph.base.links.get_outgoing(link_type=LinkType.CALL_WORK).all()
        if link.link_label == task_name
    ]
    return max(children, key=lambda node: node.pk) if children else None


def _slab_output(graphs, label, name):
    """``misc`` or ``structure`` of a slab relaxation from the first graph that has it.

    aiida-workgraph attaches no graph-level outputs while any task is
    unfinished or failed, so the work chain of the task is the fallback.
    """

    for graph in graphs:
        if f"{label}_{name}" in graph.outputs:
            return graph.outputs[f"{label}_{name}"]
        child = _child(graph, f"{label}_vasp")
        if child is not None and child.exit_status == 0 and name in child.outputs:
            return child.outputs[name]
    raise LookupError(f"no finished relaxation of {label!r} ({name}) in graphs {[graph.pk for graph in graphs]}")


def _slab_energy_ev(misc) -> float:
    """sigma -> 0 energy of a relaxation: the quantity psteros uses for the references too."""

    energies = misc.get_dict()["total_energies"]
    return float(energies["energy_extrapolated"])


def reference_dft(refs_pks):
    """``{label: (static energy in eV, relaxed structure)}`` from the reference graphs."""

    found = {}
    for pk in refs_pks:
        for label, blocks in psteros.reference_results(pk).items():
            static, relax = blocks["static"], blocks["relax"]
            if label not in found and static["energy"] is not None and relax["structure"] is not None:
                found[label] = (static["energy"], relax["structure"].get_pymatgen_structure())
    missing = {"o2", "sno2", "sn"}.difference(found)
    if missing:
        raise SystemExit(f"no finished relax + static for {sorted(missing)} in graphs {list(refs_pks)}")
    return found


def thermochemistry(refs_pks):
    """``{label: IdealGasMolecule | HarmonicSolid}`` from the first graph that finished all three."""

    errors = []
    for pk in refs_pks:
        try:
            return psteros.reference_thermochemistry(pk)
        except ValueError as error:
            errors.append(f"{pk}: {error}")
    raise SystemExit("no reference graph with finished vibrations:\n  " + "\n  ".join(errors))


def write_free_energies(systems, path: Path) -> None:
    rows = []
    for temperature in TEMPERATURES_K:
        for label, terms in psteros.free_energies(systems, temperature, PRESSURE_BAR).items():
            rows.append({"label": label, **terms.as_dict()})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"free energies: {path}")
    for row in rows:
        print(
            f"  {row['label']:5s} T={row['temperature_K']:7.2f} K  E={row['electronic_energy_eV']:12.5f}  "
            f"ZPE={row['zero_point_energy_eV']:.4f}  G={row['free_energy_eV']:12.5f} eV"
        )


def references_at(dft, systems, temperature_k):
    """Stability window references: DFT energies at ``temperature_k=None``, else free energies."""

    (e_sno2, sno2), (e_sn, sn), (e_o2, _) = dft["sno2"], dft["sn"], dft["o2"]
    if temperature_k is not None:
        free = psteros.free_energies(systems, temperature_k, PRESSURE_BAR)
        e_sno2, e_sn = free["sno2"].free_energy_ev, free["sn"].free_energy_ev
    return psteros.BinaryOxideReferences(
        bulk_energy_ev=e_sno2,
        bulk_composition=sno2.composition,
        oxygen_molecule_energy_ev=e_o2,  # the axis keeps the bare DFT energy of O2
        metal_energy_per_atom_ev=e_sn / len(sn),
    )


def write_oxygen_conditions(oxygen, path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["temperature_K", "pressure_bar", "delta_mu_O_eV", "delta_mu_O_without_zero_point_eV"])
        for temperature, pressure in CONDITIONS:
            writer.writerow([
                temperature, pressure,
                psteros.delta_mu_oxygen_ev(oxygen, temperature, pressure),
                psteros.delta_mu_oxygen_ev(oxygen, temperature, pressure, include_zero_point=False),
            ])
            print(
                f"  O2 at {temperature:7.2f} K, {pressure:g} bar: Delta mu_O = "
                f"{psteros.delta_mu_oxygen_ev(oxygen, temperature, pressure):.3f} eV "
                f"({psteros.delta_mu_oxygen_ev(oxygen, temperature, pressure, include_zero_point=False):.3f} without ZPE)"
            )


def write_transitions(diagram, oxygen, path: Path) -> None:
    header = ["delta_mu_O_eV", "from", "to"] + [f"p_O2_bar_at_{int(t)}K" for t in PRESSURE_READING_K]
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for delta_mu, below, above in diagram.transitions:
            pressures = [psteros.oxygen_pressure_bar(oxygen, t, delta_mu) for t in PRESSURE_READING_K]
            writer.writerow([delta_mu, below, above] + pressures)
            readings = ", ".join(f"p(O2) = {p:.2e} bar at {int(t)} K" for t, p in zip(PRESSURE_READING_K, pressures))
            print(f"  transition at Delta mu_O = {delta_mu:.3f} eV: {below} -> {above} ({readings})")


def compare_with_qe(diagram) -> None:
    if not QE_RESULT.exists():
        print("QE example results not found; skipping the comparison")
        return
    with QE_RESULT.open() as handle:
        qe_rows = [row for row in csv.DictReader(handle) if row["in_stability_window"] == "True"]
    qe_low, qe_high = float(qe_rows[0]["delta_mu_O_eV"]), float(qe_rows[-1]["delta_mu_O_eV"])
    low = diagram.references.oxygen_poor_limit_ev
    print(f"stable terminations (VASP window {low:.2f}..0 eV, QE window {qe_low:.2f}..{qe_high:.2f} eV):")
    print("  fraction of the window   VASP        QE")
    for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
        vasp_index = min(range(len(diagram.delta_mu_oxygen_ev)),
                         key=lambda i: abs(diagram.delta_mu_oxygen_ev[i] - (low + fraction * -low)))
        qe_target = qe_low + fraction * (qe_high - qe_low)
        qe_row = min(qe_rows, key=lambda row: abs(float(row["delta_mu_O_eV"]) - qe_target))
        print(f"  {fraction:5.2f} (0 = O-poor)       {diagram.stable[vasp_index]:10s}  {qe_row['stable_termination']}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--profile", default="psteros_sno2_vibrations")
    parser.add_argument("--refs-pk", type=int, nargs="+", required=True, help="reference graph(s), first match wins")
    parser.add_argument("--slabs-pk", type=int, nargs="+", required=True)
    parser.add_argument("--out", default="results/sno2_110", help="output path without suffix")
    parser.add_argument("--format", default="png", choices=("png", "pdf", "svg"))
    args = parser.parse_args(argv)

    from aiida import load_profile, orm

    load_profile(args.profile)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    systems = thermochemistry(args.refs_pk)
    oxygen = systems["o2"]
    write_free_energies(systems, out.parent / "reference_free_energies.csv")
    print("oxygen chemical potential:")
    write_oxygen_conditions(oxygen, out.parent / "oxygen_conditions.csv")

    dft = reference_dft(args.refs_pk)
    references = references_at(dft, systems, None)
    print(f"Delta H_f({references.formula}) = {references.formation_enthalpy_ev:.3f} eV per formula unit (DFT, 0 K)")
    graphs = [orm.load_node(pk) for pk in args.slabs_pk]
    terminations = [
        psteros.SlabTermination.from_structure(
            label,
            _slab_energy_ev(_slab_output(graphs, label, "misc")),
            _slab_output(graphs, label, "structure"),
        )
        for label in SLABS
    ]
    diagram = psteros.surface_phase_diagram(terminations, references)
    print(f"phase diagram (DFT energies), stability window {references.oxygen_poor_limit_ev:.3f} <= Delta mu_O <= 0 eV:")
    write_transitions(diagram, oxygen, out.parent / f"{out.name}_transitions.csv")
    figure = diagram.plot(out.with_suffix(f".{args.format}"), title="SnO$_2$(110), VASP-PBE")
    print(f"figure: {figure}\ndata:   {diagram.to_csv(out.with_suffix('.csv'))}")

    print("O-poor limit of the stability window (bare E(O2) on the axis):")
    print(f"  0 K, DFT energies:  {references.oxygen_poor_limit_ev:.3f} eV")
    for temperature in (600.0, 1000.0):
        window = references_at(dft, systems, temperature)
        print(f"  {temperature:6.1f} K, free energies:  {window.oxygen_poor_limit_ev:.3f} eV "
              f"(Delta G_f = {window.formation_enthalpy_ev:.3f} eV per SnO2)")
    compare_with_qe(diagram)


if __name__ == "__main__":
    main()

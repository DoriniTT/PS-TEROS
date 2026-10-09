"""SnO2(110) surface phase diagram of one finished campaign graph: the analysis.

Reads the graph written by ``campaign.py`` (``--pk``) and writes, with ``--out``
as the path without suffix:

* ``<out>.png`` and ``<out>.csv``  the surface phase diagram of the three slabs, from
  the static energies of the slabs and of the references (DFT, 0 K);
* ``<out>_transitions.csv``        the transitions as Delta mu_O (eV) and as the O2
  pressure at 600 K and 1000 K, from ``psteros.reference_thermochemistry`` on the same graph.

    python analysis.py --pk 1234 --out results/sno2_110
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import psteros

PROFILE = "psteros_sno2_vibrations"
READING_K = (600.0, 1000.0)


def binary_references(results: dict) -> psteros.BinaryOxideReferences:
    """Stability window of SnO2: static energies, and the relaxed cells for composition and atom count."""

    refs = results["references"]

    def static_energy(label: str) -> float:
        energy = refs[label]["static"]["energy"]
        if energy is None:
            raise SystemExit(f"reference {label!r} has no static energy yet")
        return energy

    def relaxed_cell(label: str):
        structure = refs[label]["relax"]["structure"]
        if structure is None:
            raise SystemExit(f"reference {label!r} has no relaxed structure yet")
        return structure.get_pymatgen_structure()

    sno2, sn = relaxed_cell("sno2"), relaxed_cell("sn")
    return psteros.BinaryOxideReferences(
        bulk_energy_ev=static_energy("sno2"),
        bulk_composition=sno2.composition,
        oxygen_molecule_energy_ev=static_energy("o2"),  # the axis keeps the bare DFT energy of O2
        metal_energy_per_atom_ev=static_energy("sn") / len(sn),
    )


def write_transitions(diagram, oxygen, path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["delta_mu_O_eV", "from", "to"] + [f"p_O2_bar_at_{int(t)}K" for t in READING_K])
        for delta_mu, below, above in diagram.transitions:
            pressures = [psteros.oxygen_pressure_bar(oxygen, t, delta_mu) for t in READING_K]
            writer.writerow([delta_mu, below, above, *pressures])
            readings = ", ".join(f"p(O2) = {p:.2e} bar at {int(t)} K" for t, p in zip(READING_K, pressures))
            print(f"  {below} -> {above} at Delta mu_O = {delta_mu:.3f} eV ({readings})")
    print(f"transitions: {path}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--profile", default=PROFILE)
    parser.add_argument("--pk", type=int, required=True, help="PK of the finished campaign graph")
    parser.add_argument("--out", default="results/sno2_110", help="output path without suffix")
    parser.add_argument("--format", default="png", choices=("png", "pdf", "svg"))
    args = parser.parse_args(argv)

    from aiida import load_profile

    load_profile(args.profile)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    results = psteros.campaign_results(args.pk)
    if not {"o2", "sno2", "sn"} <= set(results["references"]) or not results["slabs"]:
        raise SystemExit(f"graph {args.pk} is not a full campaign (O2, SnO2, alpha-Sn and the slabs); --smoke has no slabs")
    references = binary_references(results)
    print(f"Delta H_f(SnO2) = {references.formation_enthalpy_ev:.3f} eV per formula unit (DFT, 0 K)")
    print(f"stability window: {references.oxygen_poor_limit_ev:.3f} <= Delta mu_O <= 0 eV")

    diagram = psteros.surface_phase_diagram(psteros.campaign_terminations(args.pk), references)
    figure = diagram.plot(out.with_suffix(f".{args.format}"), title="SnO$_2$(110), VASP-PBE")
    data = diagram.to_csv(out.with_suffix(".csv"))
    print(f"figure: {figure}\ndata:   {data}")

    oxygen = psteros.reference_thermochemistry(args.pk)["o2"]
    write_transitions(diagram, oxygen, out.parent / f"{out.name}_transitions.csv")


if __name__ == "__main__":
    main()

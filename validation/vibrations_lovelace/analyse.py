"""Analysis of the Lovelace test: diagrams, tables, frequencies and the automated checks.

    python analyse.py --refs-pk R --slabs-pk S --vib-pk V [--ibrion5-pk I] \\
        [--legacy-src <worktree of the commit before the vibrations>]

Writes to ``results/``:

* ``phase_diagram_total_energy.{png,csv}``  the 0 K diagram from the static total energies;
* ``phase_diagram_T<T>K.{png,csv}``         the diagrams with E + F_vib(T) for T = 0 (ZPE only), 300, 600, 900 K;
* ``gamma_table.csv``                       gamma at Delta mu_O = 0, Delta G_f, the O-poor limit and the
                                            transition Delta mu_O for the static energies and each T;
* ``frequencies_<label>.csv``, ``mode_summary.csv``  every frequency per label (imaginary modes negative);
* ``vibrational_change.png``, ``gamma_vs_delta_mu.png``  figures of the above;
* ``checks.json``                           the automated checks (4 to 10 of PLAN.md) with their numbers.

Only the public ``psteros`` API is used for the physics; the independent recomputation (check 8) parses the
retrieved ``vasprun.xml`` files with pymatgen.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

import psteros

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
REPOSITORY = HERE.parents[1]

PROFILE = "psteros_vibrations_lovelace"
REFERENCES = ("sno2_bulk", "alpha_sn", "o2")
SLABS = ("slab_o", "slab_sn2o")
LABELS = REFERENCES + SLABS
TEMPERATURES = (0.0, 300.0, 600.0, 900.0)
DISPLACEMENT = 0.01
EXPECTED_MODES = {"slab_o": 36, "slab_sn2o": 24, "sno2_bulk": 33, "alpha_sn": 21, "o2": 1}
EXPECTED_JOBS = {"slab_o": 72, "slab_sn2o": 48, "sno2_bulk": 72, "alpha_sn": 48, "o2": 12}
EV_A2_TO_MEV_A2 = 1000.0

# Categorical slots 1 and 2 of the reference palette for the two terminations; temperature is ordered
# (one hue, light to dark), so it is never encoded by an unrelated colour.
COLOR = {"slab_o": "#2a78d6", "slab_sn2o": "#eb6834"}
INK, INK_SECONDARY, GRID = "#0b0b0b", "#52514e", "#e1e0d9"
TEMPERATURE_SHADE = {0.0: "#9fb8d9", 300.0: "#6f95c9", 600.0: "#3d6fb0", 900.0: "#1d457f"}


# ------------------------------------------------------------------------------------------ reading


class ModeGraph:
    """The vibrations graph seen through its ``harmonic_modes`` calls.

    ``outputs[<label>_vibrations]`` is read from the ``harmonic_modes`` call of ``graph`` or, for a label whose
    task was skipped (RESULTS.md 5.5), of ``reruns``; every other attribute is that of ``graph``. This works for a
    graph that did not finish cleanly, whose own outputs may be missing.
    """

    def __init__(self, graph, reruns=()):
        from aiida.common.links import LinkType

        self._graph, self.outputs, self.source = graph, {}, {}
        for origin in (graph, *reruns):
            for link in origin.base.links.get_outgoing(link_type=LinkType.CALL_CALC).all():
                if link.node.process_label == "harmonic_modes" and link.node.exit_status == 0:
                    name = link.link_label if link.link_label.endswith("_vibrations") else link.link_label.split("__")[0]
                    if name not in self.outputs:
                        self.outputs[name] = link.node.outputs.result
                        self.source[name] = origin.pk

    def __getattr__(self, name):
        return getattr(self._graph, name)


def read_inputs(refs_graph, slabs_graph, vib_graph):
    """Energies, relaxed structures and the harmonic modes (imaginary modes kept as negative numbers)."""

    energies, relaxed = psteros.read_vasp_results(refs_graph, REFERENCES)
    slab_energies, slab_structures = psteros.read_vasp_results(slabs_graph, SLABS)
    energies.update(slab_energies)
    relaxed.update(slab_structures)
    modes = psteros.read_vibrations(vib_graph, LABELS, imaginary_modes="drop") if vib_graph is not None else None
    return energies, relaxed, modes


def references_and_terminations(energies, relaxed, modes=None, temperature_k=None):
    """The inputs of ``surface_phase_diagram``; with ``modes`` the energies become free energies."""

    energy = dict(energies)
    if modes is not None:
        bulk = modes["sno2_bulk"]
        energy = {
            "sno2_bulk": psteros.solid_free_energy_ev(energies["sno2_bulk"], bulk, temperature_k),
            "alpha_sn": psteros.solid_free_energy_ev(energies["alpha_sn"], modes["alpha_sn"], temperature_k),
            "o2": psteros.molecule_reference_energy_ev(energies["o2"], modes["o2"]),
            **{
                label: psteros.solid_free_energy_ev(energies[label], modes[label], temperature_k, bulk=bulk)
                for label in SLABS
            },
        }
    bulk_structure = relaxed["sno2_bulk"].get_pymatgen_structure()
    metal = relaxed["alpha_sn"].get_pymatgen_structure()
    references = psteros.BinaryOxideReferences(
        bulk_energy_ev=energy["sno2_bulk"],
        bulk_composition=bulk_structure.composition,
        oxygen_molecule_energy_ev=energy["o2"],
        metal_energy_per_atom_ev=energy["alpha_sn"] / len(metal),
    )
    terminations = [psteros.SlabTermination.from_structure(label, energy[label], relaxed[label]) for label in SLABS]
    return references, terminations


def diagram_for(energies, relaxed, modes=None, temperature_k=None):
    references, terminations = references_and_terminations(energies, relaxed, modes, temperature_k)
    return psteros.surface_phase_diagram(terminations, references)


def gamma_line(diagram, label):
    """``(gamma at Delta mu = 0, slope)`` in eV/A^2: gamma is linear in Delta mu."""

    points = diagram.curves[label]
    x0, x1 = points[0].delta_mu_ev, points[-1].delta_mu_ev
    g0, g1 = points[0].gamma_ev_per_angstrom2, points[-1].gamma_ev_per_angstrom2
    assert abs(x1) < 1e-12, "the grid must end at Delta mu_O = 0"
    return g1, (g1 - g0) / (x1 - x0)


def transition(diagram):
    """Delta mu_O where the two terminations have the same gamma (exact, window or not)."""

    (g_o, s_o), (g_r, s_r) = (gamma_line(diagram, label) for label in SLABS)
    return (g_r - g_o) / (s_o - s_r) if s_o != s_r else float("nan")


# ------------------------------------------------------------------------------------------ outputs


def write_csv(path: Path, header, rows) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def diagrams_and_table(energies, relaxed, modes):
    """Write the diagrams for the static energies and every T; return the table rows."""

    RESULTS.mkdir(exist_ok=True)
    rows, diagrams = [], {}
    cases = [("total_energy", None)] + [(f"T{int(t)}K", t) for t in TEMPERATURES]
    for stem, temperature in cases:
        if temperature is not None and modes is None:
            continue
        diagram = diagram_for(energies, relaxed, modes if temperature is not None else None, temperature)
        diagrams[stem] = diagram
        references = diagram.references
        title = "SnO$_2$(110), static total energies" if temperature is None else (
            f"SnO$_2$(110), E + F$_{{vib}}$ at {temperature:g} K" + (" (ZPE only)" if temperature == 0 else "")
        )
        diagram.plot(RESULTS / f"phase_diagram_{stem}.png", title=title)
        diagram.to_csv(RESULTS / f"phase_diagram_{stem}.csv")
        row = {
            "case": stem,
            "temperature_K": "" if temperature is None else temperature,
            "formation_free_energy_eV": references.formation_enthalpy_ev,
            "oxygen_poor_limit_eV": references.oxygen_poor_limit_ev,
        }
        for label in SLABS:
            gamma0, _ = gamma_line(diagram, label)
            row[f"gamma_{label}_at_dmu0_eV_per_A2"] = gamma0
            row[f"gamma_{label}_at_dmu0_J_per_m2"] = gamma0 * psteros.EV_PER_ANGSTROM2_TO_J_PER_M2
        crossing = transition(diagram)
        row["transition_dmu_O_eV"] = crossing
        row["transition_in_stability_window"] = bool(
            references.oxygen_poor_limit_ev <= crossing <= 0.0
        )
        row["stable_at_dmu0"] = diagram.stable[-1]
        rows.append(row)
    write_csv(RESULTS / "gamma_table.csv", list(rows[0]), [list(row.values()) for row in rows])
    return diagrams, rows


def frequency_tables(vib_graph):
    """One CSV per label with every frequency, and a summary of the modes."""

    summary = []
    for label in LABELS:
        data = vib_graph.outputs[f"{label}_vibrations"].get_dict()
        frequencies = sorted(data["frequencies_cm1"])
        write_csv(
            RESULTS / f"frequencies_{label}.csv",
            ["mode", "frequency_cm1"],
            [[index, f"{value:.4f}"] for index, value in enumerate(frequencies)],
        )
        modes = psteros.HarmonicVibrations.from_dict(data, imaginary_modes="drop")
        summary.append({
            "label": label,
            "modes": len(frequencies),
            "zero_modes_removed": data["zero_modes"],
            "displaced_sites": len(data["displaced_sites"]),
            "supercell_size": data["supercell_size"],
            "lowest_cm1": frequencies[0],
            "highest_cm1": frequencies[-1],
            "imaginary_modes": len(modes.imaginary_frequencies_cm1),
            "zpe_eV": modes.zero_point_energy_ev,
            "f_vib_300K_eV": modes.free_energy_ev(300.0),
            "f_vib_900K_eV": modes.free_energy_ev(900.0),
        })
    write_csv(RESULTS / "mode_summary.csv", list(summary[0]), [list(row.values()) for row in summary])
    return summary


def figures(diagrams, rows):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 10, "axes.edgecolor": GRID, "axes.labelcolor": INK_SECONDARY,
                         "xtick.color": INK_SECONDARY, "ytick.color": INK_SECONDARY})
    temperatures = [t for t in TEMPERATURES if f"T{int(t)}K" in diagrams]
    if not temperatures:
        return

    # gamma(Delta mu_O) of each termination at every T, one panel per termination, shared axes.
    window = (min(d.references.oxygen_poor_limit_ev for d in diagrams.values()), 0.0)
    xs = np.linspace(*window, 50)
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), sharey=True, constrained_layout=True)
    for axis, label in zip(axes, SLABS):
        for stem, diagram in diagrams.items():
            g0, slope = gamma_line(diagram, label)
            kinds = {"color": TEMPERATURE_SHADE[float(stem[1:-1])]} if stem != "total_energy" else {
                "color": INK, "linestyle": (0, (4, 3)), "linewidth": 1.0}
            name = "static E" if stem == "total_energy" else f"{stem[1:-1]} K"
            axis.plot(xs, (g0 + slope * xs) * psteros.EV_PER_ANGSTROM2_TO_J_PER_M2, linewidth=1.6, **kinds)
            axis.annotate(name, (xs[-1], (g0 + slope * xs[-1]) * psteros.EV_PER_ANGSTROM2_TO_J_PER_M2),
                          xytext=(4, 0), textcoords="offset points", fontsize=8, color=INK_SECONDARY, va="center")
        axis.set_title(f"{label} ({'Sn6O12' if label == 'slab_o' else 'Sn6O8'})", color=INK, fontsize=10)
        axis.set_xlabel(r"$\Delta\mu_\mathrm{O}$ (eV)")
        axis.grid(color=GRID, linewidth=0.6)
        for side in ("top", "right"):
            axis.spines[side].set_visible(False)
        axis.margins(x=0.0)
        axis.set_xlim(window[0], 0.55)
    axes[0].set_ylabel(r"$\gamma$ (J m$^{-2}$)")
    fig.savefig(RESULTS / "gamma_vs_delta_mu.png", dpi=160)
    plt.close(fig)

    # What the vibrations change: gamma at Delta mu_O = 0 (meV/A^2) and the transition / poor limit (eV) vs T.
    by_case = {row["case"]: row for row in rows}
    static = by_case["total_energy"]
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), constrained_layout=True)
    temps = [t for t in temperatures]
    for label in SLABS:
        key = f"gamma_{label}_at_dmu0_eV_per_A2"
        change = [(by_case[f"T{int(t)}K"][key] - static[key]) * EV_A2_TO_MEV_A2 for t in temps]
        axes[0].plot(temps, change, marker="o", markersize=5, linewidth=1.8, color=COLOR[label], label=label)
        axes[0].annotate(label, (temps[-1], change[-1]), xytext=(5, 0), textcoords="offset points",
                         fontsize=8, color=INK_SECONDARY, va="center")
    axes[0].set_ylabel(r"$\gamma(T) - \gamma_\mathrm{static}$ at $\Delta\mu_\mathrm{O}=0$ (meV Å$^{-2}$)")
    for key, name, color in (("transition_dmu_O_eV", "transition", "#2a78d6"), ("oxygen_poor_limit_eV", "O-poor limit", "#eb6834")):
        values = [by_case[f"T{int(t)}K"][key] - static[key] for t in temps]
        axes[1].plot(temps, values, marker="o", markersize=5, linewidth=1.8, color=color, label=name)
        axes[1].annotate(name, (temps[-1], values[-1]), xytext=(5, 0), textcoords="offset points",
                         fontsize=8, color=INK_SECONDARY, va="center")
    axes[1].set_ylabel(r"shift of $\Delta\mu_\mathrm{O}$ relative to static (eV)")
    for axis in axes:
        axis.set_xlabel("T (K)")
        axis.grid(color=GRID, linewidth=0.6)
        axis.legend(frameon=False, fontsize=8, loc="best")
        axis.set_xlim(-30, 1100)
        for side in ("top", "right"):
            axis.spines[side].set_visible(False)
    fig.savefig(RESULTS / "vibrational_change.png", dpi=160)
    plt.close(fig)


# ------------------------------------------------------------------------------------------ checks


def check_mode_counts(summary):
    detail = {row["label"]: row["modes"] for row in summary}
    passed = all(detail[label] == count for label, count in EXPECTED_MODES.items())
    zero = {row["label"]: row["zero_modes_removed"] for row in summary}
    passed &= zero["o2"] == 5 and zero["sno2_bulk"] == 3 and zero["alpha_sn"] == 3 and zero["slab_o"] == 0
    return passed, {"modes": detail, "expected": EXPECTED_MODES, "zero_modes_removed": zero}


def check_o2_stretch(summary, ibrion5_frequency=None):
    o2 = next(row for row in summary if row["label"] == "o2")
    frequency = o2["highest_cm1"]
    passed = 1500.0 <= frequency <= 1650.0
    detail = {"o2_stretch_cm1": frequency}
    if ibrion5_frequency is not None:
        detail["ibrion5_cm1"] = ibrion5_frequency
        detail["difference_cm1"] = frequency - ibrion5_frequency
        passed &= abs(frequency - ibrion5_frequency) <= 10.0
    return passed, detail


def check_solids(vib_graph, summary):
    """Bulk modes in 0-800 cm-1, imaginary modes < 50i; read_vibrations raises by default and drops on request."""

    detail = {}
    passed = True
    bulk = next(row for row in summary if row["label"] == "sno2_bulk")
    detail["sno2_bulk_range_cm1"] = [bulk["lowest_cm1"], bulk["highest_cm1"]]
    passed &= bulk["highest_cm1"] <= 800.0
    for row in summary:
        label = row["label"]
        data = vib_graph.outputs[f"{label}_vibrations"].get_dict()
        imaginary = sorted(-value for value in data["frequencies_cm1"] if value < 0)
        entry = {"imaginary_modes_cm1": imaginary}
        try:
            psteros.read_vibrations(vib_graph, [label])
            entry["default_read"] = "ok"
            default_raised = False
        except ValueError as error:
            entry["default_read"] = f"raises ValueError: {str(error)[:90]}"
            default_raised = True
        dropped = psteros.read_vibrations(vib_graph, [label], imaginary_modes="drop")[label]
        entry["drop_modes"] = len(dropped.mode_energies_ev)
        passed &= default_raised == bool(imaginary)  # raises exactly when there are imaginary modes
        if imaginary:
            passed &= max(imaginary) < 50.0
        detail[label] = entry
    return passed, detail


_KEY = re.compile(r"s\d+_[xyz]_(?:plus|minus)")


def retrieved_folders(vib_graph, label):
    """``{displacement key: FolderData}`` of the calculations feeding ``<label>_vibrations``."""

    from aiida import orm
    from aiida.common.links import LinkType

    for link in vib_graph.base.links.get_outgoing(link_type=LinkType.CALL_CALC).all():
        node = link.node
        if node.process_label == "harmonic_modes" and link.link_label.startswith(f"{label}_vibrations"):
            folders = {}
            for incoming in node.base.links.get_incoming(node_class=orm.FolderData).all():
                match = _KEY.search(incoming.link_label)
                if match:
                    folders[match.group(0)] = incoming.node
            return node, folders
    raise LookupError(f"no harmonic_modes call found for {label}")


def check_independent_recomputation(vib_graph, label="slab_sn2o"):
    """Recompute ``label`` from the retrieved vasprun.xml files with pymatgen: forces, masses, Hessian."""

    from pymatgen.core import Element
    from pymatgen.io.vasp import Vasprun

    node, folders = retrieved_folders(vib_graph, label)
    options = node.inputs.settings.get_dict()
    structure = node.inputs.structure
    sites = list(options["displaced_sites"])
    symbols = [structure.get_kind(site.kind_name).symbol for site in structure.sites]
    masses = [Element(symbol).atomic_mass.real for symbol in symbols]

    def forces_of(key):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "vasprun.xml"
            path.write_text(folders[key].get_object_content("vasprun.xml"))
            run = Vasprun(str(path), parse_dos=False, parse_eigen=False, parse_potcar_file=False)
            assert [str(s) for s in run.final_structure.species] == symbols
            return np.array(run.ionic_steps[-1]["forces"])

    plus = [forces_of(f"s{site}_{axis}_plus") for site in sites for axis in "xyz"]
    minus = [forces_of(f"s{site}_{axis}_minus") for site in sites for axis in "xyz"]
    modes = psteros.harmonic_vibrations_from_forces(
        masses_amu=masses, displaced_sites=sites, displacement_angstrom=float(options["displacement_angstrom"]),
        forces_plus=plus, forces_minus=minus, zero_modes=int(vib_graph.outputs[f"{label}_vibrations"].get_dict()["zero_modes"]),
        imaginary_modes="drop",
    )
    graph_frequencies = np.sort(np.array(vib_graph.outputs[f"{label}_vibrations"].get_dict()["frequencies_cm1"]))
    difference = float(np.max(np.abs(np.sort(np.array(modes.frequencies_cm1)) - graph_frequencies)))
    return difference < 0.1, {"label": label, "modes": len(graph_frequencies), "max_abs_difference_cm1": difference}


def check_displacements(vib_graph, relaxed, expected_incar):
    """Every displaced structure moves one coordinate of one site by +-0.01 A; no selective dynamics; INCAR = recipe."""

    from aiida import orm
    from aiida.common.links import LinkType

    problems, displaced, vasp_checked = [], 0, 0
    calls = vib_graph.base.links.get_outgoing(link_type=LinkType.CALL_CALC).all()
    by_pk = {}
    for link in calls:
        node = link.node
        if node.process_label != "displace_site":
            continue
        displaced += 1
        before, after = node.inputs.structure, node.outputs.result
        site, vector = node.inputs.site.value, node.inputs.vector.get_list()
        a, b = np.array([s.position for s in before.sites]), np.array([s.position for s in after.sites])
        delta = b - a
        moved = np.argwhere(np.abs(delta) > 1e-9)
        ok = len(moved) == 1 and moved[0][0] == site and abs(abs(delta[tuple(moved[0])]) - DISPLACEMENT) < 1e-9
        ok &= np.allclose(np.array(before.cell), np.array(after.cell)) and [s.kind_name for s in before.sites] == [
            s.kind_name for s in after.sites]
        if not ok:
            problems.append(f"{link.link_label}: moved {moved.tolist()}")
        by_pk[after.pk] = (link.link_label, before)
    for link in vib_graph.base.links.get_outgoing(link_type=LinkType.CALL_WORK).all():
        chain = link.node
        if "_vib_" not in link.link_label:
            continue
        label = link.link_label.split("_vib_")[0]
        structure = chain.inputs.structure
        if structure.pk not in by_pk:
            problems.append(f"{link.link_label}: structure is not an output of displace_site")
        parameters = chain.inputs.parameters.get_dict()
        if "dynamics" in parameters:
            problems.append(f"{link.link_label}: selective dynamics in the parameters")
        incar = {key.lower(): value for key, value in parameters["incar"].items()}
        if incar != {key.lower(): value for key, value in expected_incar(label).items()}:
            problems.append(f"{link.link_label}: INCAR {incar}")
        if chain.inputs.options.get_dict().get("import_sys_environment", True) is not False:
            problems.append(f"{link.link_label}: import_sys_environment not False")
        for calc in chain.called:
            if not hasattr(calc, "outputs") or calc.process_label != "VaspCalculation":
                continue
            poscar = calc.base.repository.get_object_content("POSCAR")
            if "selective" in poscar.lower():
                problems.append(f"{link.link_label}: Selective dynamics in the POSCAR")
            written = calc.base.repository.get_object_content("INCAR")
            tags = {line.split("=")[0].strip().lower() for line in written.splitlines() if "=" in line}
            if "nsw" not in tags or not re.search(r"(?im)^\s*NSW\s*=\s*0\b", written):
                problems.append(f"{link.link_label}: NSW != 0 in the INCAR")
            vasp_checked += 1
    return not problems, {"displace_site_calls": displaced, "vasp_inputs_checked": vasp_checked, "problems": problems[:20]}


def check_backward_compatibility(refs_graph, slabs_graph, energies, relaxed, legacy_src):
    """The 0 K diagram from total energies equals the one of the pre-feature code and of the example's logic."""

    ours = RESULTS / "phase_diagram_total_energy.csv"
    detail = {}
    # (a) the example's phase_diagram.py logic without --vib-pk, restricted to the two terminations
    sys.path.insert(0, str(REPOSITORY / "examples" / "vasp_surface_phase_diagram"))
    import phase_diagram as example

    example.SLABS = SLABS
    reference = example.build_diagram(refs_graph, slabs_graph)
    example_csv = RESULTS / "_example_total_energy.csv"
    reference.to_csv(example_csv)
    detail["identical_to_example_phase_diagram_py"] = example_csv.read_bytes() == ours.read_bytes()
    example_csv.unlink()
    # (b) the code of the commit before the vibrations, in a separate process, fed with the same numbers
    if legacy_src is not None:
        payload = {
            "energies": energies,
            "bulk_composition": dict(relaxed["sno2_bulk"].get_pymatgen_structure().composition.get_el_amt_dict()),
            "metal_atoms": len(relaxed["alpha_sn"].get_pymatgen_structure()),
            "slabs": {
                label: {
                    "composition": dict(relaxed[label].get_pymatgen_structure().composition.get_el_amt_dict()),
                    "area": float(np.linalg.norm(np.cross(*relaxed[label].get_pymatgen_structure().lattice.matrix[:2]))),
                }
                for label in SLABS
            },
        }
        outputs = {}
        for name, source in (("legacy", Path(legacy_src)), ("current", REPOSITORY)):
            target = RESULTS / f"_{name}_total_energy.csv"
            subprocess.run(
                [sys.executable, str(HERE / "legacy_diagram.py"), str(target)], input=json.dumps(payload).encode(),
                env={**__import__("os").environ, "PYTHONPATH": str(source)}, check=True, capture_output=True,
            )
            outputs[name] = target.read_bytes()
            target.unlink()
        detail["legacy_src"] = str(legacy_src)
        detail["identical_to_pre_feature_code"] = outputs["legacy"] == outputs["current"]
        detail["legacy_equals_analysis_csv"] = outputs["legacy"] == ours.read_bytes()
    passed = all(value for key, value in detail.items() if key.startswith(("identical", "legacy_equals")))
    return passed, detail


def check_physical_sanity(rows):
    by_case = {row["case"]: row for row in rows}
    static = by_case["total_energy"]
    detail = {}
    worst_gamma = max(
        abs(by_case[f"T{int(t)}K"]["gamma_slab_o_at_dmu0_eV_per_A2"] - static["gamma_slab_o_at_dmu0_eV_per_A2"])
        * EV_A2_TO_MEV_A2
        for t in TEMPERATURES
    )
    shift = max(
        abs(by_case[f"T{int(t)}K"]["transition_dmu_O_eV"] - static["transition_dmu_O_eV"]) for t in TEMPERATURES
    )
    detail["max_gamma_slab_o_change_meV_per_A2_up_to_900K"] = worst_gamma
    detail["max_transition_shift_eV"] = shift
    detail["oxygen_poor_limit_eV"] = {case: row["oxygen_poor_limit_eV"] for case, row in by_case.items()}
    detail["formation_free_energy_eV"] = {case: row["formation_free_energy_eV"] for case, row in by_case.items()}
    # DeltaG_f(T) moves the poor limit: the limit must differ between the static case and T > 0
    moves = abs(by_case["T900K"]["oxygen_poor_limit_eV"] - static["oxygen_poor_limit_eV"]) > 1e-6
    detail["poor_limit_moves_with_T"] = bool(moves)
    return bool(worst_gamma <= 10.0 and shift <= 0.2 and moves), detail


def ibrion5_frequency(graph):
    """The O2 stretch (cm-1) of a finished IBRION=5 job, from its OUTCAR or the dynamical matrix of vasprun.xml."""

    folder = graph.outputs.o2_retrieved
    names = folder.list_object_names()
    if "OUTCAR" in names:
        text = folder.get_object_content("OUTCAR")
        values = [float(m.group(1)) for m in re.finditer(r"f\s+=\s+[\d.]+ THz\s+[\d.]+ 2PiTHz\s+([\d.]+) cm-1", text)]
        if values:
            return max(values)
    from pymatgen.io.vasp import Vasprun

    with tempfile.TemporaryDirectory() as folder_name:
        path = Path(folder_name) / "vasprun.xml"
        path.write_text(folder.get_object_content("vasprun.xml"))
        run = Vasprun(str(path), parse_dos=False, parse_eigen=False, parse_potcar_file=False)
        eigenvalues = np.array(run.normalmode_eigenvals)
        return float(np.max(np.sqrt(np.abs(eigenvalues)) * 521.47))
    raise ValueError("no frequencies in the IBRION=5 job")


# ------------------------------------------------------------------------------------------ main


def expected_incar_factory():
    sys.path.insert(0, str(HERE))
    import run_test

    def expected(label):
        incar = {**run_test.ELECTRONIC, **run_test.STATIC, "EDIFF": run_test.VIBRATIONS_EDIFF}
        if label == "o2":
            incar.update({"ISPIN": 2, "MAGMOM": [1.0, 1.0]})
        return incar

    return expected


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--profile", default=PROFILE)
    parser.add_argument("--refs-pk", type=int, required=True)
    parser.add_argument("--slabs-pk", type=int, required=True)
    parser.add_argument("--vib-pk", type=int, help="omit to analyse the total-energy diagram only")
    parser.add_argument("--ibrion5-pk", type=int)
    parser.add_argument("--rerun-pk", type=int, nargs="*", default=[],
                        help="graphs from rerun_alpha_sn_modes.py that supply the skipped mode tasks")
    parser.add_argument("--legacy-src", help="checkout of the commit before the vibrations (check 9)")
    args = parser.parse_args(argv)

    from aiida import load_profile, orm

    load_profile(args.profile)
    refs, slabs = orm.load_node(args.refs_pk), orm.load_node(args.slabs_pk)
    vib = ModeGraph(orm.load_node(args.vib_pk), [orm.load_node(pk) for pk in args.rerun_pk]) if args.vib_pk else None
    if vib is not None:
        print("mode results from:", vib.source)
    RESULTS.mkdir(exist_ok=True)

    energies, relaxed, modes = read_inputs(refs, slabs, vib)
    print("static energies (eV):", {label: round(value, 6) for label, value in energies.items()})
    diagrams, rows = diagrams_and_table(energies, relaxed, modes)
    checks = {}
    checks["9_backward_compatibility"] = check_backward_compatibility(refs, slabs, energies, relaxed, args.legacy_src)
    if vib is not None:
        summary = frequency_tables(vib)
        figures(diagrams, rows)
        frequency = ibrion5_frequency(orm.load_node(args.ibrion5_pk)) if args.ibrion5_pk else None
        checks["4_displacements"] = check_displacements(vib, relaxed, expected_incar_factory())
        checks["5_mode_counts"] = check_mode_counts(summary)
        checks["6_o2_stretch"] = check_o2_stretch(summary, frequency)
        checks["7_solids"] = check_solids(vib, summary)
        checks["8_independent_recomputation"] = check_independent_recomputation(vib)
        checks["10_physical_sanity"] = check_physical_sanity(rows)
    for name, (passed, detail) in checks.items():
        print(f"{'PASS' if passed else 'FAIL'}  {name}: {json.dumps(detail, default=float)[:600]}")
    (RESULTS / "checks.json").write_text(
        json.dumps({name: {"passed": bool(p), "detail": d} for name, (p, d) in checks.items()}, indent=2, default=float)
    )
    for row in rows:
        print({key: (round(value, 5) if isinstance(value, float) else value) for key, value in row.items()})
    return 0 if all(passed for passed, _ in checks.values()) else 1


if __name__ == "__main__":
    sys.exit(main())

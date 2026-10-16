"""Run a failed or skipped ``<label>_vibrations`` task of a vibrations graph as its own small graph.

The WorkGraph engine can start (or skip) a mode task although every displaced calculation finished with exit status 0:
after a host restart it marked one displacement FAILED and skipped the mode task (SnO2 test, graph 3400), and with two daemon
workers it started ``srtio3_bulk_vibrations`` with 58 of 60 retrieved folders (ternary test, graph 12807). This script takes
the inputs the task should have had (the structure and the settings of the task, the retrieved folder of every
displacement, taken from the finished work chains) and runs the same ``harmonic_modes`` calcfunction in a new graph, so the
modes keep their provenance.  It refuses to run unless every displaced calculation of the label finished with exit status 0.

    python rerun_modes.py --vib-pk 12807 --label srtio3_bulk [--submit]
"""

from __future__ import annotations

import argparse
import re

PROFILE = "psteros_vibrations_lovelace"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vib-pk", type=int, required=True)
    parser.add_argument("--label", required=True, help="structure label, e.g. srtio3_bulk")
    parser.add_argument("--profile", default=PROFILE)
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args(argv)

    from aiida import load_profile, orm
    from aiida.common.links import LinkType
    from aiida_workgraph import WorkGraph

    from psteros.backends.vibration_tasks import harmonic_modes

    label = args.label
    load_profile(args.profile)
    graph = orm.load_node(args.vib_pk)
    original = WorkGraph.load(args.vib_pk).tasks[f"{label}_vibrations"]
    structure = original.inputs.structure.value
    settings = original.inputs.settings.value
    key = re.compile(rf"^{label}_vib_(s\d+_[xyz]_(?:plus|minus))_vasp$")
    retrieved = {}
    for link in graph.base.links.get_outgoing(link_type=LinkType.CALL_WORK).all():
        match = key.match(link.link_label)
        if match:
            if link.node.exit_status != 0:
                raise SystemExit(f"{link.link_label} (PK {link.node.pk}) did not finish with exit status 0")
            retrieved[match.group(1)] = link.node.outputs.retrieved
    expected = 6 * len(settings.get_dict()["displaced_sites"])
    print(f"{label}: {len(retrieved)} retrieved folders of {expected}; structure PK {structure.pk}, settings PK {settings.pk}")
    if len(retrieved) != expected:
        raise SystemExit("not every displaced calculation of the label has finished with exit status 0")

    rerun = WorkGraph(name=f"{label}_vibrations_rerun")
    task = rerun.add_task(
        harmonic_modes, name=f"{label}_vibrations", structure=structure, settings=settings, retrieved=retrieved
    )
    rerun.outputs.__setattr__(f"{label}_vibrations", task.outputs.result)
    if args.submit:
        rerun.submit()
        print(f"submitted: PK={rerun.pk}")
    else:
        print("built; use --submit")


if __name__ == "__main__":
    main()

"""Run the skipped ``alpha_sn_vibrations`` task of the vibrations graph as its own small graph.

After a host restart the WorkGraph engine marked the task ``alpha_sn_vib_s0_x_minus_vasp`` FAILED although its work chain
and VASP calculation finished with exit status 0, so the engine skipped ``alpha_sn_vibrations`` (RESULTS.md 5.5).  All 48
alpha-Sn calculations are fine.  This script takes the inputs the skipped task would have had (the structure and the
settings of the task, the retrieved folder of every displacement) and runs the same ``harmonic_modes`` calcfunction in a
new graph, so the modes keep their provenance.

    python rerun_alpha_sn_modes.py --vib-pk 3400 [--submit]
"""

from __future__ import annotations

import argparse
import re

PROFILE = "psteros_vibrations_lovelace"
LABEL = "alpha_sn"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vib-pk", type=int, required=True)
    parser.add_argument("--profile", default=PROFILE)
    parser.add_argument("--submit", action="store_true")
    args = parser.parse_args(argv)

    from aiida import load_profile, orm
    from aiida.common.links import LinkType
    from aiida_workgraph import WorkGraph

    from psteros.backends.vibration_tasks import harmonic_modes

    load_profile(args.profile)
    graph = orm.load_node(args.vib_pk)
    original = WorkGraph.load(args.vib_pk).tasks[f"{LABEL}_vibrations"]
    structure = original.inputs.structure.value
    settings = original.inputs.settings.value
    key = re.compile(rf"^{LABEL}_vib_(s\d+_[xyz]_(?:plus|minus))_vasp$")
    retrieved = {}
    for link in graph.base.links.get_outgoing(link_type=LinkType.CALL_WORK).all():
        match = key.match(link.link_label)
        if match:
            if link.node.exit_status != 0:
                raise SystemExit(f"{link.link_label} (PK {link.node.pk}) did not finish with exit status 0")
            retrieved[match.group(1)] = link.node.outputs.retrieved
    expected = 6 * len(settings.get_dict()["displaced_sites"])
    print(f"{LABEL}: {len(retrieved)} retrieved folders of {expected}; structure PK {structure.pk}, settings PK {settings.pk}")
    if len(retrieved) != expected:
        raise SystemExit("not every displaced calculation of alpha_sn has finished with exit status 0")

    rerun = WorkGraph(name="alpha_sn_vibrations_rerun")
    task = rerun.add_task(
        harmonic_modes, name=f"{LABEL}_vibrations", structure=structure, settings=settings, retrieved=retrieved
    )
    rerun.outputs.__setattr__(f"{LABEL}_vibrations", task.outputs.result)
    if args.submit:
        rerun.submit()
        print(f"submitted: PK={rerun.pk}")
    else:
        print("built; use --submit")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""
Post-hoc Perturbative Expansion Analysis

This example demonstrates how to run perturbative expansion analysis with
different probe charge (q) and electron transfer (DeltaN) values using
results from a previous Phase 2 workflow.

This is useful for:
- Exploring different adsorbate charges without re-running VASP
- Comparing cationic vs anionic adsorbate behavior
- Generating energy maps for different electron transfer scenarios

Prerequisites:
- Completed Phase 2 workflow with LOCPOT and LOCPOT_FUKUI.vasp outputs
- Or manually provide these files from FukuiGrid calculations

Usage:
    source ~/envs/aiida/bin/activate
    python run_posthoc_analysis.py <previous_workflow_pk>

Example:
    python run_posthoc_analysis.py 12345

    Or interactively:
    >>> from aiida import orm, load_profile
    >>> load_profile()
    >>> from teros.core.fukui import run_perturbative_expansion_calcfunc
    >>> locpot = orm.load_node(12345)
    >>> fukui_pot = orm.load_node(12346)
    >>> result = run_perturbative_expansion_calcfunc(
    ...     locpot_neutral=locpot,
    ...     fukui_potential=fukui_pot,
    ...     probe_charge=orm.Float(0.5),
    ...     electron_transfer=orm.Float(-0.5),
    ... )
"""

import sys
from pathlib import Path
from aiida import orm, load_profile
from teros.core.fukui import (
    get_fukui_results,
    print_fukui_summary,
    run_perturbative_expansion_calcfunc,
)


def run_posthoc_analysis(workflow_pk: int):
    """
    Run post-hoc perturbative expansion analysis with different parameters.

    Args:
        workflow_pk: PK of a completed Phase 2 (or higher) Fukui workflow
    """
    print("\n" + "=" * 70)
    print("POST-HOC PERTURBATIVE EXPANSION ANALYSIS")
    print("=" * 70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='presto')
    print("   Profile loaded")

    # Load results from previous workflow
    print(f"\n2. Loading results from workflow PK={workflow_pk}...")
    try:
        results = get_fukui_results(workflow_pk)
    except Exception as e:
        print(f"   ERROR: Could not load workflow: {e}")
        return 1

    # Check required outputs
    if results.get('fukui_potential') is None:
        print("   ERROR: Workflow does not have fukui_potential output")
        print("   Make sure the workflow was run with compute_fukui_potential=True")
        return 1

    # Get LOCPOT - either from locpot_neutral or we need to explain how to get it
    locpot_neutral = results.get('locpot_neutral')
    fukui_potential = results['fukui_potential']

    if locpot_neutral is None:
        print("   WARNING: locpot_neutral not found in workflow outputs")
        print("   This workflow was run without compute_perturbative_expansion=True")
        print("   You need to manually provide the LOCPOT file from the neutral calculation")
        print("\n   To get the LOCPOT, check the neutral (delta_n=0.0) calculation:")
        print("   >>> verdi process report <workflow_pk>")
        print("   >>> # Find the delta_0_00 calculation PK")
        print("   >>> verdi calcjob inputcat <calcjob_pk> -p INCAR")
        print("\n   If LOCPOT was not retrieved, you'll need to re-run with:")
        print("   >>> compute_perturbative_expansion=True")
        return 1

    print(f"   LOCPOT (neutral): PK={locpot_neutral.pk}")
    print(f"   Fukui potential: PK={fukui_potential.pk}")

    # Print original workflow summary
    print("\n3. Original workflow summary:")
    print("-" * 50)
    print_fukui_summary(workflow_pk)

    # Define different q/DeltaN scenarios to explore
    print("\n4. Running perturbative expansion with different parameters...")
    print("-" * 70)

    scenarios = [
        # (q, DeltaN, description)
        (0.3, -0.3, "Cationic adsorbate (e.g., Li+, Na+)"),
        (-0.3, 0.3, "Anionic adsorbate (e.g., Cl-, OH-)"),
        (0.5, -0.5, "Strongly cationic (e.g., Ca2+)"),
        (-0.5, 0.5, "Strongly anionic (e.g., O2-)"),
        (0.1, -0.1, "Weakly charged cation"),
        (1.0, -1.0, "Full electron transfer (oxidation)"),
    ]

    output_files = []

    for q, delta_n, description in scenarios:
        print(f"\n   Scenario: {description}")
        print(f"   q = {q}, DeltaN = {delta_n}")

        try:
            result = run_perturbative_expansion_calcfunc(
                locpot_neutral=locpot_neutral,
                fukui_potential=fukui_potential,
                probe_charge=orm.Float(q),
                electron_transfer=orm.Float(delta_n),
            )
            print(f"   Output: PK={result.pk}, filename={result.filename}")

            # Generate output filename
            q_str = f"{q:.1f}".replace('.', 'p').replace('-', 'm')
            dn_str = f"{delta_n:.1f}".replace('.', 'p').replace('-', 'm')
            output_name = f"MODELPOT_q{q_str}_dN{dn_str}.vasp"
            output_files.append((result, output_name, description))

        except Exception as e:
            print(f"   ERROR: {e}")

    # Export results
    print("\n5. Exporting model potential files...")
    print("-" * 50)

    output_dir = Path(__file__).parent / 'posthoc_results'
    output_dir.mkdir(exist_ok=True)

    for result, output_name, description in output_files:
        output_path = output_dir / output_name
        try:
            content = result.get_content()
            mode = 'wb' if isinstance(content, bytes) else 'w'
            with open(output_path, mode) as f:
                f.write(content)
            print(f"   Saved: {output_path.name}")
        except Exception as e:
            print(f"   ERROR saving {output_name}: {e}")

    print(f"\n   All files saved to: {output_dir}")

    # Print visualization instructions
    print("\n" + "=" * 70)
    print("VISUALIZATION INSTRUCTIONS")
    print("=" * 70)
    print("""
Open the MODELPOT_*.vasp files in VESTA:

1. File > Open > select MODELPOT_*.vasp
2. Objects > Properties > Isosurfaces
3. Add isosurface with appropriate isovalue (try 0.1-0.5 eV)
4. Use different colors for positive (red) and negative (blue) isosurfaces

Interpretation:
- Negative values (blue): Favorable adsorption sites
  The adsorbate gains energy by being at these locations
- Positive values (red): Unfavorable adsorption sites
  The adsorbate loses energy (repelled) from these locations

Compare different scenarios:
- Cationic adsorbates prefer different sites than anionic ones
- Stronger charges show more pronounced site preferences
- The interplay of q and DeltaN determines the energy landscape
""")
    print("=" * 70)

    return 0


def main():
    """Main entry point."""
    if len(sys.argv) < 2:
        print("Usage: python run_posthoc_analysis.py <workflow_pk>")
        print("\nExample:")
        print("  python run_posthoc_analysis.py 12345")
        print("\nThe workflow must have been run with compute_fukui_potential=True")
        print("Ideally also with compute_perturbative_expansion=True to have LOCPOT available")

        # Try to read PK from last_pk.txt in parent directories
        for parent in [Path(__file__).parent, Path(__file__).parent.parent / "02_Electrodes"]:
            pk_file = parent / 'last_pk.txt'
            if pk_file.exists():
                pk = int(pk_file.read_text().strip())
                print(f"\nFound last_pk.txt with PK={pk}")
                print(f"Run: python run_posthoc_analysis.py {pk}")
                break

        return 1

    workflow_pk = int(sys.argv[1])
    return run_posthoc_analysis(workflow_pk)


if __name__ == '__main__':
    sys.exit(main())

#!/home/thiagotd/envs/aiida/bin/python
"""
Example demonstrating electronic properties (DOS and bands) calculation for relaxed bulk.

This script demonstrates the ELECTRONIC PROPERTIES workflow for Ag2O:
1. Uses default builders for bulk parameters
2. Uses electronic_properties builder for DOS/bands parameters
3. Relaxes bulk structure
4. Computes DOS and band structure for relaxed bulk using vasp.v2.bands

The workflow automatically:
- Performs SCF calculation with charge density output
- Computes band structure along high-symmetry paths
- Computes density of states with tetrahedron method

Usage:
    source ~/envs/psteros/bin/activate && python bulk_dos_bands_ag2o.py
"""

import sys
import os
# Add PS-TEROS root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from aiida import load_profile
from teros.core.workgraph import build_core_workgraph
from teros.core.builders.default_ag2o_builders import get_ag2o_defaults
from teros.core.builders import get_electronic_properties_defaults


def main():
    """Main function to run the Ag2O bulk electronic properties workflow."""

    # Load AiiDA profile
    print("Loading AiiDA profile...")
    load_profile(profile='psteros')

    # Define directories
    structures_dir = "/home/thiagotd/git/PS-TEROS/examples/structures"

    # ===== GET DEFAULT BUILDERS FOR Ag2O =====
    print(f"\n{'='*80}")
    print("GETTING DEFAULT PARAMETERS")
    print(f"{'='*80}")

    bulk_defaults = get_ag2o_defaults(
        structures_dir=structures_dir,
        code_label="VASP-VTST-6.4.3@bohr",
        potential_family="PBE",
    )

    print("✓ Bulk parameters loaded from get_ag2o_defaults()")

    # ===== GET ELECTRONIC PROPERTIES DEFAULTS =====
    ep_defaults = get_electronic_properties_defaults(
        energy_cutoff=bulk_defaults['bulk_parameters']['ENCUT'],  # Match bulk ENCUT
        electronic_convergence=1e-5,
        ncore=4,
        ispin=2,  # Spin-polarized for Ag2O
        lasph=True,
        lreal="Auto",
        kpoints_mesh_density=0.3,  # SCF k-mesh density
        band_kpoints_distance=0.2,  # Band path density
        dos_kpoints_distance=0.2,  # DOS k-mesh density
        line_density=0.2,  # Points along high-symmetry lines
        nedos=2000,  # DOS grid points
        sigma_bands=0.01,  # Smearing for bands (eV)
        symprec=1e-4,  # Symmetry precision
        band_mode="seekpath-aiida",  # Use seekpath for band paths
    )

    print("✓ Electronic properties parameters loaded from get_electronic_properties_defaults()")

    print(f"\n{'='*80}")
    print("ELECTRONIC PROPERTIES CONFIGURATION")
    print(f"{'='*80}")
    print(f"Band mode: {ep_defaults['band_settings']['band_mode']}")
    print(f"Band k-points distance: {ep_defaults['band_settings']['band_kpoints_distance']}")
    print(f"DOS k-points distance: {ep_defaults['band_settings']['dos_kpoints_distance']}")
    print(f"Line density: {ep_defaults['band_settings']['line_density']}")
    print(f"NEDOS (DOS grid points): {ep_defaults['dos']['NEDOS']}")
    print(f"SCF k-mesh density: {ep_defaults['scf_kpoints_distance']}")
    print(f"\nKey SCF parameters (for bands/DOS restart):")
    print(f"  LWAVE: {ep_defaults['scf']['LWAVE']} (must be True)")
    print(f"  LCHARG: {ep_defaults['scf']['LCHARG']} (must be True)")
    print(f"\nBands calculation:")
    print(f"  ISMEAR: {ep_defaults['bands']['ISMEAR']} (Gaussian smearing)")
    print(f"  SIGMA: {ep_defaults['bands']['SIGMA']} eV")
    print(f"\nDOS calculation:")
    print(f"  ISMEAR: {ep_defaults['dos']['ISMEAR']} (Tetrahedron method)")

    # ===== OPTIONAL: OVERRIDE PARAMETERS =====
    # Example: Change NEDOS for finer DOS grid
    # ep_defaults['dos']['NEDOS'] = 3000

    # Example: Use denser k-points for DOS
    # ep_defaults['band_settings']['dos_kpoints_distance'] = 0.15

    # ===== CREATE WORKGRAPH =====
    print(f"\n{'='*80}")
    print("Creating WorkGraph with electronic properties calculation...")
    print(f"{'='*80}")

    wg = build_core_workgraph(
        structures_dir=structures_dir,
        bulk_name="ag2o.cif",
        metal_name="Ag.cif",
        oxygen_name="O2.cif",
        code_label="VASP-VTST-6.4.3@bohr",
        potential_family="PBE",
        bulk_potential_mapping=bulk_defaults['bulk_potential_mapping'],
        metal_potential_mapping=bulk_defaults['metal_potential_mapping'],
        oxygen_potential_mapping=bulk_defaults['oxygen_potential_mapping'],
        kpoints_spacing=bulk_defaults['kpoints_spacing'],
        bulk_parameters=bulk_defaults['bulk_parameters'],
        bulk_options=bulk_defaults['bulk_options'],
        metal_parameters=bulk_defaults['metal_parameters'],
        metal_options=bulk_defaults['metal_options'],
        oxygen_parameters=bulk_defaults['oxygen_parameters'],
        oxygen_options=bulk_defaults['oxygen_options'],
        clean_workdir=False,
        # Disable slab-related calculations (we only want bulk electronic properties)
        compute_cleavage=False,
        # Electronic properties flags
        compute_electronic_properties_bulk=True,  # Enable DOS and bands
        bands_parameters=ep_defaults,  # Pass all electronic properties params
        band_settings=ep_defaults['band_settings'],  # Pass band workflow settings
        bands_options=bulk_defaults['bulk_options'],  # Use same resources as bulk
        name="Ag2O_Bulk_Electronic_Properties",
    )

    print("✓ WorkGraph created with electronic properties enabled")

    # Optional: Export to HTML
    try:
        html_file = "ag2o_bulk_electronic_properties.html"
        wg.to_html(html_file)
        print(f"\n✓ WorkGraph visualization saved to: {html_file}")
    except Exception as e:
        print(f"\n✗ Could not generate HTML visualization: {e}")

    # ===== SUBMIT WORKGRAPH =====
    print(f"\n{'='*80}")
    print("Submitting WorkGraph...")
    print(f"{'='*80}")
    wg.submit(wait=False)

    print(f"\n✓ WorkGraph submitted successfully!")
    print(f"\nWorkGraph PK: {wg.pk}")

    print(f"\n{'-'*80}")
    print(f"MONITORING")
    print(f"{'-'*80}")
    print(f"  Check status:     verdi process show {wg.pk}")
    print(f"  Check report:     verdi process report {wg.pk}")
    print(f"  List processes:   verdi process list -a -p1")
    print(f"\n  Wait ~15-30 seconds for workflow to initialize, then:")
    print(f"  verdi process show {wg.pk}")

    print(f"\n{'-'*80}")
    print(f"EXPECTED OUTPUTS")
    print(f"{'-'*80}")
    print(f"  Bulk relaxation:")
    print(f"    - bulk_structure (StructureData) - Relaxed Ag2O structure")
    print(f"    - bulk_energy (Float) - Total energy of relaxed bulk")
    print(f"\n  Reference energies:")
    print(f"    - metal_energy (Float) - Energy of Ag reference")
    print(f"    - oxygen_energy (Float) - Energy of O2 molecule")
    print(f"\n  Formation enthalpy:")
    print(f"    - formation_enthalpy (Dict) - ΔH_f for Ag2O")
    print(f"\n  Electronic properties (NEW!):")
    print(f"    - bulk_bands (BandsData) - Band structure along high-symmetry paths")
    print(f"    - bulk_dos (BandsData) - Density of states")
    print(f"    - bulk_primitive_structure (StructureData) - Primitive cell used for bands")
    print(f"    - bulk_seekpath_parameters (Dict) - Seekpath symmetry analysis")

    print(f"\n{'-'*80}")
    print(f"ACCESSING RESULTS (after completion)")
    print(f"{'-'*80}")
    print(f"  In Python:")
    print(f"    from aiida import load_node")
    print(f"    wg = load_node({wg.pk})")
    print(f"    ")
    print(f"    # Get relaxed bulk structure")
    print(f"    relaxed_bulk = wg.outputs.bulk_structure")
    print(f"    print(f'Relaxed bulk energy: {{wg.outputs.bulk_energy.value}} eV')")
    print(f"    ")
    print(f"    # Get band structure")
    print(f"    bands = wg.outputs.bulk_bands")
    print(f"    print(f'Band structure available: {{bands}}')")
    print(f"    print(f'Number of bands: {{len(bands.get_bands())}}')")
    print(f"    ")
    print(f"    # Get DOS")
    print(f"    dos = wg.outputs.bulk_dos")
    print(f"    dos_dict = dos.get_dict()")
    print(f"    print(f'DOS arrays: {{list(dos_dict.keys())}}')")
    print(f"    energies = dos_dict['energy']  # Energy grid")
    print(f"    total_dos = dos_dict['dos']  # Total DOS")
    print(f"    ")
    print(f"    # Export relaxed bulk to file")
    print(f"    atoms = relaxed_bulk.get_ase()")
    print(f"    atoms.write('ag2o_bulk_relaxed.cif')")
    print(f"    ")
    print(f"    # Plot band structure (requires matplotlib)")
    print(f"    # from aiida.tools.visualization import plot_bands")
    print(f"    # plot_bands(bands)")

    print(f"\n{'-'*80}")
    print(f"WORKFLOW DETAILS")
    print(f"{'-'*80}")
    print(f"  The vasp.v2.bands workchain internally performs:")
    print(f"  1. SCF calculation with LWAVE=True, LCHARG=True")
    print(f"  2. Band structure calculation (non-SCF, ICHARG=11)")
    print(f"  3. DOS calculation with tetrahedron method (ISMEAR=-5)")
    print(f"\n  High-symmetry k-path is determined automatically using:")
    print(f"  - Band mode: {ep_defaults['band_settings']['band_mode']}")
    print(f"  - Symmetry precision: {ep_defaults['band_settings']['symprec']}")

    print(f"\n{'-'*80}")
    print(f"IMPORTANT NOTES")
    print(f"{'-'*80}")
    print(f"  • This workflow computes DOS and bands for BULK Ag2O only")
    print(f"  • For slab DOS/bands, you'll need to add similar functionality for slabs")
    print(f"  • SCF calculation uses dense k-mesh ({ep_defaults['scf_kpoints_distance']} spacing)")
    print(f"  • DOS uses tetrahedron method (ISMEAR=-5) for accurate integration")
    print(f"  • Bands use Gaussian smearing (ISMEAR=0) for smooth curves")
    print(f"  • Expected runtime: ~5-15 minutes depending on system size and resources")

    print(f"\n{'='*80}\n")

    return wg


if __name__ == "__main__":
    """
    Run the bulk electronic properties calculation workflow for Ag2O.

    Before running:
    1. Make sure AiiDA profile 'psteros' is set as default:
       verdi profile set-default psteros

    2. Check AiiDA status:
       verdi status

    3. Start daemon if not running:
       verdi daemon start

    4. Ensure structure files exist:
       ls /home/thiagotd/git/PS-TEROS/examples/structures/
       Should contain: ag2o.cif, Ag.cif, O2.cif

    5. Clear Python cache if you made changes:
       find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete

    6. Restart daemon after code changes:
       verdi daemon restart

    7. Run this script:
       source ~/envs/psteros/bin/activate && python bulk_dos_bands_ag2o.py
    """
    try:
        wg = main()
    except Exception as e:
        print(f"\n{'='*80}")
        print(f"ERROR")
        print(f"{'='*80}")
        print(f"{e}")
        print(f"\nFull traceback:")
        import traceback
        traceback.print_exc()
        import sys
        sys.exit(1)

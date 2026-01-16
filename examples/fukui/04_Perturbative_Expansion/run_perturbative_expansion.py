#!/usr/bin/env python
"""
Perturbative Expansion Model - SnO2 (110) Surface Example

This example demonstrates the full Phase 1 + Phase 2 + Phase 4 Fukui workflow:

Phase 1 (Interpolation):
    1. Runs 4 parallel VASP calculations at different charge states
    2. Collects CHGCAR files and computes Fukui function via FukuiGrid

Phase 2 (Electrodes):
    3. Runs DFPT calculation on bulk SnO2 to compute dielectric constant
    4. Applies electrodes method to compute Fukui potential (LOCPOT_FUKUI.vasp)

Phase 4 (Perturbative Expansion):
    5. Retrieves LOCPOT from neutral slab calculation (electrostatic potential)
    6. Computes interaction energy map using c-DFT perturbative expansion:
       DeltaU(r) = q * Phi(r) - q * DeltaN * vf(r)
    7. Outputs MODELPOT_LOCPOT.vasp showing favorable/unfavorable adsorption sites

The output MODELPOT_LOCPOT.vasp can be visualized in VESTA:
- Negative values (blue) = favorable adsorption sites
- Positive values (red) = unfavorable adsorption sites

Usage:
    source ~/envs/aiida/bin/activate
    python run_perturbative_expansion.py

Requirements:
    - AiiDA profile configured
    - VASP code registered
    - POTCAR files available
    - FukuiGrid cloned to teros/external/FukuiGrid/
"""

import os
import sys
from pathlib import Path
from aiida import orm, load_profile
from teros.core.fukui import (
    build_fukui_workgraph,
    get_fukui_results,
    print_fukui_summary,
    calculate_nelect,
    print_nelect_breakdown,
)


def create_sno2_bulk_structure() -> orm.StructureData:
    """
    Create bulk SnO2 rutile structure (mp-856 equivalent).

    Returns:
        AiiDA StructureData for bulk SnO2 rutile
    """
    from pymatgen.core import Structure, Lattice

    # SnO2 rutile structure (P4_2/mnm, space group 136)
    # Lattice parameters from Materials Project mp-856
    a = 4.7374
    c = 3.1864

    lattice = Lattice.tetragonal(a, c)

    # Atomic positions in fractional coordinates
    species = ["Sn", "Sn", "O", "O", "O", "O"]
    coords = [
        [0.0, 0.0, 0.0],           # Sn
        [0.5, 0.5, 0.5],           # Sn
        [0.3056, 0.3056, 0.0],     # O
        [0.6944, 0.6944, 0.0],     # O
        [0.1944, 0.8056, 0.5],     # O
        [0.8056, 0.1944, 0.5],     # O
    ]

    structure_pmg = Structure(lattice, species, coords)
    bulk_structure = orm.StructureData(pymatgen=structure_pmg)

    print(f"   Created SnO2 rutile (mp-856 equivalent)")
    print(f"   Lattice: a = {a:.4f} A, c = {c:.4f} A")
    print(f"   Space group: P4_2/mnm (136)")
    print(f"   Number of atoms: {len(bulk_structure.sites)}")

    return bulk_structure


def fetch_bulk_structure(material_id: str = "mp-856") -> orm.StructureData:
    """
    Fetch bulk structure from Materials Project, or create it if API unavailable.

    Args:
        material_id: Materials Project ID (default: mp-856 for SnO2 rutile)

    Returns:
        AiiDA StructureData for the bulk material
    """
    try:
        from mp_api.client import MPRester

        print(f"   Fetching {material_id} from Materials Project...")

        # Check for API key
        api_key = os.environ.get('MP_API_KEY')
        if api_key:
            mpr = MPRester(api_key)
        else:
            # Will use key from ~/.config/.pmgrc.yaml if available
            mpr = MPRester()

        with mpr:
            # Get structure from Materials Project
            structure_mp = mpr.get_structure_by_material_id(material_id)

        # Convert pymatgen Structure to AiiDA StructureData
        bulk_structure = orm.StructureData(pymatgen=structure_mp)

        print(f"   Composition: {bulk_structure.get_composition()}")
        print(f"   Space group: {structure_mp.get_space_group_info()[0]}")
        print(f"   Number of atoms: {len(bulk_structure.sites)}")

        return bulk_structure

    except ImportError:
        print("   mp-api not installed, creating structure from scratch...")
        return create_sno2_bulk_structure()
    except Exception as e:
        print(f"   Materials Project API error: {e}")
        print("   Creating structure from scratch instead...")
        return create_sno2_bulk_structure()


def main():
    """Run perturbative expansion (Phase 4) calculation example."""

    print("\n" + "=" * 70)
    print("PERTURBATIVE EXPANSION MODEL - SnO2 (110) Surface")
    print("Phase 1 (Interpolation) + Phase 2 (Electrodes) + Phase 4 (Perturbative)")
    print("=" * 70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='presto')
    print("   Profile loaded")

    # Check daemon status
    from aiida.engine.daemon.client import get_daemon_client
    client = get_daemon_client()
    if not client.is_daemon_running:
        print("\n   WARNING: AiiDA daemon is not running!")
        print("   Start with: verdi daemon start")
        return 1
    print("   Daemon is running")

    # Load slab structure (from Phase 1 example)
    print("\n2. Loading SnO2 (110) slab structure...")
    structure_file = Path(__file__).parent.parent / "01_Interpolation" / "sno2_110_paper.vasp"
    if structure_file.exists():
        from ase.io import read
        atoms = read(str(structure_file))
        slab_structure = orm.StructureData(ase=atoms)
        print(f"   Loaded from file: {structure_file}")
    else:
        print(f"   ERROR: Slab structure file not found: {structure_file}")
        print("   Please ensure the Phase 1 example directory exists with sno2_110_paper.vasp")
        return 1

    print(f"   Slab composition: {slab_structure.get_composition()}")
    print(f"   Number of atoms: {len(slab_structure.sites)}")

    # Fetch bulk structure from Materials Project
    print("\n3. Fetching bulk SnO2 structure from Materials Project...")
    try:
        bulk_structure = fetch_bulk_structure("mp-856")
    except Exception as e:
        print(f"   ERROR: Could not fetch from Materials Project: {e}")
        print("   Make sure you have set MP_API_KEY or have ~/.config/.pmgrc.yaml configured")
        return 1

    # Define POTCAR family and mapping
    potential_family = 'PBE'
    potential_mapping = {
        'Sn': 'Sn_d',  # Sn_d has 14 valence electrons
        'O': 'O',      # O has 6 valence electrons
    }

    # Calculate NELECT for slab
    print("\n4. Calculating NELECT for slab from POTCAR files...")
    nelect_neutral = calculate_nelect(
        structure=slab_structure,
        potential_family=potential_family,
        potential_mapping=potential_mapping,
    )
    print_nelect_breakdown(slab_structure, potential_family, potential_mapping)

    # Define perturbative expansion parameters
    print("\n5. Defining perturbative expansion parameters...")

    # Probe charge and electron transfer for a typical cationic adsorbate
    # Example: Li+ approaching the surface
    # q = +0.3 (partial positive charge)
    # DeltaN = -0.3 (surface donates electrons to adsorbate)
    probe_charge = 0.3
    electron_transfer = -0.3

    print(f"   Probe charge (q): {probe_charge} |e|")
    print(f"   Electron transfer (DeltaN): {electron_transfer}")
    print(f"   Interpretation: Partially charged cation approaching surface")
    print(f"   Formula: DeltaU(r) = q*Phi(r) - q*DeltaN*vf(r)")

    # Define builder inputs
    print("\n6. Defining VASP calculation inputs...")

    # Code label for obelix cluster
    code_label = 'VASP-6.5.1-idefix@obelix'

    # Common VASP settings for charge density calculations
    builder_inputs = {
        'parameters': {
            'incar': {
                # Precision and cutoff
                'prec': 'Accurate',
                'encut': 500,

                # SCF convergence - robust settings for fractional charges
                'algo': 'All',
                'nelm': 300,
                'ediff': 1e-6,

                # Charge density initialization
                'icharg': 2,

                # Mixing parameters for better convergence
                'amix': 0.2,
                'bmix': 0.0001,
                'amix_mag': 0.8,
                'bmix_mag': 0.0001,
                'lmaxmix': 4,

                # Smearing
                'ismear': 0,
                'sigma': 0.05,

                # Real-space projection
                'lreal': 'Auto',

                # Spin settings
                'ispin': 2,

                # Parallelization
                'ncore': 2,
                'kpar': 1,
            }
        },
        'options': {
            'resources': {
                'num_machines': 1,
                'num_mpiprocs_per_machine': 4,
            },
            'custom_scheduler_commands': '''#PBS -l cput=90000:00:00
#PBS -l nodes=1:ppn=88:skylake
#PBS -j oe
#PBS -N Perturbative_SnO2''',
        },
        'kpoints_spacing': 0.03,
        'potential_family': potential_family,
        'potential_mapping': potential_mapping,
        'clean_workdir': False,
    }

    # DFPT-specific settings
    dfpt_builder_inputs = {
        'parameters': {
            'incar': {
                'prec': 'Accurate',
                'encut': 500,
                'ediff': 1e-8,      # Tighter convergence for DFPT
                'ismear': 0,
                'sigma': 0.05,
                'ispin': 1,         # Non-spin-polarized for bulk DFPT
                'ncore': 2,
                'kpar': 1,
            }
        },
        'options': {
            'resources': {
                'num_machines': 1,
                'num_mpiprocs_per_machine': 4,
            },
            'custom_scheduler_commands': '''#PBS -l cput=90000:00:00
#PBS -l nodes=1:ppn=88:skylake
#PBS -j oe
#PBS -N DFPT_SnO2_bulk''',
        },
        'kpoints_spacing': 0.02,    # Finer k-mesh for DFPT
        'potential_family': potential_family,
        'potential_mapping': potential_mapping,
        'clean_workdir': False,
    }

    delta_n_values = [0.0, 0.05, 0.10, 0.15]

    print(f"   Code: {code_label}")
    print(f"   ENCUT: {builder_inputs['parameters']['incar']['encut']} eV")
    print(f"   Slab k-points spacing: {builder_inputs['kpoints_spacing']} A^-1")
    print(f"   DFPT k-points spacing: {dfpt_builder_inputs['kpoints_spacing']} A^-1")
    print(f"   Delta N values: {delta_n_values}")

    # Build WorkGraph with all phases enabled
    print("\n7. Building Fukui WorkGraph (Phase 1 + Phase 2 + Phase 4)...")
    wg = build_fukui_workgraph(
        structure=slab_structure,
        nelect_neutral=nelect_neutral,
        delta_n_values=delta_n_values,
        code_label=code_label,
        builder_inputs=builder_inputs,
        relax_first=False,
        fukui_type='plus',
        # Phase 1: Interpolation
        compute_fukui=True,
        # Phase 2: Electrodes
        compute_fukui_potential=True,
        bulk_structure=bulk_structure,
        dfpt_builder_inputs=dfpt_builder_inputs,
        # Phase 4: Perturbative Expansion (NEW!)
        compute_perturbative_expansion=True,
        probe_charge=probe_charge,
        electron_transfer=electron_transfer,
        max_concurrent_jobs=5,
        name='Perturbative_Expansion_SnO2_110',
    )

    print(f"   WorkGraph created: {wg.name}")
    print(f"   compute_fukui=True: Phase 1 (FukuiGrid interpolation)")
    print(f"   compute_fukui_potential=True: Phase 2 (Electrodes method)")
    print(f"   compute_perturbative_expansion=True: Phase 4 (Perturbative expansion)")

    # Print expected calculations
    print("\n8. Expected VASP calculations:")
    print("-" * 50)
    print("   Phase 1 (Slab charge states):")
    for delta_n in delta_n_values:
        nelect = nelect_neutral - delta_n
        locpot_note = " + LOCPOT retrieval" if delta_n == 0.0 else ""
        print(f"      delta_N = {delta_n:.2f} -> NELECT = {nelect:.2f}{locpot_note}")
    print("   Phase 2 (DFPT for dielectric constant):")
    print(f"      Bulk SnO2 (mp-856) with LEPSILON=.TRUE.")
    print("   Phase 4 (Post-processing):")
    print(f"      Perturbative expansion with q={probe_charge}, DeltaN={electron_transfer}")

    # Submit
    print("\n9. Submitting WorkGraph...")
    wg.submit(wait=False)
    print(f"   Submitted!")
    print(f"   PK: {wg.pk}")

    # Save PK to file
    pk_file = Path(__file__).parent / 'last_pk.txt'
    with open(pk_file, 'w') as f:
        f.write(f"{wg.pk}\n")
    print(f"   PK saved to: {pk_file}")

    print("\n" + "=" * 70)
    print("NEXT STEPS")
    print("=" * 70)
    print(f"""
Monitor progress:
    verdi process show {wg.pk}
    verdi process report {wg.pk}

After completion, extract results:
    >>> from teros.core.fukui import get_fukui_results, print_fukui_summary
    >>> results = get_fukui_results({wg.pk})
    >>> print_fukui_summary({wg.pk})

Phase 1 Outputs:
    - results['chgcar_folder']: FolderData with CHGCAR files
    - results['fukui_chgcar']: SinglefileData with CHGCAR_FUKUI.vasp

Phase 2 Outputs:
    - results['dielectric_constant']: Float with computed epsilon
    - results['fukui_potential']: SinglefileData with LOCPOT_FUKUI.vasp

Phase 4 Outputs (NEW!):
    - results['locpot_neutral']: SinglefileData with LOCPOT (electrostatic potential)
    - results['modelpot']: SinglefileData with MODELPOT_LOCPOT.vasp (interaction energy)

Export model potential for visualization in VESTA:
    >>> modelpot_file = results['modelpot']
    >>> content = modelpot_file.get_content()
    >>> with open('MODELPOT_LOCPOT.vasp', 'wb' if isinstance(content, bytes) else 'w') as f:
    ...     f.write(content)

Interpretation of MODELPOT_LOCPOT.vasp:
    - Negative values (blue in VESTA) = favorable adsorption sites
    - Positive values (red in VESTA) = unfavorable adsorption sites

Post-hoc analysis with different q/DeltaN values:
    >>> from teros.core.fukui import run_perturbative_expansion_calcfunc
    >>> from aiida import orm
    >>> locpot = results['locpot_neutral']
    >>> fukui_pot = results['fukui_potential']
    >>> # Try different probe charge / electron transfer
    >>> new_result = run_perturbative_expansion_calcfunc(
    ...     locpot_neutral=locpot,
    ...     fukui_potential=fukui_pot,
    ...     probe_charge=orm.Float(0.5),      # Different charge
    ...     electron_transfer=orm.Float(-0.5), # Different transfer
    ... )
""")
    print("=" * 70)

    return 0


if __name__ == '__main__':
    sys.exit(main())

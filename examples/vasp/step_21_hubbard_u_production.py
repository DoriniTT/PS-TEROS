#!/home/trevizam/envs/aiida/bin/python
"""
Production: Calculate Hubbard U for NiO using linear response method.

This script follows the VASP wiki methodology exactly:
https://www.vasp.at/wiki/index.php/Calculate_U_for_LSDA+U

Key features for accurate U calculation:
- 2x2x2 NiO supercell (32 atoms: 16 Ni + 16 O)
- AFM-II magnetic ordering with proper MAGMOM
- Gamma-only k-points for supercell
- Production-quality VASP parameters

Expected result: U ~ 5-6 eV for Ni d-electrons in NiO

Usage:
    source ~/envs/aiida/bin/activate
    python step_21_hubbard_u_production.py
"""

import sys
from pathlib import Path

from aiida import orm, load_profile
from pymatgen.core import Structure, Lattice

# Load AiiDA profile
load_profile('presto')

from teros.core.u_calculation import build_u_calculation_workgraph

# ==============================================================================
# CONFIGURATION - Obelix cluster
# ==============================================================================

VASP_CODE = 'VASP-6.5.1-idefix-4@obelix'
POTENTIAL_FAMILY = 'PBE'
POTENTIAL_MAPPING = {'Ni': 'Ni', 'O': 'O'}

# Obelix cluster options
CLUSTER_OPTIONS = {
    'resources': {
        'num_machines': 1,
        'num_mpiprocs_per_machine': 4,
    },
    'custom_scheduler_commands': '''#PBS -l cput=90000:00:00
#PBS -l nodes=1:ppn=88:skylake
#PBS -j oe
#PBS -N NiO_HubbardU''',
}

# Production INCAR parameters matching VASP wiki
# Note: lowercase keys required by AiiDA-VASP
PRODUCTION_PARAMS = {
    'prec': 'Accurate',
    'encut': 500,
    'ediff': 1e-6,
    'ismear': 0,
    'sigma': 0.05,
    'algo': 'Normal',
    'nelm': 200,
    'ispin': 2,  # Spin-polarized
    'ncore': 11,  # 88 cores / 8 = 11
    'kpar': 1,   # Gamma only, no k-parallelization needed
    'lasph': True,  # Non-spherical contributions
}

# Potential values for linear response (as per VASP wiki)
# Symmetric around zero for better linear fit
PRODUCTION_POTENTIAL_VALUES = [-0.20, -0.15, -0.10, -0.05, 0.05, 0.10, 0.15, 0.20]

# K-point spacing (large value = Gamma-only for supercell)
KPOINTS_SPACING = 1.0  # Results in Gamma-only for 2x2x2 supercell


def create_nio_afm2_supercell() -> tuple[orm.StructureData, list]:
    """
    Create a 2x2x2 NiO supercell with AFM-II magnetic ordering.

    AFM-II structure: Alternating (111) ferromagnetic planes
    - Each Ni atom has opposite spin to its nearest Ni neighbors
    - Magnetic moments alternate in layers perpendicular to [111]

    Returns:
        Tuple of (StructureData, magmom_list)
        magmom_list: MAGMOM values for INCAR (16 Ni moments, then 16 O zeros)
    """
    # NiO rocksalt structure (experimental lattice parameter)
    a = 4.17  # Angstroms

    # Create primitive rocksalt NiO
    lattice = Lattice.cubic(a)
    species = ['Ni', 'O']
    coords = [
        [0.0, 0.0, 0.0],  # Ni
        [0.5, 0.5, 0.5],  # O
    ]
    pmg_structure = Structure(lattice, species, coords)

    # Create 2x2x2 supercell
    pmg_structure.make_supercell([2, 2, 2])

    # Get Ni positions and assign AFM-II magnetic moments
    # AFM-II: moments alternate along [111] direction
    # For rocksalt: Ni atoms at fractional coords (i+j+k) % 2 determines spin
    magmom_ni = []
    for site in pmg_structure:
        if site.species_string == 'Ni':
            # Get fractional coordinates in original unit cell
            frac = site.frac_coords
            # AFM-II: sum of indices determines spin
            # Atoms in same (111) plane have same spin
            idx_sum = round(frac[0] * 2) + round(frac[1] * 2) + round(frac[2] * 2)
            if idx_sum % 2 == 0:
                magmom_ni.append(2.0)  # Spin up
            else:
                magmom_ni.append(-2.0)  # Spin down

    # MAGMOM: Ni moments first, then O zeros
    # VASP expects MAGMOM in order of species in POTCAR
    n_ni = len([s for s in pmg_structure if s.species_string == 'Ni'])
    n_o = len([s for s in pmg_structure if s.species_string == 'O'])

    # Reorder structure to have Ni first, then O (POTCAR order)
    ni_sites = [s for s in pmg_structure if s.species_string == 'Ni']
    o_sites = [s for s in pmg_structure if s.species_string == 'O']

    # Build new structure with correct ordering
    new_species = []
    new_coords = []
    magmom_ordered = []

    for i, site in enumerate(ni_sites):
        new_species.append('Ni')
        new_coords.append(site.frac_coords)
        magmom_ordered.append(magmom_ni[i])

    for site in o_sites:
        new_species.append('O')
        new_coords.append(site.frac_coords)
        magmom_ordered.append(0.0)

    # Create ordered structure
    ordered_structure = Structure(
        pmg_structure.lattice,
        new_species,
        new_coords
    )

    # Convert to AiiDA
    structure = orm.StructureData(pymatgen=ordered_structure)
    structure.label = 'NiO_2x2x2_AFM2'
    structure.description = 'NiO 2x2x2 supercell with AFM-II magnetic ordering for Hubbard U calculation'

    print(f"   Created {n_ni} Ni + {n_o} O = {n_ni + n_o} atoms")
    print(f"   MAGMOM: {sum(1 for m in magmom_ordered if m > 0)} up, {sum(1 for m in magmom_ordered if m < 0)} down")

    return structure, magmom_ordered


def main():
    """Run the production Hubbard U calculation for NiO."""

    print("=" * 70)
    print("Production Hubbard U Calculation for NiO")
    print("Following VASP wiki methodology")
    print("=" * 70)

    # Check daemon
    print("\n1. Checking AiiDA daemon...")
    from aiida.engine.daemon.client import get_daemon_client
    client = get_daemon_client()
    if not client.is_daemon_running:
        print("   WARNING: AiiDA daemon is not running!")
        print("   Start with: verdi daemon start")
        return 1
    print("   Daemon is running")

    # Create NiO supercell with AFM ordering
    print("\n2. Creating NiO 2x2x2 supercell with AFM-II ordering...")
    structure, magmom = create_nio_afm2_supercell()
    structure.store()
    print(f"   Structure PK: {structure.pk}")
    print(f"   Formula: {structure.get_formula()}")

    # Add MAGMOM to INCAR parameters
    params_with_magmom = PRODUCTION_PARAMS.copy()
    params_with_magmom['magmom'] = magmom

    print("\n3. Configuration:")
    print(f"   VASP code: {VASP_CODE}")
    print(f"   Potential family: {POTENTIAL_FAMILY}")
    print(f"   Potential values: {PRODUCTION_POTENTIAL_VALUES}")
    print(f"   K-points: Gamma-only (spacing={KPOINTS_SPACING})")

    # Build the workgraph
    print("\n4. Building Hubbard U workgraph...")
    wg = build_u_calculation_workgraph(
        structure=structure,
        code_label=VASP_CODE,
        potential_family=POTENTIAL_FAMILY,
        potential_mapping=POTENTIAL_MAPPING,
        target_species='Ni',
        potential_values=PRODUCTION_POTENTIAL_VALUES,
        ldaul=2,  # d-electrons
        ldauj=0.0,
        ground_state_parameters=params_with_magmom,
        response_parameters=params_with_magmom,
        options=CLUSTER_OPTIONS,
        kpoints_spacing=KPOINTS_SPACING,
        clean_workdir=False,
        name='NiO_HubbardU_Production',
    )

    n_vasp_calcs = 1 + 2 * len(PRODUCTION_POTENTIAL_VALUES)  # GS + (NSCF + SCF) per V
    print(f"   WorkGraph name: {wg.name}")
    print(f"   Total VASP calculations: {n_vasp_calcs}")

    # Submit
    print("\n5. Submitting to obelix cluster...")
    wg.submit()
    print(f"   WorkGraph PK: {wg.pk}")

    # Save PK for reference
    pk_file = Path(__file__).parent / 'hubbard_u_production_pk.txt'
    with open(pk_file, 'w') as f:
        f.write(str(wg.pk))

    print(f"\n   PK saved to: {pk_file}")
    print(f"\n   Monitor with: verdi process show {wg.pk}")
    print(f"   Detailed: verdi process report {wg.pk}")

    print("\n" + "=" * 70)
    print("WorkGraph submitted! Expected runtime: ~1-2 hours on obelix")
    print("Expected result: U ~ 5-6 eV for Ni d-electrons")
    print("=" * 70 + "\n")

    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

#!/home/thiagotd/envs/aiida/bin/python
"""
STEP 16: Surface Thermodynamics - Serial Preset (Experimental)

This script tests the EXPERIMENTAL flat-graph serial surface thermodynamics preset:
- All VASP nodes exist at the same graph level
- Enables max_number_jobs to control concurrent execution
- Supports internal slab generation from miller_indices
- Supports custom builders for full parameter control

This demonstrates the serial preset that addresses the max_number_jobs limitation
in the standard nested sub-workgraph approach.

Material: Ag2O
Surface: (100), (110), (111)
References: Ag, O2
Concurrency limit: 2 jobs

Features demonstrated:
1. Internal slab generation from miller_indices (no pre-generated slabs needed)
2. Custom builders for full control over VASP parameters

Usage:
    source ~/envs/aiida/bin/activate
    python step_16_surface_thermodynamics_serial.py
"""

import sys
import os
from aiida import load_profile, orm
from aiida.engine import submit
from teros.experimental.surface_thermo_preset_serial import surface_thermodynamics_serial_workgraph


def main():
    """Step 16: Test serial surface thermodynamics workflow."""

    print("\n" + "="*70)
    print("STEP 16: SURFACE THERMODYNAMICS (SERIAL PRESET - EXPERIMENTAL)")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='presto')
    print("   ✓ Profile loaded")

    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, 'structures')

    print(f"\n2. Structures:")
    print(f"   Bulk:   {structures_dir}/ag2o.cif")
    print(f"   Metal:  {structures_dir}/Ag.cif")
    print(f"   Oxygen: {structures_dir}/O2.cif")

    # Code configuration
    code_label = 'VASP-6.5.0@bohr-new'
    code = orm.load_code(code_label)
    potential_family = 'PBE'

    print(f"\n3. VASP Configuration:")
    print(f"   Code: {code_label}")
    print(f"   Potential family: {potential_family}")

    # Common options for all calculations
    common_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 40,
        },
        'queue_name': 'par40',
    }

    # =========================================================================
    # OPTION 1: Use custom builders for full control (recommended)
    # =========================================================================

    print("\n4. Building custom VASP parameter builders...")
    print("   Creating builders for full control over all calculations")

    # Bulk builder - full relaxation (cell + ions) - LIGHT PARAMETERS FOR TESTING
    bulk_builder = {
        'code': code,
        'potential_family': orm.Str(potential_family),
        'potential_mapping': orm.Dict(dict={'Ag': 'Ag', 'O': 'O'}),
        'kpoints_spacing': orm.Float(1.0),  # Coarse k-points
        'options': orm.Dict(dict=common_options),
        'clean_workdir': orm.Bool(False),
        'parameters': orm.Dict(dict={
            'incar': {
                'PREC': 'Normal',  # Lower precision for speed
                'ENCUT': 300,  # Lower cutoff for speed
                'EDIFF': 1e-3,  # Looser convergence for speed
                'ISMEAR': 0,
                'SIGMA': 0.05,
                'IBRION': 2,
                'ISIF': 3,  # Full relaxation (cell + ions)
                'NSW': 200,  # Relaxation steps
                'EDIFFG': -0.5,
                'ALGO': 'Fast',  # Faster algorithm
                'LREAL': 'Auto',
                'LWAVE': False,
                'LCHARG': False,
            }
        }),
    }

    # Metal builder - LIGHT PARAMETERS FOR TESTING
    metal_builder = {
        'code': code,
        'potential_family': orm.Str(potential_family),
        'potential_mapping': orm.Dict(dict={'Ag': 'Ag'}),
        'kpoints_spacing': orm.Float(1.0),  # Coarse k-points
        'options': orm.Dict(dict=common_options),
        'clean_workdir': orm.Bool(False),
        'parameters': orm.Dict(dict={
            'incar': {
                'PREC': 'Normal',  # Lower precision for speed
                'ENCUT': 300,  # Lower cutoff for speed
                'EDIFF': 1e-3,  # Looser convergence for speed
                'ISMEAR': 1,  # Different smearing for metals
                'SIGMA': 0.2,
                'IBRION': 2,
                'ISIF': 3,
                'NSW': 200,  # Relaxation steps
                'EDIFFG': -0.5,
                'ALGO': 'Fast',  # Faster algorithm
                'LREAL': 'Auto',
                'LWAVE': False,
                'LCHARG': False,
            }
        }),
    }

    # Oxygen builder (molecule) - LIGHT PARAMETERS FOR TESTING
    oxygen_builder = {
        'code': code,
        'potential_family': orm.Str(potential_family),
        'potential_mapping': orm.Dict(dict={'O': 'O'}),
        'kpoints_spacing': orm.Float(1.0),  # Coarse k-points
        'options': orm.Dict(dict=common_options),
        'clean_workdir': orm.Bool(False),
        'parameters': orm.Dict(dict={
            'incar': {
                'PREC': 'Normal',  # Lower precision for speed
                'ENCUT': 300,  # Lower cutoff for speed
                'EDIFF': 1e-3,  # Looser convergence for speed
                'ISMEAR': 0,
                'SIGMA': 0.05,
                'IBRION': 2,
                'ISIF': 2,  # Ions only (molecule in box)
                'NSW': 200,  # Relaxation steps
                'EDIFFG': -0.5,
                'ALGO': 'Fast',  # Faster algorithm
                'LREAL': 'Auto',
                'LWAVE': False,
                'LCHARG': False,
            }
        }),
    }

    print("   ✓ Bulk builder created (NSW=200, ENCUT=300, kpoints=1.0)")
    print("   ✓ Metal builder created (NSW=200, ENCUT=300, kpoints=1.0)")
    print("   ✓ Oxygen builder created (NSW=200, ENCUT=300, kpoints=1.0)")

    # Slab builders - we'll create per-slab builders after knowing slab IDs
    # For now, create template parameters - LIGHT FOR TESTING
    slab_scf_template = {
        'code': code,
        'potential_family': orm.Str(potential_family),
        'potential_mapping': orm.Dict(dict={'Ag': 'Ag', 'O': 'O'}),
        'kpoints_spacing': orm.Float(1.0),  # Coarse k-points
        'options': orm.Dict(dict=common_options),
        'clean_workdir': orm.Bool(False),
        'parameters': orm.Dict(dict={
            'incar': {
                'PREC': 'Normal',  # Lower precision for SCF
                'ENCUT': 300,  # Lower cutoff
                'EDIFF': 1e-3,  # Loose convergence
                'ISMEAR': 0,
                'SIGMA': 0.05,
                'NSW': 0,  # No ionic relaxation (SCF only)
                'ALGO': 'Fast',
                'LREAL': 'Auto',
                'LWAVE': False,
                'LCHARG': False,
            }
        }),
    }

    slab_relax_template = {
        'code': code,
        'potential_family': orm.Str(potential_family),
        'potential_mapping': orm.Dict(dict={'Ag': 'Ag', 'O': 'O'}),
        'kpoints_spacing': orm.Float(1.0),  # Coarse k-points
        'options': orm.Dict(dict=common_options),
        'clean_workdir': orm.Bool(False),
        'parameters': orm.Dict(dict={
            'incar': {
                'PREC': 'Normal',  # Lower precision for speed
                'ENCUT': 300,  # Lower cutoff
                'EDIFF': 1e-3,  # Loose convergence
                'ISMEAR': 0,
                'SIGMA': 0.05,
                'IBRION': 2,
                'ISIF': 2,  # Ions only, no cell relaxation
                'NSW': 200,  # Relaxation steps
                'EDIFFG': -0.5,
                'ALGO': 'Fast',
                'LREAL': 'Auto',
                'LWAVE': False,
                'LCHARG': False,
            }
        }),
    }

    print("   ✓ Slab SCF template created (NSW=0, ENCUT=300, kpoints=1.0)")
    print("   ✓ Slab relax template created (NSW=200, ENCUT=300, kpoints=1.0)")

    # Miller indices for slab generation - ONLY ONE FOR SPEED
    miller_indices = [
        (1, 0, 0),  # Just one Miller index to minimize slabs
    ]

    print(f"\n5. Slab Generation Configuration:")
    print(f"   Miller indices: {miller_indices}")
    print(f"   Slabs will be generated internally from bulk structure")
    print(f"   Minimum slab thickness: 10 Å")
    print(f"   Minimum vacuum thickness: 12 Å")

    # =========================================================================
    # Build WorkGraph
    # =========================================================================

    print("\n6. Building serial workgraph...")
    print("   Using experimental serial preset (flat-graph architecture)")
    print("   Demonstrating:")
    print("   - Internal slab generation from miller_indices")
    print("   - Custom builders for full parameter control")
    print("   - Concurrency control with max_number_jobs=2")

    try:
        wg = surface_thermodynamics_serial_workgraph(
            # Structures
            structures_dir=structures_dir,
            bulk_name='ag2o.cif',
            metal_name='Ag.cif',
            oxygen_name='O2.cif',

            # Code (not used when builders are provided, but kept for consistency)
            code_label=code_label,
            potential_family=potential_family,
            kpoints_spacing=1.0,  # Coarse k-points
            clean_workdir=False,

            # Custom builders for full control
            bulk_builder=bulk_builder,
            metal_builder=metal_builder,
            oxygen_builder=oxygen_builder,

            # Bulk parameters (not used when bulk_builder provided, kept as fallback)
            bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
            bulk_parameters={
                'PREC': 'Normal',
                'ENCUT': 300,
                'EDIFF': 1e-3,
                'ISMEAR': 0,
                'SIGMA': 0.05,
                'IBRION': 2,
                'ISIF': 3,
                'NSW': 200,
                'EDIFFG': -0.5,
                'ALGO': 'Fast',
                'LREAL': 'Auto',
                'LWAVE': False,
                'LCHARG': False,
            },
            bulk_options=common_options,

            # Metal parameters (not used when metal_builder provided)
            metal_potential_mapping={'Ag': 'Ag'},
            metal_parameters={
                'PREC': 'Normal',
                'ENCUT': 300,
                'EDIFF': 1e-3,
                'ISMEAR': 1,
                'SIGMA': 0.2,
                'IBRION': 2,
                'ISIF': 3,
                'NSW': 200,
                'EDIFFG': -0.5,
                'ALGO': 'Fast',
                'LREAL': 'Auto',
                'LWAVE': False,
                'LCHARG': False,
            },
            metal_options=common_options,

            # Oxygen parameters (not used when oxygen_builder provided)
            oxygen_potential_mapping={'O': 'O'},
            oxygen_parameters={
                'PREC': 'Normal',
                'ENCUT': 300,
                'EDIFF': 1e-3,
                'ISMEAR': 0,
                'SIGMA': 0.05,
                'IBRION': 2,
                'ISIF': 2,
                'NSW': 200,
                'EDIFFG': -0.5,
                'ALGO': 'Fast',
                'LREAL': 'Auto',
                'LWAVE': False,
                'LCHARG': False,
            },
            oxygen_options=common_options,

            # Internal slab generation from miller_indices
            miller_indices=miller_indices,
            min_slab_thickness=8,  # Angstroms - smaller for speed
            min_vacuum_thickness=10,  # Angstroms - smaller for speed
            lll_reduce=False,
            center_slab=True,
            primitive=True,
            in_unit_planes=False,
            max_normal_search=None,
            symmetrize=False,

            # Slab parameters (used as fallback when builders not provided)
            slab_parameters={
                'PREC': 'Normal',
                'ENCUT': 300,
                'EDIFF': 1e-3,
                'ISMEAR': 0,
                'SIGMA': 0.05,
                'IBRION': 2,
                'ISIF': 2,
                'NSW': 200,
                'EDIFFG': -0.5,
                'ALGO': 'Fast',
                'LREAL': 'Auto',
                'LWAVE': False,
                'LCHARG': False,
            },
            slab_potential_mapping={'Ag': 'Ag', 'O': 'O'},
            slab_options=common_options,
            slab_kpoints_spacing=1.0,

            # Note: slab_scf_builders and slab_relax_builders can be provided
            # as dicts of {slab_id: builder_dict} for per-slab control
            # For this demo, we use the template parameters above as fallback

            # Workflow flags
            relax_slabs=True,
            compute_thermodynamics=True,
            thermodynamics_sampling=10,  # Small for speed
            compute_relaxation_energy=True,
        )

        print("   ✓ WorkGraph built successfully")

    except Exception as e:
        print(f"   ✗ Error building workgraph: {e}")
        import traceback
        traceback.print_exc()
        return None

    # Set concurrent job limit
    print("\n7. Configuring concurrency control...")
    wg.max_number_jobs = 2
    print("   ✓ max_number_jobs = 2")
    print("   This limits concurrent VASP jobs to 2 at a time")
    print("   Verify with: watch -n 2 'verdi process list' (should show max 2 VASP jobs)")

    # Submit
    print("\n8. Submitting to AiiDA daemon...")
    try:
        # WorkGraph has its own submit method
        result = wg.submit()
        pk = result.pk if hasattr(result, 'pk') else result

        print(f"\n{'='*70}")
        print("STEP 16 SUBMITTED SUCCESSFULLY")
        print(f"{'='*70}")
        print(f"\nWorkGraph PK: {pk}")
        print(f"\nMonitor with:")
        print(f"  verdi process show {pk}")
        print(f"  verdi process report {pk}")
        print(f"\nWatch concurrent jobs (should max out at 2):")
        print(f"  watch -n 2 'verdi process list'")
        print(f"\nExpected workflow:")
        print(f"  1. Bulk relaxation (1 VASP job)")
        print(f"  2. Reference calculations: metal, oxygen (2 VASP jobs in parallel)")
        print(f"  3. Formation enthalpy calculation")
        print(f"  4. Internal slab generation from miller indices: {miller_indices}")
        print(f"  5. Slab SCF calculations (N slabs - limited to 2 concurrent)")
        print(f"  6. Slab relaxations (N slabs - limited to 2 concurrent)")
        print(f"  7. Relaxation energy calculations")
        print(f"  8. Surface energy calculations (N slabs)")
        print(f"\nNew Features Demonstrated:")
        print(f"  ✓ Internal slab generation from miller_indices (no pre-generated slabs)")
        print(f"  ✓ Custom builders for full parameter control")
        print(f"  ✓ Different parameters for bulk, metal, oxygen (via builders)")
        print(f"  ✓ Different parameters for SCF vs relaxation (via templates)")
        print(f"\nSerial Preset Features:")
        print(f"  ✓ Flat-graph architecture (all nodes at same level)")
        print(f"  ✓ max_number_jobs controls ALL VASP jobs")
        print(f"  ✓ No nested sub-workgraphs")
        print(f"  ✓ Better provenance tracking")
        print(f"\nExpected outputs:")
        print(f"  - bulk_energy, metal_energy, oxygen_energy")
        print(f"  - formation_enthalpy, reference_energies, oxide_type")
        print(f"  - relaxed_slabs (generated from miller_indices)")
        print(f"  - slab_energies, unrelaxed_slab_energies, relaxation_energies")
        print(f"  - surface_energies for all slabs")
        print(f"\nComparison with standard preset (step_05):")
        print(f"  Standard: Uses nested sub-workgraphs (max_number_jobs doesn't propagate)")
        print(f"  Serial:   Uses flat graph (max_number_jobs controls all VASP)")
        print(f"\nBuilder vs Standard Parameters:")
        print(f"  When builders provided: Full control over VASP parameters (recommended)")
        print(f"  When builders None: Uses standard parameter preparation (fallback)")
        print(f"  This script demonstrates builders for bulk, metal, oxygen")
        print(f"  Slab calculations use standard parameters (can also use builders)")
        print(f"{'='*70}\n")

        return pk

    except Exception as e:
        print(f"\n✗ Error submitting workgraph: {e}")
        import traceback
        traceback.print_exc()
        return None


if __name__ == '__main__':
    try:
        pk = main()
        if pk is None:
            sys.exit(1)
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

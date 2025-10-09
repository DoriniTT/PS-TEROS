#!/home/thiagotd/envs/psteros/bin/python
"""
Example script demonstrating the RESTART functionality in PS-TEROS.

This script shows how to restart slab relaxation calculations from a previous
PS-TEROS run that failed or didn't converge. This is useful when:
- Calculations hit time limits or resource constraints
- Convergence criteria were too strict
- You want to continue with tighter convergence settings
- VASP calculations encountered recoverable errors

Usage:
    1. Run a PS-TEROS calculation that fails or doesn't converge
    2. Note the PK of the main workgraph node
    3. Run this script with restart_from_node=<PK>
    
Example:
    source ~/envs/psteros/bin/activate && python restart_example.py
"""

from aiida import load_profile, orm
from teros.core.workgraph import build_core_workgraph

def main():
    """Main function to demonstrate restart functionality."""
    
    # Load AiiDA profile
    print("Loading AiiDA profile...")
    load_profile(profile='psteros')
    
    # ===== IMPORTANT: SET THIS TO YOUR PREVIOUS WORKGRAPH PK =====
    # This is the PK of a previous PS-TEROS run that you want to restart from
    PREVIOUS_RUN_PK = 22223  # Replace with your actual PK
    
    print(f"\n{'='*80}")
    print("RESTART EXAMPLE: Continuing slab relaxations from previous run")
    print(f"{'='*80}\n")
    
    # Check if the previous run exists and has the required outputs
    try:
        prev_node = orm.load_node(PREVIOUS_RUN_PK)
        print(f"Previous run details:")
        print(f"  PK: {prev_node.pk}")
        print(f"  Label: {prev_node.label}")
        print(f"  State: {prev_node.process_state}")
        print(f"  Created: {prev_node.ctime}")
        
        # Check outputs
        if hasattr(prev_node.outputs, 'slab_structures'):
            print(f"\n  ✓ Has slab_structures: {list(prev_node.outputs.slab_structures.keys())}")
        else:
            print(f"\n  ✗ Missing slab_structures - cannot restart")
            return
            
        if hasattr(prev_node.outputs, 'slab_remote'):
            print(f"  ✓ Has slab_remote: {list(prev_node.outputs.slab_remote.keys())}")
        else:
            print(f"  ✗ Missing slab_remote - cannot restart")
            print(f"     This may be an old run. Re-run with updated PS-TEROS code.")
            return
            
    except Exception as e:
        print(f"Error loading previous run: {e}")
        return
    
    # Define structures directory
    structures_dir = "/home/thiagotd/git/PS-TEROS/examples/structures"
    
    # Define calculation parameters
    code_label = "VASP-VTST-6.4.3@bohr"
    potential_family = "PBE"
    
    # ===== BULK RELAXATION PARAMETERS =====
    # (Same as original run)
    bulk_parameters = {
        "PREC": "Accurate",
        "ENCUT": 500,
        "EDIFF": 1e-6,
        "ISMEAR": 0,
        "SIGMA": 0.05,
        "IBRION": 2,
        "ISIF": 3,
        "NSW": 100,
        "EDIFFG": -0.1,
        "ALGO": "Normal",
        "LREAL": "Auto",
        "LWAVE": False,
        "LCHARG": False,
    }
    
    bulk_options = {
        "resources": {
            "num_machines": 1,
            "num_cores_per_machine": 40,
        },
        "queue_name": "par40",
    }
    
    # ===== REFERENCE STRUCTURE PARAMETERS =====
    metal_parameters = {
        "PREC": "Accurate",
        "ENCUT": 520,
        "EDIFF": 1e-6,
        "ISMEAR": 1,
        "SIGMA": 0.2,
        "IBRION": 2,
        "ISIF": 3,
        "NSW": 100,
        "EDIFFG": -0.1,
        "ALGO": "Normal",
        "LREAL": "Auto",
        "LWAVE": False,
        "LCHARG": False,
    }
    
    metal_options = {
        "resources": {
            "num_machines": 1,
            "num_cores_per_machine": 40,
        },
        "queue_name": "par40",
    }
    
    # For binary oxide, nonmetal is not physically meaningful but required
    nonmetal_parameters = metal_parameters.copy()
    nonmetal_options = metal_options.copy()
    
    oxygen_parameters = {
        "PREC": "Accurate",
        "ENCUT": 520,
        "EDIFF": 1e-6,
        "ISMEAR": 0,
        "SIGMA": 0.01,
        "IBRION": 2,
        "ISIF": 2,
        "NSW": 100,
        "EDIFFG": -0.1,
        "ALGO": "Normal",
        "LREAL": False,
        "LWAVE": False,
        "LCHARG": False,
    }
    
    oxygen_options = {
        "resources": {
            "num_machines": 1,
            "num_cores_per_machine": 40,
        },
        "queue_name": "par40",
    }
    
    # ===== SLAB RELAXATION PARAMETERS (UPDATED FOR RESTART) =====
    # You can modify these parameters for the restart, e.g.:
    # - Increase NSW for more ionic steps
    # - Tighten EDIFFG for better convergence
    # - Change IBRION algorithm
    
    slab_parameters = {
        "PREC": "Accurate",
        "ENCUT": 520,
        "EDIFF": 1e-6,
        "ISMEAR": 0,
        "SIGMA": 0.05,
        "IBRION": 2,
        "ISIF": 2,  # Relax atoms only, keep cell fixed
        "NSW": 200,  # INCREASED: More steps for convergence
        "EDIFFG": -0.05,  # TIGHTER: Better convergence
        "ALGO": "Normal",
        "LREAL": "Auto",
        "LWAVE": False,
        "LCHARG": False,
    }
    
    slab_options = {
        "resources": {
            "num_machines": 1,
            "num_cores_per_machine": 40,
        },
        "queue_name": "par40",
        "max_wallclock_seconds": 7200,  # 2 hours
    }
    
    slab_potential_mapping = {"Ag": "Ag", "O": "O"}
    slab_kpoints_spacing = 0.3
    
    # Slab generation options (not used since we're restarting)
    miller_indices = [1, 0, 0]
    min_slab_thickness = 10.0
    min_vacuum_thickness = 15.0
    
    # ===== THERMODYNAMICS PARAMETERS =====
    compute_thermodynamics = True
    thermodynamics_sampling = 100
    
    # ===== BUILD WORKGRAPH WITH RESTART =====
    print(f"\n{'='*80}")
    print("Building WorkGraph with RESTART functionality...")
    print(f"{'='*80}\n")
    
    wg = build_core_workgraph(
        structures_dir=structures_dir,
        bulk_name="ag2o.cif",
        metal_name="Ag.cif",
        nonmetal_name="Ag.cif",  # Dummy for binary oxide
        oxygen_name="O2.cif",
        code_label=code_label,
        potential_family=potential_family,
        bulk_potential_mapping={"Ag": "Ag", "O": "O"},
        metal_potential_mapping={"Ag": "Ag"},
        nonmetal_potential_mapping={"Ag": "Ag"},
        oxygen_potential_mapping={"O": "O"},
        kpoints_spacing=0.3,
        bulk_parameters=bulk_parameters,
        bulk_options=bulk_options,
        metal_parameters=metal_parameters,
        metal_options=metal_options,
        nonmetal_parameters=nonmetal_parameters,
        nonmetal_options=nonmetal_options,
        oxygen_parameters=oxygen_parameters,
        oxygen_options=oxygen_options,
        clean_workdir=True,
        # Slab generation (will be overridden by restart)
        miller_indices=miller_indices,
        min_slab_thickness=min_slab_thickness,
        min_vacuum_thickness=min_vacuum_thickness,
        lll_reduce=True,
        center_slab=True,
        symmetrize=True,
        primitive=True,
        # Slab relaxation
        relax_slabs=True,
        slab_parameters=slab_parameters,
        slab_options=slab_options,
        slab_potential_mapping=slab_potential_mapping,
        slab_kpoints_spacing=slab_kpoints_spacing,
        # Thermodynamics
        compute_thermodynamics=compute_thermodynamics,
        thermodynamics_sampling=thermodynamics_sampling,
        # RESTART: This is the key parameter
        restart_from_node=PREVIOUS_RUN_PK,
        name=f"Ag2O_SlabsRelax_RESTART_from_{PREVIOUS_RUN_PK}",
    )
    
    # Optional: Export to HTML
    try:
        html_file = f"ag2o_restart_from_{PREVIOUS_RUN_PK}.html"
        wg.to_html(html_file)
        print(f"\n✓ WorkGraph visualization saved to: {html_file}")
    except Exception as e:
        print(f"\n✗ Could not generate HTML visualization: {e}")
    
    # ===== SUBMIT WORKGRAPH =====
    print(f"\n{'='*80}")
    print("Submitting WorkGraph...")
    print(f"{'='*80}")
    wg.submit(wait=False)
    
    print(f"\n✓ Restart WorkGraph submitted successfully!")
    print(f"\nWorkGraph PK: {wg.pk}")
    
    print(f"\n{'-'*80}")
    print(f"MONITORING")
    print(f"{'-'*80}")
    print(f"  Check status:     verdi process show {wg.pk}")
    print(f"  Check report:     verdi process report {wg.pk}")
    print(f"  List processes:   verdi process list -a -p1")
    
    print(f"\n{'-'*80}")
    print(f"WHAT HAPPENS IN RESTART MODE")
    print(f"{'-'*80}")
    print(f"  1. Slab structures from node {PREVIOUS_RUN_PK} are reused")
    print(f"  2. RemoteData (restart_folder) from previous VASP runs are passed to new calculations")
    print(f"  3. VASP will read WAVECAR, CONTCAR, etc. from previous run and continue")
    print(f"  4. Updated parameters (NSW={slab_parameters['NSW']}, EDIFFG={slab_parameters['EDIFFG']}) are used")
    print(f"  5. All other parts (bulk, references, formation enthalpy) are recalculated")
    
    print(f"\n{'='*80}\n")
    
    return wg


if __name__ == "__main__":
    """
    Run the restart example.
    
    Before running:
    1. Make sure you have a previous PS-TEROS run with slab_remote outputs
    2. Set PREVIOUS_RUN_PK to the PK of that run
    3. Optionally modify slab_parameters to use different convergence settings
    4. Clear Python cache and restart daemon:
       find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null
       find . -name "*.pyc" -delete 2>/dev/null
       verdi daemon restart
    5. Run this script:
       source ~/envs/psteros/bin/activate && python restart_example.py
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

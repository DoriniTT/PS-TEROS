#!/usr/bin/env python
"""
How to access B-based results from a completed workflow

After running your workflow again with the updated code, you can access
both A-based and B-based results as shown below.
"""

import aiida
aiida.load_profile('psteros')

from aiida import orm

# Replace with your new workflow PK after rerunning
WORKFLOW_PK = 29082  # <-- Change this to your NEW workflow PK

def show_both_formulations(workflow_pk):
    """Show how to access both A-based and B-based results."""
    
    # Load the workflow
    node = orm.load_node(workflow_pk)
    
    print(f"\n{'='*70}")
    print(f"Workflow PK: {workflow_pk}")
    print(f"State: {node.process_state}")
    print(f"{'='*70}\n")
    
    # Get surface energy outputs
    surface_energies = node.outputs.surface_energies
    
    # Get termination keys
    term_keys = [key for key in dir(surface_energies) if key.startswith('term_')]
    
    # Iterate through all terminations
    for term_key in term_keys:
        print(f"\n{'-'*70}")
        print(f"Termination: {term_key}")
        print(f"{'-'*70}")
        
        term_node = getattr(surface_energies, term_key)
        term_dict = term_node.get_dict()
        
        # Check if new format with A_based and B_based
        if 'A_based' in term_dict and 'B_based' in term_dict:
            print("\n✓ New format detected with both formulations!")
            
            # A-based results
            a_based = term_dict['A_based']
            print(f"\n  A-based formulation:")
            print(f"    Element A (independent): {a_based['element_A_independent']}")
            print(f"    Element B (reference):   {a_based['element_B_reference']}")
            print(f"    Gamma_A:  {a_based['Gamma_A']:+.6f} atoms/Å²")
            print(f"    Gamma_O:  {a_based['Gamma_O']:+.6f} atoms/Å²")
            print(f"    phi:      {a_based['phi']:.6f} eV/Å²")
            print(f"    Δμ_A range: [{min(a_based['delta_mu_A_range']):.3f}, {max(a_based['delta_mu_A_range']):.3f}] eV")
            
            # B-based results
            b_based = term_dict['B_based']
            print(f"\n  B-based formulation:")
            print(f"    Element B (independent): {b_based['element_B_independent']}")
            print(f"    Element A (reference):   {b_based['element_A_reference']}")
            print(f"    Gamma_B:  {b_based['Gamma_B']:+.6f} atoms/Å²")
            print(f"    Gamma_O:  {b_based['Gamma_O']:+.6f} atoms/Å²")
            print(f"    phi:      {b_based['phi']:.6f} eV/Å²")
            print(f"    Δμ_B range: [{min(b_based['delta_mu_B_range']):.3f}, {max(b_based['delta_mu_B_range']):.3f}] eV")
            
            # Common info
            print(f"\n  System information:")
            print(f"    Area:       {term_dict['area_A2']:.3f} Å²")
            print(f"    Stoich:     {term_dict['bulk_stoichiometry_AxByOz']}")
            print(f"    Slab atoms: {term_dict['slab_atom_counts']}")
            
        else:
            print("\n✗ Old format detected (workflow ran before code update)")
            print("   Please rerun the workflow to get B-based results.")
            print(f"\n  Current results (A-based only):")
            print(f"    Element M: {term_dict.get('element_M_independent', 'N/A')}")
            print(f"    Gamma_M:   {term_dict.get('Gamma_M_vs_Nref', 0):+.6f} atoms/Å²")
            print(f"    Gamma_O:   {term_dict.get('Gamma_O_vs_Nref', 0):+.6f} atoms/Å²")
            print(f"    phi:       {term_dict.get('phi', 0):.6f} eV/Å²")
    
    print(f"\n{'='*70}\n")


def export_for_plotting(workflow_pk, termination='term_0'):
    """
    Export data for plotting phase diagrams.
    
    Returns dictionaries ready for use with matplotlib/plotly.
    """
    node = orm.load_node(workflow_pk)
    term_node = getattr(node.outputs.surface_energies, termination)
    term_dict = term_node.get_dict()
    
    if 'A_based' not in term_dict or 'B_based' not in term_dict:
        print("Warning: Old format detected. Please rerun workflow.")
        return None, None
    
    # Prepare A-based data for plotting
    a_data = {
        'delta_mu_A': term_dict['A_based']['delta_mu_A_range'],
        'delta_mu_O': term_dict['A_based']['delta_mu_O_range'],
        'gamma_grid': term_dict['A_based']['gamma_grid'],
        'element_A': term_dict['A_based']['element_A_independent'],
    }
    
    # Prepare B-based data for plotting
    b_data = {
        'delta_mu_B': term_dict['B_based']['delta_mu_B_range'],
        'delta_mu_O': term_dict['B_based']['delta_mu_O_range'],
        'gamma_grid': term_dict['B_based']['gamma_grid'],
        'element_B': term_dict['B_based']['element_B_independent'],
    }
    
    return a_data, b_data


if __name__ == '__main__':
    print("\n" + "="*70)
    print("ACCESSING B-BASED RESULTS FROM WORKFLOW")
    print("="*70)
    
    print("\nIMPORTANT:")
    print("  Your current workflow (PK 29082) ran BEFORE the code update.")
    print("  To see B-based results, you need to:")
    print("    1. Run the workflow again:")
    print("       python examples/complete/complete_ag3po4_example.py")
    print("    2. Wait for completion")
    print("    3. Update WORKFLOW_PK in this script to the new PK")
    print("    4. Run this script again")
    print("\n  The daemon has been restarted and will use the new code.")
    
    print("\n" + "="*70)
    print("CHECKING CURRENT WORKFLOW (OLD VERSION)")
    print("="*70)
    
    try:
        show_both_formulations(WORKFLOW_PK)
    except Exception as e:
        print(f"\nError: {e}")
        print("\nMake sure the workflow PK is correct and the workflow has finished.")
    
    print("\nFor plotting, after rerunning the workflow:")
    print("  a_data, b_data = export_for_plotting(NEW_PK, 'term_0')")
    print("  # Then use example_usage_B_formulation.py functions to plot")

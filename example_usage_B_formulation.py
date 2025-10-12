#!/usr/bin/env python
"""
Example: Using both A-based and B-based formulations for ternary oxides

This example shows how to access and plot both formulations from the
calculate_surface_energy_ternary function results.
"""

import numpy as np
import matplotlib.pyplot as plt

def plot_both_formulations(result_dict):
    """
    Create side-by-side plots of A-based and B-based surface energy phase diagrams.
    
    Args:
        result_dict: Dictionary returned from calculate_surface_energy_ternary
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # ===== A-based formulation =====
    a_based = result_dict['A_based']
    
    # Convert to numpy arrays
    delta_mu_A = np.array(a_based['delta_mu_A_range'])
    delta_mu_O = np.array(a_based['delta_mu_O_range'])
    gamma_A = np.array(a_based['gamma_grid'])
    
    # Create meshgrid for contour plot
    MU_O, MU_A = np.meshgrid(delta_mu_O, delta_mu_A)
    
    # Plot A-based phase diagram
    contour_A = ax1.contourf(MU_O, MU_A, gamma_A, levels=20, cmap='viridis')
    ax1.set_xlabel(f'Δμ_O (eV)', fontsize=12)
    ax1.set_ylabel(f'Δμ_{a_based["element_A_independent"]} (eV)', fontsize=12)
    ax1.set_title(f'A-based: γ(Δμ_{a_based["element_A_independent"]}, Δμ_O)', fontsize=13, fontweight='bold')
    cbar1 = plt.colorbar(contour_A, ax=ax1)
    cbar1.set_label('Surface Energy γ (eV/Å²)', fontsize=11)
    ax1.grid(True, alpha=0.3)
    
    # ===== B-based formulation =====
    b_based = result_dict['B_based']
    
    # Convert to numpy arrays
    delta_mu_B = np.array(b_based['delta_mu_B_range'])
    delta_mu_O_B = np.array(b_based['delta_mu_O_range'])
    gamma_B = np.array(b_based['gamma_grid'])
    
    # Create meshgrid for contour plot
    MU_O_B, MU_B = np.meshgrid(delta_mu_O_B, delta_mu_B)
    
    # Plot B-based phase diagram
    contour_B = ax2.contourf(MU_O_B, MU_B, gamma_B, levels=20, cmap='viridis')
    ax2.set_xlabel(f'Δμ_O (eV)', fontsize=12)
    ax2.set_ylabel(f'Δμ_{b_based["element_B_independent"]} (eV)', fontsize=12)
    ax2.set_title(f'B-based: γ(Δμ_{b_based["element_B_independent"]}, Δμ_O)', fontsize=13, fontweight='bold')
    cbar2 = plt.colorbar(contour_B, ax=ax2)
    cbar2.set_label('Surface Energy γ (eV/Å²)', fontsize=11)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Add system information
    stoich = result_dict['bulk_stoichiometry_AxByOz']
    formula = f"{a_based['element_A_independent']}_{stoich['x_'+a_based['element_A_independent']]}"
    formula += f"{b_based['element_B_independent']}_{stoich['y_'+b_based['element_B_independent']]}"
    formula += f"O_{stoich['z_O']}"
    
    fig.suptitle(f'Surface Energy Phase Diagrams: {formula}', 
                 fontsize=14, fontweight='bold', y=1.02)
    
    return fig


def plot_1d_slices(result_dict):
    """
    Plot 1D slices of surface energy at specific conditions.
    
    Args:
        result_dict: Dictionary returned from calculate_surface_energy_ternary
    """
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))
    
    a_based = result_dict['A_based']
    b_based = result_dict['B_based']
    
    # A-based: γ vs Δμ_O at Δμ_A = 0
    ax1.plot(a_based['delta_mu_O_range'], a_based['gamma_at_muA_zero'], 
             'b-', linewidth=2, label=f'Δμ_{a_based["element_A_independent"]} = 0')
    ax1.set_xlabel('Δμ_O (eV)', fontsize=11)
    ax1.set_ylabel('Surface Energy γ (eV/Å²)', fontsize=11)
    ax1.set_title(f'A-based: γ vs Δμ_O at {a_based["element_A_independent"]}-rich limit', fontsize=12)
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # A-based: γ vs Δμ_A at Δμ_O = 0
    ax2.plot(a_based['delta_mu_A_range'], a_based['gamma_at_muO_zero'], 
             'r-', linewidth=2, label='Δμ_O = 0')
    ax2.set_xlabel(f'Δμ_{a_based["element_A_independent"]} (eV)', fontsize=11)
    ax2.set_ylabel('Surface Energy γ (eV/Å²)', fontsize=11)
    ax2.set_title('A-based: γ vs Δμ_A at O-rich limit', fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # B-based: γ vs Δμ_O at Δμ_B = 0
    ax3.plot(b_based['delta_mu_O_range'], b_based['gamma_at_muB_zero'], 
             'g-', linewidth=2, label=f'Δμ_{b_based["element_B_independent"]} = 0')
    ax3.set_xlabel('Δμ_O (eV)', fontsize=11)
    ax3.set_ylabel('Surface Energy γ (eV/Å²)', fontsize=11)
    ax3.set_title(f'B-based: γ vs Δμ_O at {b_based["element_B_independent"]}-rich limit', fontsize=12)
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    
    # B-based: γ vs Δμ_B at Δμ_O = 0
    ax4.plot(b_based['delta_mu_B_range'], b_based['gamma_at_muO_zero'], 
             'm-', linewidth=2, label='Δμ_O = 0')
    ax4.set_xlabel(f'Δμ_{b_based["element_B_independent"]} (eV)', fontsize=11)
    ax4.set_ylabel('Surface Energy γ (eV/Å²)', fontsize=11)
    ax4.set_title('B-based: γ vs Δμ_B at O-rich limit', fontsize=12)
    ax4.grid(True, alpha=0.3)
    ax4.legend()
    
    plt.tight_layout()
    
    return fig


def compare_surface_excesses(result_dict):
    """
    Print comparison of surface excesses for both formulations.
    
    Args:
        result_dict: Dictionary returned from calculate_surface_energy_ternary
    """
    print("\n" + "="*70)
    print("SURFACE EXCESS COMPARISON")
    print("="*70)
    
    a_based = result_dict['A_based']
    b_based = result_dict['B_based']
    
    print(f"\nA-based formulation ({a_based['element_A_independent']} as independent):")
    print(f"  Γ_{a_based['element_A_independent']} = {a_based['Gamma_A']:+.6f} atoms/Å²")
    print(f"  Γ_O = {a_based['Gamma_O']:+.6f} atoms/Å²")
    print(f"  φ = {a_based['phi']:.6f} eV/Å²")
    
    print(f"\nB-based formulation ({b_based['element_B_independent']} as independent):")
    print(f"  Γ_{b_based['element_B_independent']} = {b_based['Gamma_B']:+.6f} atoms/Å²")
    print(f"  Γ_O = {b_based['Gamma_O']:+.6f} atoms/Å²")
    print(f"  φ = {b_based['phi']:.6f} eV/Å²")
    
    print("\nInterpretation:")
    if a_based['Gamma_A'] > 0:
        print(f"  • Surface is {a_based['element_A_independent']}-deficient (relative to bulk stoichiometry)")
    else:
        print(f"  • Surface is {a_based['element_A_independent']}-rich (relative to bulk stoichiometry)")
    
    if b_based['Gamma_B'] > 0:
        print(f"  • Surface is {b_based['element_B_independent']}-deficient (relative to bulk stoichiometry)")
    else:
        print(f"  • Surface is {b_based['element_B_independent']}-rich (relative to bulk stoichiometry)")
    
    print("="*70)


# ===== Example usage =====
if __name__ == '__main__':
    print("This is an example script showing how to use both formulations.")
    print("To use with actual data, call the functions above with your result_dict:")
    print()
    print("  result = calculate_surface_energy_ternary(...)")
    print("  result_dict = result.result.get_dict()  # If using from workflow")
    print("  # or")
    print("  result_dict = result.get_dict()  # If calling calcfunction directly")
    print()
    print("  # Then plot:")
    print("  fig1 = plot_both_formulations(result_dict)")
    print("  fig2 = plot_1d_slices(result_dict)")
    print("  compare_surface_excesses(result_dict)")
    print("  plt.show()  # or fig1.savefig('phase_diagram.png')")

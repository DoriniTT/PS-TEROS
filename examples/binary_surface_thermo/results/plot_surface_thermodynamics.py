#!/usr/bin/env python3
"""
Surface Thermodynamics Analysis for SnO2 (110) - Binary Oxide
==============================================================

Generates surface energy plots from PS-TEROS workflow outputs:
- gamma vs Delta_mu_O (1D plot with temperature scale)

Material: SnO2 (rutile tin dioxide)
Surface: (110) - the most stable rutile surface

For binary oxides, the surface energy depends only on the oxygen
chemical potential:
    gamma(Delta_mu_O) = phi - Gamma_O * Delta_mu_O

Usage:
    source ~/envs/aiida/bin/activate
    python plot_surface_thermodynamics.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from aiida import load_profile
from aiida.orm import load_node

# Set seaborn style for better aesthetics
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.2)

# =============================================================================
# CONFIGURATION
# =============================================================================

AIIDA_PROFILE = 'presto'

# Workflow PK or UUID
WORKFLOW_PK = 20372

# Material info
MATERIAL_NAME = 'SnO2'
SURFACE_MILLER = '(110)'

# Output configuration
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_FORMAT = 'pdf'

# Unit conversion
EV_PER_A2_TO_J_PER_M2 = 16.0217663

# Plot aesthetics
FIGURE_DPI = 150
FONTSIZE_LABEL = 14
FONTSIZE_TICK = 12
FONTSIZE_LEGEND = 10
FONTSIZE_TITLE = 16

# Color palette for terminations
TERMINATION_COLORS = sns.color_palette("deep", 10)

# =============================================================================
# THERMODYNAMIC CONSTANTS FOR T-P DEPENDENCE
# =============================================================================

# Boltzmann constant in eV/K
KB_EV = 8.617333262e-5

# Standard entropy of O2 at 298K (from NIST-JANAF)
# S(O2, 298K) = 205.15 J/(mol*K) = 0.002126 eV/K (per molecule)
S_O2_298K = 0.002126  # eV/K per O2 molecule

# Reference Delta_mu_O at 298K, 1 atm (typical value)
DELTA_MU_O_REF = -0.27  # eV


def T_from_delta_mu_O(delta_mu_O, P_atm=1.0, delta_mu_O_ref=DELTA_MU_O_REF):
    """
    Calculate temperature from Delta_mu_O at given pressure.

    For ideal gas O2:
        Delta_mu_O(T,P) = Delta_mu_O(T_ref) - (T - T_ref) * S_O / 2 + (kT/2) * ln(P/P_ref)

    Simplified inversion (neglects pressure term's T dependence):
        T = T_ref - 2 * (Delta_mu_O - Delta_mu_O_ref) / S_O2

    Args:
        delta_mu_O: Chemical potential deviation (eV)
        P_atm: Pressure in atm
        delta_mu_O_ref: Reference Delta_mu_O at 298K, 1 atm (eV)

    Returns:
        Temperature in K

    Note:
        With the Reuter & Scheffler convention:
        - Delta_mu_O = 0 corresponds to low T (O-rich, equilibrium with O2)
        - Delta_mu_O < 0 corresponds to high T (O-poor, reducing conditions)
    """
    T_ref = 298.15
    T = T_ref - 2 * (delta_mu_O - delta_mu_O_ref) / S_O2_298K
    return max(T, 100)  # Minimum 100K


# =============================================================================
# DATA LOADING FUNCTIONS
# =============================================================================

def load_surface_energies(workflow_pk: int) -> dict:
    """
    Load surface energy data from a completed PS-TEROS workflow.

    Args:
        workflow_pk: AiiDA PK of the workflow

    Returns:
        Dictionary with termination data
    """
    wg = load_node(workflow_pk)

    if wg.process_state.value != 'finished':
        raise RuntimeError(f"Workflow {workflow_pk} not finished: {wg.process_state.value}")

    if 'surface_energies' not in wg.outputs:
        raise KeyError(f"Workflow has no 'surface_energies' output")

    result = {}
    for term_label in wg.outputs.surface_energies:
        term_data = wg.outputs.surface_energies[term_label].get_dict()

        # For binary oxides, data is already computed
        # Use 'primary' dict if available, otherwise use top-level keys
        data_source = term_data.get('primary', term_data)

        result[term_label] = {
            'delta_mu_O_range': np.array(data_source['delta_mu_O_range']),
            'gamma_array': np.array(data_source['gamma_array']),
            'phi': data_source['phi'],
            'Gamma_O': data_source['Gamma_O'],
            'gamma_O_rich': data_source['gamma_O_rich'],
            'gamma_O_poor': data_source['gamma_O_poor'],
            'element_M': data_source.get('element_M', 'Sn'),
            # Additional info
            'area_A2': term_data.get('area_A2'),
            'E_slab_eV': term_data.get('E_slab_eV'),
            'slab_atom_counts': term_data.get('slab_atom_counts', {}),
        }

    return result


# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

def plot_gamma_vs_mu_O(
    material_name: str,
    surface_miller: str,
    all_term_data: dict,
    output_path: str
) -> None:
    """
    Create 1D plot of surface energy vs Delta_mu_O for binary oxide.

    Includes secondary axis for Temperature @ 1 atm.

    Args:
        material_name: Name of the material (e.g., 'SnO2')
        surface_miller: Miller indices string (e.g., '(110)')
        all_term_data: Dictionary of termination data
        output_path: Path to save the figure
    """
    fig, ax = plt.subplots(figsize=(8, 6), dpi=FIGURE_DPI)

    term_labels = sorted(all_term_data.keys())

    for i, term_label in enumerate(term_labels):
        data = all_term_data[term_label]

        delta_mu_O = data['delta_mu_O_range']
        gamma = data['gamma_array']

        color = TERMINATION_COLORS[i % len(TERMINATION_COLORS)]
        # Label format: term_0 -> T1, term_1 -> T2, etc.
        term_number = int(term_label.split('_')[1]) + 1
        label = f'T{term_number}'
        ax.plot(delta_mu_O, gamma, '-', color=color, linewidth=2.0, label=label)

    ax.set_xlabel(r'$\Delta\mu_{\rm O}$ (eV)', fontsize=FONTSIZE_LABEL)
    ax.set_ylabel(r'$\gamma$ (J/m$^2$)', fontsize=FONTSIZE_LABEL)
    ax.tick_params(axis='both', labelsize=FONTSIZE_TICK)

    # Get x-axis limits from data
    x_min, x_max = ax.get_xlim()

    # Define sensible tick positions for Delta_mu_O
    # With corrected convention: range is approximately -3 to 0 eV
    mu_ticks = np.arange(np.ceil(x_min * 2) / 2, np.floor(x_max * 2) / 2 + 0.25, 0.5)
    mu_ticks = mu_ticks[(mu_ticks >= x_min - 0.1) & (mu_ticks <= x_max + 0.1)]

    # --- Secondary axis: Temperature @ 1 atm ---
    ax_temp = ax.secondary_xaxis('top')
    ax_temp.set_xlabel('Temperature (K) @ 1 atm', fontsize=FONTSIZE_LABEL - 2)

    # Calculate temperatures using the thermodynamic relationship
    # With Reuter & Scheffler convention:
    # - Delta_mu_O = 0 corresponds to low T (O-rich)
    # - Delta_mu_O < 0 corresponds to high T (O-poor)
    t_labels = []
    for mu in mu_ticks:
        T = T_from_delta_mu_O(mu, P_atm=1.0)
        t_labels.append(f'{int(T)}')

    ax_temp.set_xticks(mu_ticks)
    ax_temp.set_xticklabels(t_labels, fontsize=FONTSIZE_TICK - 1)

    # Legend
    ax.legend(fontsize=FONTSIZE_LEGEND, loc='best', framealpha=0.95,
              edgecolor='gray', fancybox=False, title=f'{material_name} {surface_miller}')

    # Remove top and right spines (open axes style)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.8)
    ax.spines['bottom'].set_linewidth(0.8)

    plt.tight_layout()
    fig.subplots_adjust(top=0.88)  # Make room for temperature axis
    plt.savefig(output_path, format=OUTPUT_FORMAT, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {output_path}")


def plot_stability_diagram(
    material_name: str,
    surface_miller: str,
    all_term_data: dict,
    output_path: str
) -> None:
    """
    Create stability diagram showing which termination is most stable
    as a function of Delta_mu_O.

    Args:
        material_name: Name of the material
        surface_miller: Miller indices string
        all_term_data: Dictionary of termination data
        output_path: Path to save the figure
    """
    fig, ax = plt.subplots(figsize=(10, 4), dpi=FIGURE_DPI)

    term_labels = sorted(all_term_data.keys())
    n_terms = len(term_labels)

    # Get common x-axis (all should have same range)
    first_data = all_term_data[term_labels[0]]
    delta_mu_O = first_data['delta_mu_O_range']

    # Build gamma array for all terminations
    gamma_all = np.zeros((n_terms, len(delta_mu_O)))
    for i, term_label in enumerate(term_labels):
        gamma_all[i] = all_term_data[term_label]['gamma_array']

    # Find minimum gamma at each point
    stable_idx = np.argmin(gamma_all, axis=0)

    # Create colored regions
    colors = TERMINATION_COLORS[:n_terms]

    # Find transition points
    transitions = [0]
    for i in range(1, len(stable_idx)):
        if stable_idx[i] != stable_idx[i-1]:
            transitions.append(i)
    transitions.append(len(stable_idx) - 1)

    # Plot colored bars for each stable region
    y_height = 1.0
    for j in range(len(transitions) - 1):
        i_start = transitions[j]
        i_end = transitions[j + 1]
        term_idx = stable_idx[i_start]
        term_number = int(term_labels[term_idx].split('_')[1]) + 1

        x_start = delta_mu_O[i_start]
        x_end = delta_mu_O[i_end]

        ax.fill_between([x_start, x_end], 0, y_height,
                        color=colors[term_idx], alpha=0.7,
                        edgecolor='white', linewidth=1)

        # Add label in center
        x_center = (x_start + x_end) / 2
        ax.text(x_center, y_height / 2, f'T{term_number}',
                ha='center', va='center', fontsize=FONTSIZE_LEGEND + 2,
                fontweight='bold', color='white')

    ax.set_xlim(delta_mu_O[0], delta_mu_O[-1])
    ax.set_ylim(0, y_height)
    ax.set_xlabel(r'$\Delta\mu_{\rm O}$ (eV)', fontsize=FONTSIZE_LABEL)
    ax.set_yticks([])
    ax.set_ylabel('Stable\nTermination', fontsize=FONTSIZE_LABEL - 2, rotation=0,
                  ha='right', va='center')

    # Add O-poor / O-rich labels
    # With Reuter & Scheffler convention:
    # - Left (negative Delta_mu_O) = O-poor = high T = reducing
    # - Right (Delta_mu_O = 0) = O-rich = low T = oxidizing
    ax.text(delta_mu_O[0], -0.15, 'O-poor\n(high T)', ha='left', va='top',
            fontsize=FONTSIZE_TICK - 2, transform=ax.get_xaxis_transform())
    ax.text(delta_mu_O[-1], -0.15, 'O-rich\n(low T)', ha='right', va='top',
            fontsize=FONTSIZE_TICK - 2, transform=ax.get_xaxis_transform())

    ax.set_title(f'{material_name} {surface_miller} - Termination Stability',
                 fontsize=FONTSIZE_TITLE)

    # Style
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_path, format=OUTPUT_FORMAT, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {output_path}")


def print_summary_table(all_term_data: dict) -> None:
    """
    Print a summary table of surface energies for all terminations.
    """
    print("\n" + "=" * 70)
    print("SURFACE ENERGY SUMMARY")
    print("=" * 70)
    print(f"{'Term':<8} {'phi (J/m2)':<12} {'gamma_O-poor':<14} {'gamma_O-rich':<14} {'Gamma_O':<12}")
    print("-" * 70)

    term_labels = sorted(all_term_data.keys())
    for term_label in term_labels:
        data = all_term_data[term_label]
        term_number = int(term_label.split('_')[1]) + 1
        print(f"T{term_number:<7} {data['phi']:<12.4f} {data['gamma_O_poor']:<14.4f} "
              f"{data['gamma_O_rich']:<14.4f} {data['Gamma_O']:<12.6f}")

    print("-" * 70)
    print("Units: gamma in J/m^2, Gamma_O in atoms/A^2")
    print("O-poor: high temperature / reducing conditions")
    print("O-rich: low temperature / oxidizing conditions")
    print("=" * 70 + "\n")


def export_data_to_csv(all_term_data: dict, output_path: str) -> None:
    """
    Export surface energy data to CSV file.
    """
    term_labels = sorted(all_term_data.keys())

    # Get common x-axis
    delta_mu_O = all_term_data[term_labels[0]]['delta_mu_O_range']

    # Build header
    header = ['delta_mu_O_eV']
    for term_label in term_labels:
        term_number = int(term_label.split('_')[1]) + 1
        header.append(f'gamma_T{term_number}_Jm2')

    # Build data array
    data = np.zeros((len(delta_mu_O), len(term_labels) + 1))
    data[:, 0] = delta_mu_O
    for i, term_label in enumerate(term_labels):
        data[:, i + 1] = all_term_data[term_label]['gamma_array']

    # Save
    np.savetxt(output_path, data, delimiter=',', header=','.join(header), comments='')
    print(f"  Exported data to: {output_path}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    """Main execution: load data and generate plots."""

    print("\n" + "=" * 70)
    print(f"SURFACE THERMODYNAMICS ANALYSIS: {MATERIAL_NAME} {SURFACE_MILLER}")
    print("=" * 70)

    print(f"\nLoading AiiDA profile: {AIIDA_PROFILE}")
    load_profile(profile=AIIDA_PROFILE)
    print("  Profile loaded successfully")

    print(f"\nLoading workflow data (PK: {WORKFLOW_PK})...")
    try:
        surf_energies = load_surface_energies(WORKFLOW_PK)
        n_terms = len(surf_energies)
        print(f"  Found {n_terms} termination(s): {list(surf_energies.keys())}")

        # Print summary
        print_summary_table(surf_energies)

        # Generate plots
        safe_name = MATERIAL_NAME.lower().replace(' ', '_')
        miller_safe = SURFACE_MILLER.replace('(', '').replace(')', '')

        print("Generating plots...")

        # 1D gamma vs mu_O plot
        output_1d = os.path.join(OUTPUT_DIR, f'{safe_name}_{miller_safe}_gamma_vs_muO.{OUTPUT_FORMAT}')
        plot_gamma_vs_mu_O(MATERIAL_NAME, SURFACE_MILLER, surf_energies, output_1d)

        # Stability diagram
        output_stab = os.path.join(OUTPUT_DIR, f'{safe_name}_{miller_safe}_stability.{OUTPUT_FORMAT}')
        plot_stability_diagram(MATERIAL_NAME, SURFACE_MILLER, surf_energies, output_stab)

        # Export data to CSV
        output_csv = os.path.join(OUTPUT_DIR, f'{safe_name}_{miller_safe}_surface_energies.csv')
        export_data_to_csv(surf_energies, output_csv)

    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1

    print("\n" + "=" * 70)
    print("Analysis complete!")
    print("=" * 70 + "\n")

    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())

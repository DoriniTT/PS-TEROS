#!/usr/bin/env python
"""
Plot the planar average of the Fukui function in the z-direction.

Reads the pre-computed planar_average output directly from the Fukui WorkGraph.
"""

from aiida import orm, load_profile
import matplotlib.pyplot as plt
import numpy as np

# Load AiiDA profile
load_profile('presto')

# Load the Fukui WorkGraph
#WORKGRAPH_PK = 19982 # Fukui plus
WORKGRAPH_PK = 20163 # Fukui minus
print(f"Loading Fukui WorkGraph (PK {WORKGRAPH_PK})...")
wg_node = orm.load_node(WORKGRAPH_PK)

# Get the planar_average Dict output
planar_avg_dict = wg_node.outputs.planar_average.get_dict()

# Extract data
z_coords = np.array(planar_avg_dict['z_coordinates'])
z_avg = np.array(planar_avg_dict['planar_average'])
structure_info = planar_avg_dict['structure_info']
grid_dims = planar_avg_dict['grid_dimensions']

# Print info
print(f"Structure: {structure_info['formula']} ({structure_info['n_atoms']} atoms)")
print(f"Cell dimensions: a={structure_info['lattice_a']:.3f}, "
      f"b={structure_info['lattice_b']:.3f}, c={structure_info['lattice_c']:.3f} Å")
print(f"Grid dimensions: {grid_dims}")
print(f"Z-range: {planar_avg_dict['z_range']['min']:.3f} to {planar_avg_dict['z_range']['max']:.3f} Å")
print(f"Fukui function range: {planar_avg_dict['value_range']['min']:.6f} to {planar_avg_dict['value_range']['max']:.6f}")

# Create the plot
fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(z_coords, z_avg, 'b-', linewidth=1.5)
ax.set_xlabel('z (Å)', fontsize=12)
ax.set_ylabel('Fukui function f⁺(r) (e/Å³)', fontsize=12)
ax.set_title(f'Planar-averaged Fukui⁺ function - {structure_info["formula"]} (110)', fontsize=14)
ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
ax.grid(True, alpha=0.3)

# Save the figure
output_file = 'fukui_planar_average.png'
plt.tight_layout()
plt.savefig(output_file, dpi=150)
print(f"\nPlot saved to: {output_file}")

# Print summary
print(f"\n=== Summary ===")
print(f"Integral of |f⁺(z)|: {planar_avg_dict['integral']:.6f}")
print(f"Number of grid points: {planar_avg_dict['n_points']}")

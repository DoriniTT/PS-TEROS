"""
Plotting Utilities Submodule for TEROS.

This submodule provides functions for generating various types of plots
relevant to surface thermodynamics and materials analysis, such as
surface energy diagrams, phase diagrams, and other visualizations
derived from TEROS workflow outputs.
"""

# Import available plotting functions to make them accessible under this namespace
# and to control 'from teros.utils.plots import *' behavior.
try:
    from teros.utils.plots.delta_oxygen_chem_vs_surface_energy import (
        plot_delta_oxygen_chem_vs_surface_energy
    )
    __all__ = ['plot_delta_oxygen_chem_vs_surface_energy']
except ImportError:
    __all__ = []
"""
Teros: Ab Initio Atomistic Thermodynamics for Surface Energetics

A package for calculating surface Gibbs free energy using ab initio atomistic 
thermodynamics with AiiDA Workgraph framework.

Main Components
--------------
- core: Core workgraph definitions
- functions: Thermodynamics and structure manipulation functions
- utils: Utility functions and plotting tools
- examples: Example workflows for various DFT codes

References
---------
For theoretical background on ab initio atomistic thermodynamics, see:

1. M. Scheffler, C. Stampfl, "Theory of Adsorption on Metal Substrates",
   in: Handbook of Surface Science, Vol. 2: Electronic Structure, 
   Elsevier, Amsterdam (2000)

2. K. Reuter, M. Scheffler, "Composition, structure, and stability of 
   RuO2(110) as a function of oxygen pressure", 
   Phys. Rev. B 65, 035406 (2001)
"""

__version__ = '0.1.0'

# Import main components for easier access
try:
    from .core.workgraph import create_teros_workgraph
except ImportError:
    # Fallback for development/testing
    pass
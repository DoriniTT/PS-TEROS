"""
Default builders for PS-TEROS materials.

This module provides default VASP calculation parameters, scheduler options,
and structure configurations for specific material systems. These defaults
can be easily overridden or extended with custom parameters.

Example usage:
    >>> from teros.default_builders import get_ag3po4_defaults, update_builder_params
    >>> 
    >>> # Get defaults with required parameters
    >>> defaults = get_ag3po4_defaults(
    ...     structures_dir="/path/to/structures",
    ...     code_label="VASP-VTST-6.4.3@bohr",
    ...     potential_family="PBE"
    ... )
    >>> 
    >>> # Override specific parameters
    >>> defaults = update_builder_params(defaults, {
    ...     'bulk_parameters': {'ENCUT': 600},  # Increase ENCUT
    ...     'slab_options': {
    ...         'resources': {'num_machines': 2}  # Use 2 machines
    ...     }
    ... })
    >>> 
    >>> # Use with build_core_workgraph
    >>> wg = build_core_workgraph(**defaults, input_slabs=my_slabs)
"""

from copy import deepcopy

def update_builder_params(defaults, overrides):
    """
    Deep merge override parameters into default parameters.
    
    This function recursively merges dictionaries, updating nested dictionaries
    rather than replacing them entirely. This allows selective parameter updates.
    
    Args:
        defaults (dict): Default parameters dictionary
        overrides (dict): Override parameters to merge into defaults
        
    Returns:
        dict: Merged parameters dictionary
        
    Example:
        >>> defaults = {'params': {'ENCUT': 520, 'EDIFF': 1e-6}}
        >>> overrides = {'params': {'ENCUT': 600}}
        >>> result = update_builder_params(defaults, overrides)
        >>> result
        {'params': {'ENCUT': 600, 'EDIFF': 1e-6}}
    """
    result = deepcopy(defaults)
    
    for key, value in overrides.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            # Recursively merge nested dictionaries
            result[key] = update_builder_params(result[key], value)
        else:
            # Replace value for non-dict or new keys
            result[key] = deepcopy(value)
    
    return result


def get_ag3po4_defaults(
    structures_dir=None,
    code_label=None,
    potential_family=None,
    **overrides
):
    """
    Get default builder parameters for Ag3PO4 material system.
    
    This function returns a complete set of default VASP calculation parameters,
    scheduler options, structure information, slab generation parameters, and 
    potential mappings for the Ag3PO4 material system, including:
    - Bulk Ag3PO4 relaxation parameters
    - Metal (Ag) reference parameters
    - Nonmetal (P) reference parameters
    - Oxygen (O2) reference parameters
    - Slab generation parameters (Miller index, thickness, vacuum)
    - Slab relaxation parameters
    
    Note: Currently only ONE Miller index can be specified at a time.
    
    Args:
        structures_dir (str, optional): Path to directory containing structure files.
            If not provided, must be set before using with build_core_workgraph.
        code_label (str, optional): VASP code label (e.g., "VASP-VTST-6.4.3@bohr").
            If not provided, must be set before using with build_core_workgraph.
        potential_family (str, optional): Potential family name (e.g., "PBE").
            If not provided, must be set before using with build_core_workgraph.
        **overrides: Any parameter to override. Will be deep-merged with defaults.
            Example: bulk_parameters={'ENCUT': 600}, miller_indices=[1,1,0]
            
    Returns:
        dict: Complete builder parameters dictionary that can be unpacked with **
            into build_core_workgraph().

    Example:
        >>> from teros.core.builders.default_ag3po4_builders import get_ag3po4_defaults
        >>> from teros.core.workgraph import build_core_workgraph
        >>>
        >>> # Get defaults
        >>> defaults = get_ag3po4_defaults(
        ...     structures_dir="/home/user/structures",
        ...     code_label="VASP-VTST-6.4.3@bohr",
        ...     potential_family="PBE",
        ...     # Override specific parameters
        ...     bulk_parameters={'ENCUT': 600},
        ...     miller_indices=[1, 1, 0]  # (110) surface
        ... )
        >>>
        >>> # Use with workgraph builder (slab will be generated automatically)
        >>> wg = build_core_workgraph(
        ...     **defaults,
        ...     name="MyAg3PO4Calculation"
        ... )
        >>>
        >>> # Or provide your own slabs
        >>> wg = build_core_workgraph(
        ...     **defaults,
        ...     input_slabs=my_input_slabs,  # Overrides slab generation
        ...     name="MyAg3PO4Calculation"
        ... )
    """
    # ===== BULK RELAXATION PARAMETERS =====
    bulk_parameters = {
        "PREC": "Accurate",
        "ENCUT": 520,
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

    # ===== METAL (Ag) REFERENCE PARAMETERS =====
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

    # ===== NONMETAL (P) REFERENCE PARAMETERS =====
    nonmetal_parameters = {
        "PREC": "Accurate",
        "ENCUT": 520,
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

    nonmetal_options = {
        "resources": {
            "num_machines": 1,
            "num_cores_per_machine": 40,
        },
        "queue_name": "par40",
    }

    # ===== OXYGEN (O2) REFERENCE PARAMETERS =====
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

    # ===== SLAB RELAXATION PARAMETERS =====
    slab_parameters = {
        "PREC": "Accurate",
        "ENCUT": 520,
        "EDIFF": 1e-6,
        "ISMEAR": 0,
        "SIGMA": 0.05,
        "IBRION": 2,
        "ISIF": 2,
        "NSW": 100,
        "EDIFFG": -0.1,
        "ALGO": "Normal",
        "LREAL": "Auto",
        "LWAVE": False,
        "LCHARG": False,
        "LASPH": True,
    }

    slab_options = {
        "resources": {
            "num_machines": 1,
            "num_cores_per_machine": 40,
        },
        "queue_name": "par40",
    }

    # ===== POTENTIAL MAPPINGS =====
    bulk_potential_mapping = {"Ag": "Ag", "P": "P", "O": "O"}
    metal_potential_mapping = {"Ag": "Ag"}
    nonmetal_potential_mapping = {"P": "P"}
    oxygen_potential_mapping = {"O": "O"}
    slab_potential_mapping = {"Ag": "Ag", "P": "P", "O": "O"}

    # ===== STRUCTURE FILES =====
    bulk_name = "ag3po4.cif"
    metal_name = "Ag.cif"
    nonmetal_name = "P.cif"
    oxygen_name = "O2.cif"

    # ===== K-POINTS SETTINGS =====
    kpoints_spacing = 0.3
    slab_kpoints_spacing = 0.3

    # ===== SLAB GENERATION PARAMETERS =====
    # Default Miller index for Ag3PO4 (100) surface
    # Note: Currently only ONE surface can be generated at a time
    # Format: [h, k, l] not [(h, k, l)]
    miller_indices = [1, 0, 0]  # (100) surface
    min_slab_thickness = 18.0  # Minimum slab thickness in Angstroms
    min_vacuum_thickness = 15.0  # Minimum vacuum thickness in Angstroms

    # ===== OTHER SETTINGS =====
    clean_workdir = False
    relax_slabs = True
    compute_thermodynamics = True
    thermodynamics_sampling = 300

    # Build defaults dictionary
    defaults = {
        "structures_dir": structures_dir,
        "bulk_name": bulk_name,
        "metal_name": metal_name,
        "nonmetal_name": nonmetal_name,
        "oxygen_name": oxygen_name,
        "code_label": code_label,
        "potential_family": potential_family,
        "bulk_potential_mapping": bulk_potential_mapping,
        "metal_potential_mapping": metal_potential_mapping,
        "nonmetal_potential_mapping": nonmetal_potential_mapping,
        "oxygen_potential_mapping": oxygen_potential_mapping,
        "kpoints_spacing": kpoints_spacing,
        "bulk_parameters": bulk_parameters,
        "bulk_options": bulk_options,
        "metal_parameters": metal_parameters,
        "metal_options": metal_options,
        "nonmetal_parameters": nonmetal_parameters,
        "nonmetal_options": nonmetal_options,
        "oxygen_parameters": oxygen_parameters,
        "oxygen_options": oxygen_options,
        "clean_workdir": clean_workdir,
        # Slab generation parameters
        "miller_indices": miller_indices,
        "min_slab_thickness": min_slab_thickness,
        "min_vacuum_thickness": min_vacuum_thickness,
        # Slab relaxation
        "relax_slabs": relax_slabs,
        "slab_parameters": slab_parameters,
        "slab_options": slab_options,
        "slab_potential_mapping": slab_potential_mapping,
        "slab_kpoints_spacing": slab_kpoints_spacing,
        # Thermodynamics
        "compute_thermodynamics": compute_thermodynamics,
        "thermodynamics_sampling": thermodynamics_sampling,
    }

    # Apply overrides if provided
    if overrides:
        defaults = update_builder_params(defaults, overrides)

    return defaults

# Template for adding new materials:
# def get_<material>_defaults(structures_dir=None, code_label=None, potential_family=None, **overrides):
#     """
#     Get default builder parameters for <Material> material system.
#     
#     Args:
#         structures_dir (str, optional): Path to directory containing structure files.
#         code_label (str, optional): VASP code label.
#         potential_family (str, optional): Potential family name.
#         **overrides: Any parameter to override.
#         
#     Returns:
#         dict: Complete builder parameters dictionary.
#     """
#     # Define parameters similar to get_ag3po4_defaults()
#     # ...
#     pass

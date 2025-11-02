"""Utility functions for the serial surface thermodynamics preset."""

from aiida import orm


def prepare_vasp_parameters(
    base_parameters: dict,
    code: orm.Code,
    potential_family: str,
    potential_mapping: dict,
    kpoints_spacing: float,
    options: dict = None,
    clean_workdir: bool = False,
) -> dict:
    """
    Prepare standardized parameters for VASP WorkChain.

    Args:
        base_parameters: INCAR-like parameters dict
        code: VASP code to use
        potential_family: Potential family name
        potential_mapping: Element to potential mapping
        kpoints_spacing: K-points spacing
        options: Computer options (num_machines, etc.)
        clean_workdir: Whether to clean working directory

    Returns:
        Dictionary ready for VaspWorkChain inputs
    """
    params = {
        'code': code,
        'parameters': orm.Dict(dict=base_parameters),
        'potential_family': orm.Str(potential_family),
        'potential_mapping': orm.Dict(dict=potential_mapping),
        'options': orm.Dict(dict=options or {}),
        'kpoints_spacing': orm.Float(kpoints_spacing),
        'clean_workdir': orm.Bool(clean_workdir),
        'settings': orm.Dict(dict={
            'parser_settings': {
                'add_energy': True,
                'add_trajectory': True,
                'add_structure': True,
                'add_kpoints': True,
            }
        }),
    }
    return params


def create_default_bulk_parameters() -> dict:
    """Create default VASP parameters for bulk relaxation."""
    return {
        'PREC': 'Accurate',
        'EDIFF': 1e-6,
        'ENCUT': 520,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'IBRION': 2,
        'ISIF': 3,
        'NSW': 100,
        'LWAVE': False,
        'LCHARG': False,
    }


def create_default_slab_parameters() -> dict:
    """Create default VASP parameters for slab calculations."""
    return {
        'PREC': 'Accurate',
        'EDIFF': 1e-6,
        'ENCUT': 520,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'IBRION': 2,
        'ISIF': 2,  # Relax ions only, not cell
        'NSW': 100,
        'LWAVE': False,
        'LCHARG': False,
    }


def create_default_scf_parameters() -> dict:
    """Create default VASP parameters for SCF (single-point) calculations."""
    return {
        'PREC': 'Accurate',
        'EDIFF': 1e-6,
        'ENCUT': 520,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'NSW': 0,  # No ionic relaxation
        'LWAVE': False,
        'LCHARG': False,
    }

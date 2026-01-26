"""Explorer module for lightweight, incremental VASP calculations.

This module provides a simple API for exploratory VASP work:
- Submit a calculation, check results, decide next step, optionally restart
- No presets - always specify INCAR manually for maximum flexibility
- Specific file retrieval - get exactly the files you need
- Non-blocking default - submit and return immediately

Example usage:

    >>> from teros.core.explorer import quick_vasp, get_results, get_status
    >>>
    >>> # Single calculation
    >>> pk = quick_vasp(
    ...     structure=my_structure,
    ...     code_label='VASP-6.5.1@localwork',
    ...     incar={'NSW': 100, 'IBRION': 2},
    ...     kpoints_spacing=0.03,
    ...     potential_family='PBE',
    ...     potential_mapping={'Sn': 'Sn_d', 'O': 'O'},
    ...     options={'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 8}},
    ...     retrieve=['CONTCAR', 'CHGCAR'],
    ...     name='sno2_relax',
    ... )
    >>>
    >>> # Check status
    >>> get_status(pk)  # -> 'waiting', 'running', 'finished', 'failed'
    >>>
    >>> # Get results when done
    >>> results = get_results(pk)
    >>> print(f"Energy: {results['energy']:.4f} eV")
    >>>
    >>> # Restart from previous calculation
    >>> pk2 = quick_vasp(
    ...     restart_from=pk,
    ...     code_label='VASP-6.5.1@localwork',
    ...     incar={'NSW': 0, 'NEDOS': 2000},
    ...     retrieve=['DOSCAR'],
    ...     name='sno2_dos',
    ... )

DOS calculation using BandsWorkChain:

    >>> from teros.core.explorer import quick_dos, get_dos_results
    >>>
    >>> # DOS calculation (SCF + DOS handled internally)
    >>> # Note: AiiDA-VASP requires lowercase INCAR keys
    >>> pk = quick_dos(
    ...     structure=my_structure,
    ...     code_label='VASP-6.5.1@localwork',
    ...     scf_incar={'encut': 400, 'ediff': 1e-6, 'ismear': 0},
    ...     dos_incar={'nedos': 2000, 'lorbit': 11, 'ismear': -5},
    ...     kpoints_spacing=0.03,
    ...     dos_kpoints_spacing=0.02,
    ...     potential_family='PBE',
    ...     potential_mapping={'Sn': 'Sn_d', 'O': 'O'},
    ...     options={'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 8}},
    ...     retrieve=['DOSCAR'],
    ... )
    >>>
    >>> # Get DOS results
    >>> results = get_dos_results(pk)
    >>> print(f"Energy: {results['energy']:.4f} eV")
"""

from .workgraph import (
    quick_vasp,
    quick_vasp_batch,
    quick_dos,
    get_batch_results_from_workgraph,
)
from .results import (
    get_results,
    get_energy,
    get_batch_results,
    get_batch_energies,
    print_results,
    get_dos_results,
    print_dos_results,
)
from .utils import (
    get_status,
    export_files,
    list_calculations,
    get_restart_info,
)


__all__ = [
    # Core functions
    'quick_vasp',
    'quick_vasp_batch',
    'quick_dos',
    # Result extraction
    'get_results',
    'get_energy',
    'get_batch_results',
    'get_batch_energies',
    'get_batch_results_from_workgraph',
    'print_results',
    'get_dos_results',
    'print_dos_results',
    # Utilities
    'get_status',
    'export_files',
    'list_calculations',
    'get_restart_info',
]

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

Batch DOS calculation (multiple structures in parallel):

    >>> from teros.core.explorer import quick_dos_batch, get_batch_dos_results
    >>>
    >>> # Compare DOS for different structures
    >>> result = quick_dos_batch(
    ...     structures={'pristine': s1, 'vacancy': s2, 'interstitial': s3},
    ...     code_label='VASP-6.5.1@localwork',
    ...     scf_incar={'encut': 400, 'ediff': 1e-6, 'ismear': 0},
    ...     dos_incar={'nedos': 2000, 'lorbit': 11, 'ismear': -5},
    ...     kpoints_spacing=0.03,
    ...     potential_family='PBE',
    ...     potential_mapping={'Sn': 'Sn_d', 'O': 'O'},
    ...     options={'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 8}},
    ...     max_concurrent_jobs=2,  # Run 2 DOS calcs in parallel
    ...     retrieve=['DOSCAR'],
    ... )
    >>> print(f"WorkGraph PK: {result['__workgraph_pk__']}")
    >>>
    >>> # Get batch results when done
    >>> batch_results = get_batch_dos_results(result)
    >>> for key, dos_result in batch_results.items():
    ...     print(f"{key}: E = {dos_result['energy']:.4f} eV")

Sequential multi-stage calculation with restart chaining:

    >>> from teros.core.explorer import quick_vasp_sequential, print_sequential_results
    >>>
    >>> # Define stages with automatic restart chaining
    >>> stages = [
    ...     {
    ...         'name': 'relax_rough',
    ...         'incar': {'NSW': 100, 'IBRION': 2, 'ISIF': 2, 'ENCUT': 400},
    ...         'kpoints_spacing': 0.06,
    ...         'retrieve': ['CONTCAR', 'OUTCAR'],
    ...     },
    ...     {
    ...         'name': 'relax_fine',
    ...         'incar': {'NSW': 100, 'IBRION': 2, 'ISIF': 2, 'ENCUT': 520},
    ...         'kpoints_spacing': 0.03,
    ...         'retrieve': ['CONTCAR', 'OUTCAR', 'WAVECAR'],
    ...         # restart='previous' is default -> restart from previous stage
    ...     },
    ...     {
    ...         'name': 'relax_supercell',
    ...         'supercell': [2, 2, 1],  # Create supercell, no restart_folder
    ...         'incar': {'NSW': 100, 'IBRION': 2, 'ISIF': 2, 'ENCUT': 520},
    ...         'kpoints_spacing': 0.03,
    ...         'retrieve': ['CONTCAR', 'OUTCAR', 'CHGCAR'],
    ...     },
    ... ]
    >>>
    >>> result = quick_vasp_sequential(
    ...     structure=my_structure,
    ...     stages=stages,
    ...     code_label='VASP-6.5.1@localwork',
    ...     potential_family='PBE',
    ...     potential_mapping={'Sn': 'Sn_d', 'O': 'O'},
    ...     options={'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 8}},
    ... )
    >>>
    >>> # Get results when done
    >>> print_sequential_results(result)
"""

from .workgraph import (
    quick_vasp,
    quick_vasp_batch,
    quick_vasp_sequential,
    quick_dos,
    quick_dos_batch,
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
    get_batch_dos_results,
    print_batch_dos_results,
    get_sequential_results,
    get_stage_results,
    print_sequential_results,
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
    'quick_vasp_sequential',
    'quick_dos',
    'quick_dos_batch',
    # Result extraction
    'get_results',
    'get_energy',
    'get_batch_results',
    'get_batch_energies',
    'get_batch_results_from_workgraph',
    'print_results',
    'get_dos_results',
    'print_dos_results',
    'get_batch_dos_results',
    'print_batch_dos_results',
    'get_sequential_results',
    'get_stage_results',
    'print_sequential_results',
    # Utilities
    'get_status',
    'export_files',
    'list_calculations',
    'get_restart_info',
]

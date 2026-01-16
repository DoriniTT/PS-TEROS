"""WorkGraph builder for Fukui function calculations.

This module implements parallel VASP calculations at different charge states
to generate CHGCAR files for Fukui function interpolation.
"""

import typing as t

from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import WorkGraph, task, namespace, dynamic

from .tasks import (
    extract_total_energy,
    collect_chgcar_files_internal,
    generate_fukui_summary_internal,
    run_fukui_interpolation_calcfunc,
)
from .utils import make_delta_label, validate_fukui_inputs, DEFAULT_DELTA_N_VALUES
from ..utils import deep_merge_dicts


def _prepare_fukui_inputs(
    builder_inputs: dict,
    nelect: float,
    structure: orm.StructureData,
    code: orm.InstalledCode,
) -> dict:
    """
    Prepare builder inputs for a single Fukui charge-state calculation.

    Converts plain dicts to AiiDA types and applies Fukui-specific overrides:
    - NSW=0 (static calculation)
    - IBRION=-1 (no ionic relaxation)
    - LCHARG=True (write CHGCAR)
    - NELECT=<value> (fractional electron count)
    - ADDITIONAL_RETRIEVE_LIST includes CHGCAR

    Args:
        builder_inputs: Base builder inputs dict
        nelect: NELECT value for this calculation
        structure: Input structure
        code: VASP code

    Returns:
        dict with AiiDA-compatible types ready for VaspTask
    """
    # Static calculation overrides
    static_overrides = {
        'parameters': {
            'incar': {
                'nsw': 0,
                'ibrion': -1,
                'lcharg': True,
                'lwave': False,
                'nelect': nelect,
            }
        },
    }

    # Merge overrides with user inputs (overrides take precedence)
    merged = deep_merge_dicts(builder_inputs, static_overrides)

    # Ensure CHGCAR is retrieved
    if 'settings' not in merged:
        merged['settings'] = {}
    if 'ADDITIONAL_RETRIEVE_LIST' not in merged['settings']:
        merged['settings']['ADDITIONAL_RETRIEVE_LIST'] = []
    if 'CHGCAR' not in merged['settings']['ADDITIONAL_RETRIEVE_LIST']:
        merged['settings']['ADDITIONAL_RETRIEVE_LIST'].append('CHGCAR')

    # Convert to AiiDA types
    prepared = {
        'structure': structure,
        'code': code,
    }

    # Convert parameters dict to orm.Dict
    if 'parameters' in merged:
        if isinstance(merged['parameters'], dict):
            prepared['parameters'] = orm.Dict(dict=merged['parameters'])
        else:
            prepared['parameters'] = merged['parameters']

    # Convert options dict to orm.Dict
    if 'options' in merged:
        if isinstance(merged['options'], dict):
            prepared['options'] = orm.Dict(dict=merged['options'])
        else:
            prepared['options'] = merged['options']

    # Convert potential_mapping dict to orm.Dict
    if 'potential_mapping' in merged:
        if isinstance(merged['potential_mapping'], dict):
            prepared['potential_mapping'] = orm.Dict(dict=merged['potential_mapping'])
        else:
            prepared['potential_mapping'] = merged['potential_mapping']

    # Convert settings dict to orm.Dict
    if 'settings' in merged:
        if isinstance(merged['settings'], dict):
            prepared['settings'] = orm.Dict(dict=merged['settings'])
        else:
            prepared['settings'] = merged['settings']

    # Handle kpoints_spacing - keep as plain Python float
    if 'kpoints_spacing' in merged:
        kps = merged['kpoints_spacing']
        if isinstance(kps, (int, float)):
            prepared['kpoints_spacing'] = float(kps)
        else:
            prepared['kpoints_spacing'] = kps

    # Copy string/bool values directly
    for key in ('potential_family', 'clean_workdir'):
        if key in merged:
            prepared[key] = merged[key]

    return prepared


@task.graph
def FukuiCalculationScatter(
    structure: orm.StructureData,
    nelect_neutral: orm.Int,
    nelect_values: orm.List,
    delta_n_values: orm.List,
    delta_n_labels: orm.List,
    fukui_type: orm.Str,
    code: orm.InstalledCode,
    builder_inputs: dict,
    max_number_jobs: int = None,
) -> t.Annotated[dict, namespace(
    chgcar_files=orm.FolderData,
    summary=orm.Dict,
)]:
    """
    Run parallel VASP static calculations at different NELECT values.

    This is the core scatter-gather task that creates independent VASP
    calculations for each charge state, collects CHGCAR files, and generates
    a summary.

    Args:
        structure: Input structure (same for all calculations)
        nelect_neutral: Number of electrons in neutral system
        nelect_values: List of NELECT values (e.g., [320.00, 319.95, 319.90, 319.85])
        delta_n_values: List of delta_n values (e.g., [0.0, 0.05, 0.10, 0.15])
        delta_n_labels: List of labels (e.g., ['delta_0_00', 'delta_0_05', ...])
        fukui_type: Type of Fukui calculation ('plus' or 'minus')
        code: VASP InstalledCode
        builder_inputs: Base builder inputs dict
        max_number_jobs: Maximum concurrent VASP calculations (optional)

    Returns:
        Dictionary with outputs:
            - chgcar_files: FolderData containing all CHGCAR files
            - summary: Dict with calculation metadata
    """
    from aiida_workgraph import get_current_graph

    # Set max_number_jobs on this workgraph to control concurrency
    if max_number_jobs is not None:
        wg = get_current_graph()
        wg.max_number_jobs = max_number_jobs

    # Get VASP workchain and wrap as task
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    VaspTask = task(VaspWorkChain)

    # Store outputs for collection
    retrieved_list = []
    misc_list = []

    nelect_list = nelect_values.get_list()
    delta_n_list = delta_n_values.get_list()
    labels = delta_n_labels.get_list()

    for i, label in enumerate(labels):
        nelect = nelect_list[i]

        # Prepare inputs with specific NELECT
        prepared_inputs = _prepare_fukui_inputs(
            builder_inputs, nelect, structure, code
        )

        # Create VASP task
        vasp_calc = VaspTask(**prepared_inputs)

        # Collect outputs (in order)
        retrieved_list.append(vasp_calc.retrieved)
        misc_list.append(vasp_calc.misc)

    # Collect CHGCAR files using internal calcfunction
    chgcar_result = collect_chgcar_files_internal(
        delta_n_list=delta_n_values,
        labels_list=delta_n_labels,
        retrieved_0=retrieved_list[0] if len(retrieved_list) > 0 else None,
        retrieved_1=retrieved_list[1] if len(retrieved_list) > 1 else None,
        retrieved_2=retrieved_list[2] if len(retrieved_list) > 2 else None,
        retrieved_3=retrieved_list[3] if len(retrieved_list) > 3 else None,
    )

    # Generate summary using internal calcfunction
    summary_result = generate_fukui_summary_internal(
        nelect_neutral=nelect_neutral,
        delta_n_list=delta_n_values,
        fukui_type=fukui_type,
        misc_0=misc_list[0] if len(misc_list) > 0 else None,
        misc_1=misc_list[1] if len(misc_list) > 1 else None,
        misc_2=misc_list[2] if len(misc_list) > 2 else None,
        misc_3=misc_list[3] if len(misc_list) > 3 else None,
    )

    return {
        'chgcar_files': chgcar_result.result,
        'summary': summary_result.result,
    }


def build_fukui_workgraph(
    structure: t.Union[orm.StructureData, int],
    nelect_neutral: int,
    delta_n_values: t.List[float] = None,
    code_label: str = None,
    builder_inputs: dict = None,
    relax_first: bool = False,
    relax_builder_inputs: dict = None,
    fix_atoms: bool = False,
    fix_type: str = 'center',
    fix_thickness: float = 0.0,
    fix_elements: t.List[str] = None,
    fukui_type: str = 'plus',
    compute_fukui: bool = False,
    max_concurrent_jobs: int = None,
    name: str = 'FukuiWorkGraph',
) -> WorkGraph:
    """
    Build a WorkGraph for Fukui function calculations.

    The Fukui function is calculated using finite differences on electron density:
    - f+(r) = [rho(N) - rho(N-dN)] / dN  (nucleophilic attack - where electrons are added)
    - f-(r) = [rho(N+dN) - rho(N)] / dN  (electrophilic attack - where electrons are removed)

    This workflow:
    1. (Optional) Relaxes the input structure
    2. Runs parallel static VASP calculations at different NELECT values
    3. Collects CHGCAR files into a single FolderData output

    Args:
        structure: Input structure (StructureData or PK)
        nelect_neutral: Number of electrons in the neutral system (required).
                       This can be determined from a previous VASP calculation
                       (check OUTCAR or vasprun.xml).
        delta_n_values: List of dN values for fractional charges.
                       Default: [0.0, 0.05, 0.10, 0.15]
                       NELECT = nelect_neutral - delta_n for Fukui+
                       NELECT = nelect_neutral + delta_n for Fukui-
        code_label: VASP code label (e.g., 'VASP-6.5.1@localwork')
        builder_inputs: VASP builder configuration dict.
                       NSW and IBRION are overridden for static calculations.
                       Format: {'parameters': {'incar': {...}}, 'options': {...}, ...}
        relax_first: If True, perform geometry relaxation before charge calculations.
                    Default: False
        relax_builder_inputs: Optional separate builder inputs for relaxation.
                             If None and relax_first=True, uses builder_inputs.
        fix_atoms: If True, apply selective dynamics during relaxation.
                  Only used when relax_first=True. Default: False
        fix_type: Where to fix atoms: 'center' (fix center, relax surfaces),
                 'bottom', or 'top'. Default: 'center'
        fix_thickness: Thickness in Angstroms for the fixing region.
                      For 'center', fixes atoms within ±fix_thickness/2 from center.
        fix_elements: Optional list of element symbols to fix (e.g., ['Sn']).
                     If None, all atoms in the region are fixed.
        fukui_type: 'plus' for nucleophilic (remove electrons),
                   'minus' for electrophilic (add electrons)
        compute_fukui: If True, run FukuiGrid interpolation to compute the
                      Fukui function from the CHGCAR files. Output will include
                      'fukui_chgcar' (SinglefileData containing CHGCAR_FUKUI.vasp).
                      Default: False
        max_concurrent_jobs: Limit parallel VASP calculations (default: unlimited)
        name: WorkGraph name

    Returns:
        WorkGraph ready for submission

    Outputs (available after completion):
        - chgcar_files: FolderData containing CHGCAR_0.00, CHGCAR_0.05, etc.
        - summary: Dict with calculation metadata (energies, NELECT values, convergence)
        - relaxed_structure: StructureData (only if relax_first=True)
        - fukui_chgcar: SinglefileData with CHGCAR_FUKUI.vasp (only if compute_fukui=True)

    Example:
        >>> wg = build_fukui_workgraph(
        ...     structure=my_structure,
        ...     nelect_neutral=192,
        ...     code_label='VASP-6.5.1@localwork',
        ...     builder_inputs={
        ...         'parameters': {'incar': {'encut': 400, 'ediff': 1e-5}},
        ...         'options': {'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 8}},
        ...         'kpoints_spacing': 0.3,
        ...         'potential_family': 'PBE',
        ...         'potential_mapping': {'Sn': 'Sn_d', 'O': 'O'},
        ...     },
        ...     fukui_type='plus',
        ... )
        >>> wg.submit(wait=False)
    """
    # Set defaults
    if delta_n_values is None:
        delta_n_values = DEFAULT_DELTA_N_VALUES.copy()

    if builder_inputs is None:
        builder_inputs = {}

    # Validate inputs
    validate_fukui_inputs(nelect_neutral, delta_n_values, fukui_type)

    # Load structure if PK provided
    if isinstance(structure, int):
        structure = orm.load_node(structure)

    # Load code
    code = orm.load_code(code_label)

    # Calculate NELECT values for each delta_n
    if fukui_type == 'plus':
        # Fukui+: remove electrons (nucleophilic attack)
        nelect_values = [nelect_neutral - dn for dn in delta_n_values]
    else:
        # Fukui-: add electrons (electrophilic attack)
        nelect_values = [nelect_neutral + dn for dn in delta_n_values]

    # Create labels for each calculation
    delta_n_labels = [make_delta_label(dn) for dn in delta_n_values]

    # Create WorkGraph
    wg = WorkGraph(name=name)

    # Track the structure to use for charge calculations
    structure_for_charges = structure

    # Optional relaxation step
    if relax_first:
        VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
        VaspTask = task(VaspWorkChain)

        # Use relax_builder_inputs if provided, otherwise use builder_inputs
        relax_inputs = relax_builder_inputs if relax_builder_inputs else builder_inputs

        # Prepare relaxation inputs (don't apply static overrides)
        prepared_relax = {}

        # Convert to AiiDA types
        if 'parameters' in relax_inputs:
            if isinstance(relax_inputs['parameters'], dict):
                prepared_relax['parameters'] = orm.Dict(dict=relax_inputs['parameters'])
            else:
                prepared_relax['parameters'] = relax_inputs['parameters']

        if 'options' in relax_inputs:
            if isinstance(relax_inputs['options'], dict):
                prepared_relax['options'] = orm.Dict(dict=relax_inputs['options'])
            else:
                prepared_relax['options'] = relax_inputs['options']

        if 'potential_mapping' in relax_inputs:
            if isinstance(relax_inputs['potential_mapping'], dict):
                prepared_relax['potential_mapping'] = orm.Dict(dict=relax_inputs['potential_mapping'])
            else:
                prepared_relax['potential_mapping'] = relax_inputs['potential_mapping']

        if 'settings' in relax_inputs:
            if isinstance(relax_inputs['settings'], dict):
                prepared_relax['settings'] = orm.Dict(dict=relax_inputs['settings'])
            else:
                prepared_relax['settings'] = relax_inputs['settings']

        if 'kpoints_spacing' in relax_inputs:
            prepared_relax['kpoints_spacing'] = float(relax_inputs['kpoints_spacing'])

        for key in ('potential_family', 'clean_workdir'):
            if key in relax_inputs:
                prepared_relax[key] = relax_inputs[key]

        # Apply selective dynamics if fix_atoms=True
        if fix_atoms and fix_type is not None and fix_thickness > 0.0:
            from ..fixed_atoms import get_fixed_atoms_list

            fixed_atoms_list = get_fixed_atoms_list(
                structure=structure,
                fix_type=fix_type,
                fix_thickness=fix_thickness,
                fix_elements=fix_elements,
            )

            if fixed_atoms_list:
                num_atoms = len(structure.sites)
                positions_dof = []
                for i in range(1, num_atoms + 1):  # 1-based indexing
                    if i in fixed_atoms_list:
                        positions_dof.append([False, False, False])  # Fixed
                    else:
                        positions_dof.append([True, True, True])     # Relax

                prepared_relax['dynamics'] = orm.Dict(dict={'positions_dof': positions_dof})

        # Add relaxation task
        relax_task = wg.add_task(
            VaspTask,
            name='relax_structure',
            structure=structure,
            code=code,
            **prepared_relax
        )

        # Use relaxed structure for charge calculations
        structure_for_charges = relax_task.outputs.structure

        # Expose relaxed structure as output
        wg.outputs.relaxed_structure = relax_task.outputs.structure

    # Add scatter task for charge calculations
    scatter_task = wg.add_task(
        FukuiCalculationScatter,
        name='fukui_scatter',
        structure=structure_for_charges,
        nelect_neutral=orm.Int(nelect_neutral),
        nelect_values=orm.List(list=nelect_values),
        delta_n_values=orm.List(list=delta_n_values),
        delta_n_labels=orm.List(list=delta_n_labels),
        fukui_type=orm.Str(fukui_type),
        code=code,
        builder_inputs=builder_inputs,
    )

    # Expose outputs from scatter task
    wg.outputs.chgcar_files = scatter_task.outputs.chgcar_files
    wg.outputs.summary = scatter_task.outputs.summary

    # Optional: Run FukuiGrid interpolation to compute Fukui function
    if compute_fukui:
        interpolation_task = wg.add_task(
            run_fukui_interpolation_calcfunc,
            name='fukui_interpolation',
            chgcar_files=scatter_task.outputs.chgcar_files,
            delta_n_values=orm.List(list=delta_n_values),
            fukui_type=orm.Str(fukui_type),
        )
        wg.outputs.fukui_chgcar = interpolation_task.outputs.result

    # Set max_concurrent_jobs at WorkGraph level
    if max_concurrent_jobs is not None:
        wg.max_number_jobs = max_concurrent_jobs

    return wg


def get_fukui_results(workgraph) -> dict:
    """
    Extract results from completed Fukui WorkGraph.

    Args:
        workgraph: Completed WorkGraph node (PK, UUID, or WorkGraphNode)

    Returns:
        dict with:
        - chgcar_folder: FolderData containing all CHGCAR files
        - summary: Dict with calculation metadata
        - file_names: List of CHGCAR file names in the folder
        - relaxed_structure: StructureData (if relax_first was True)
        - fukui_chgcar: SinglefileData with CHGCAR_FUKUI.vasp (if compute_fukui was True)

    Example:
        >>> results = get_fukui_results(wg.pk)
        >>> print(results['file_names'])
        ['CHGCAR_0.00', 'CHGCAR_0.05', 'CHGCAR_0.10', 'CHGCAR_0.15']
        >>> print(results['summary']['fukui_type'])
        'plus'
    """
    # Load workgraph if PK or UUID provided
    if isinstance(workgraph, (int, str)):
        workgraph = orm.load_node(workgraph)

    results = {}

    # Get chgcar_files output
    if hasattr(workgraph.outputs, 'chgcar_files'):
        chgcar_folder = workgraph.outputs.chgcar_files
        results['chgcar_folder'] = chgcar_folder
        results['file_names'] = sorted(chgcar_folder.list_object_names())
    else:
        results['chgcar_folder'] = None
        results['file_names'] = []

    # Get summary output
    if hasattr(workgraph.outputs, 'summary'):
        summary = workgraph.outputs.summary
        if isinstance(summary, orm.Dict):
            results['summary'] = summary.get_dict()
        else:
            results['summary'] = dict(summary)
    else:
        results['summary'] = {}

    # Get relaxed structure if available
    if hasattr(workgraph.outputs, 'relaxed_structure'):
        results['relaxed_structure'] = workgraph.outputs.relaxed_structure
    else:
        results['relaxed_structure'] = None

    # Get fukui_chgcar if available (from compute_fukui=True)
    if hasattr(workgraph.outputs, 'fukui_chgcar'):
        results['fukui_chgcar'] = workgraph.outputs.fukui_chgcar
    else:
        results['fukui_chgcar'] = None

    return results


def print_fukui_summary(workgraph) -> None:
    """
    Print formatted summary of Fukui calculation results.

    Args:
        workgraph: Completed WorkGraph node (PK, UUID, or WorkGraphNode)
    """
    results = get_fukui_results(workgraph)
    summary = results.get('summary', {})

    print("\n" + "=" * 60)
    print("FUKUI CALCULATION SUMMARY")
    print("=" * 60)

    if summary:
        print(f"\nFukui Type: {summary.get('fukui_type', 'N/A')}")
        print(f"NELECT (neutral): {summary.get('nelect_neutral', 'N/A')}")
        print(f"Delta N values: {summary.get('delta_n_values', [])}")
        print(f"Number of calculations: {summary.get('n_calculations', 'N/A')}")

        print(f"\nFormula: {summary.get('formula', 'N/A')}")

        calculations = summary.get('calculations', {})
        if calculations:
            print("\nCalculations:")
            print("-" * 50)
            print(f"{'Label':<15} {'Delta N':<10} {'NELECT':<12} {'Energy (eV)':<15}")
            print("-" * 50)
            for label, data in sorted(calculations.items()):
                delta_n = data.get('delta_n', 'N/A')
                nelect = data.get('nelect', 'N/A')
                energy = data.get('energy_eV', 'N/A')
                if energy is not None and energy != 'N/A':
                    energy = f"{energy:.6f}"
                print(f"{label:<15} {delta_n:<10} {nelect:<12} {energy:<15}")

    print("\nCHGCAR Files:")
    print("-" * 50)
    for fname in results.get('file_names', []):
        print(f"  - {fname}")

    if results.get('relaxed_structure'):
        print(f"\nRelaxed Structure: PK={results['relaxed_structure'].pk}")

    if results.get('fukui_chgcar'):
        fukui_file = results['fukui_chgcar']
        print(f"\nFukui Function: PK={fukui_file.pk}")
        print(f"  Filename: {fukui_file.filename}")
        # Get file size if possible
        try:
            content = fukui_file.get_content()
            size_mb = len(content) / (1024 * 1024)
            print(f"  Size: {size_mb:.1f} MB")
        except Exception:
            pass

    print("=" * 60)

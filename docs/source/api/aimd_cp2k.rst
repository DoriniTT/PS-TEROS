===============================
AIMD CP2K Module
===============================

.. automodule:: teros.core.aimd_cp2k
   :members:
   :undoc-members:
   :show-inheritance:

Key Functions
=============

aimd_single_stage_scatter_cp2k()
---------------------------------

.. autofunction:: teros.core.aimd_cp2k.aimd_single_stage_scatter_cp2k

Description
^^^^^^^^^^^

Run ab initio molecular dynamics on multiple slab structures in parallel using CP2K.

This function is the CP2K equivalent of ``aimd_single_stage_scatter`` from the VASP module.
It handles ONE AIMD stage (temperature/timestep configuration) for all slabs simultaneously.

**Key Differences from VASP:**

* Uses CP2K's GPW method instead of PAW
* Requires BASIS_MOLOPT and GTH_POTENTIALS files
* Different restart mechanism (``aiida-1.restart`` instead of WAVECAR)
* Outputs include CP2K-specific trajectory format

Parameters
^^^^^^^^^^

:param slabs: Dictionary mapping slab labels to StructureData nodes
:type slabs: dict[str, orm.StructureData]

:param temperature: Target temperature in Kelvin
:type temperature: float

:param steps: Number of MD steps for this stage
:type steps: int

:param code: CP2K code configured in AiiDA
:type code: orm.Code

:param aimd_parameters: CP2K input parameters (from ``get_aimd_defaults_cp2k``)
:type aimd_parameters: dict

:param basis_file: SinglefileData containing BASIS_MOLOPT
:type basis_file: orm.SinglefileData

:param pseudo_file: SinglefileData containing GTH_POTENTIALS
:type pseudo_file: orm.SinglefileData

:param options: Scheduler options for CP2K calculations
:type options: dict

:param clean_workdir: Whether to clean work directory after completion
:type clean_workdir: bool

:param restart_folders: Remote folders from previous stage (for restart)
:type restart_folders: dict[str, orm.RemoteData], optional

:param fixed_atoms_lists: Pre-computed fixed atom indices for each slab
:type fixed_atoms_lists: dict[str, list[int]], optional

:param fix_components: Components to fix ('XYZ', 'XY', 'Z')
:type fix_components: str, optional

:param fix_type: Fixing type for dynamic calculation ('bottom', 'top', 'center')
:type fix_type: str, optional

:param fix_thickness: Thickness of fixed region in Angstroms
:type fix_thickness: float, optional

:param fix_elements: List of element symbols to fix (None = all)
:type fix_elements: list[str], optional

Returns
^^^^^^^

Dictionary with namespaced outputs for each slab:

* ``structures``: Final structures from AIMD (dynamic namespace)
* ``remote_folders``: RemoteData for next stage restart (dynamic namespace)
* ``parameters``: Output parameters from CP2K (dynamic namespace)
* ``trajectories``: Trajectory data from MD (dynamic namespace)
* ``retrieved``: Retrieved folder data (dynamic namespace)

Example
^^^^^^^

Basic usage:

.. code-block:: python

    from teros.core.aimd_cp2k import aimd_single_stage_scatter_cp2k
    from teros.core.builders.aimd_builder_cp2k import get_aimd_defaults_cp2k
    from aiida import orm
    import io

    # Prepare CP2K parameters
    aimd_params = get_aimd_defaults_cp2k(
        cutoff=400,
        rel_cutoff=60,
        timestep=1.0,
        thermostat='NOSE',
    )

    # Add KIND section
    aimd_params['FORCE_EVAL']['SUBSYS']['KIND'] = [
        {"_": "Ag", "BASIS_SET": "DZVP-MOLOPT-PBE-GTH-q11", "POTENTIAL": "GTH-PBE-q11"},
        {"_": "O", "BASIS_SET": "DZVP-MOLOPT-PBE-GTH-q6", "POTENTIAL": "GTH-PBE-q6"},
    ]

    # Prepare basis and pseudo files
    from teros.core.builders.aimd_builder_cp2k import (
        get_basis_molopt_content,
        get_gth_potentials_content,
    )

    basis_file = orm.SinglefileData(
        io.BytesIO(get_basis_molopt_content().encode('utf-8')),
        filename='BASIS_MOLOPT'
    )
    pseudo_file = orm.SinglefileData(
        io.BytesIO(get_gth_potentials_content().encode('utf-8')),
        filename='GTH_POTENTIALS'
    )

    # Load slabs and code
    slabs = {
        'slab_111': orm.load_node(<STRUCTURE_PK>),
    }
    code = orm.load_code('CP2K@computer')

    # Run AIMD stage
    stage = aimd_single_stage_scatter_cp2k(
        slabs=slabs,
        temperature=300.0,
        steps=100,
        code=code,
        aimd_parameters=aimd_params,
        basis_file=basis_file,
        pseudo_file=pseudo_file,
        options={'resources': {'num_machines': 1, 'num_cores_per_machine': 48}},
        clean_workdir=False,
    )

With fixed atoms:

.. code-block:: python

    # Run AIMD with bottom 7 Å fixed
    stage = aimd_single_stage_scatter_cp2k(
        slabs=slabs,
        temperature=300.0,
        steps=100,
        code=code,
        aimd_parameters=aimd_params,
        basis_file=basis_file,
        pseudo_file=pseudo_file,
        options=options,
        clean_workdir=False,
        # Fixed atoms (calculated dynamically for each slab)
        fix_type='bottom',
        fix_thickness=7.0,
        fix_elements=None,  # All elements
        fix_components='XYZ',  # Fully rigid
    )

Notes
^^^^^

* Temperature and steps in ``aimd_parameters`` are automatically overridden
* BASIS_MOLOPT and GTH_POTENTIALS files are auto-generated by PS-TEROS
* For sequential stages, pass ``restart_folders`` from previous stage
* Fixed atoms can be pre-computed (``fixed_atoms_lists``) or calculated dynamically (``fix_type``, ``fix_thickness``)
* CP2K constraint section is automatically added when fixed atoms are specified

See Also
^^^^^^^^

* :func:`~teros.core.aimd.aimd_single_stage_scatter` - VASP equivalent
* :func:`~teros.core.builders.aimd_builder_cp2k.get_aimd_defaults_cp2k` - Default parameters
* :func:`~teros.core.fixed_atoms.get_fixed_atoms_list` - Calculate fixed atoms
* :func:`~teros.core.fixed_atoms.add_fixed_atoms_to_cp2k_parameters` - Add constraints

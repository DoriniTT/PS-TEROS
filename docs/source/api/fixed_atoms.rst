===============================
Fixed Atoms Module
===============================

.. automodule:: teros.core.fixed_atoms
   :members:
   :undoc-members:
   :show-inheritance:

Overview
========

The ``fixed_atoms`` module provides calculator-agnostic utilities for constraining atoms in slab structures during AIMD simulations. This is essential for realistic surface modeling where bottom layers should remain fixed to represent the bulk substrate.

Key Functions
=============

get_fixed_atoms_list()
-----------------------

.. autofunction:: teros.core.fixed_atoms.get_fixed_atoms_list

Description
^^^^^^^^^^^

Identify atoms to fix in a slab structure based on position criteria.

This function analyzes a slab structure and returns a list of atom indices that should be constrained during MD simulations. Atoms are selected based on their Z-coordinate (perpendicular to the slab surface).

Parameters
^^^^^^^^^^

:param structure: Slab structure to analyze
:type structure: orm.StructureData

:param fix_type: Where to fix atoms - 'bottom', 'top', 'center', or None
:type fix_type: str, optional

:param fix_thickness: Thickness in Angstroms for fixing region
:type fix_thickness: float, optional

:param fix_elements: Optional list of element symbols to fix (e.g., ['Ag', 'O'])
:type fix_elements: list[str], optional

Returns
^^^^^^^

List of 1-based atom indices to fix (sorted)

:rtype: list[int]

Examples
^^^^^^^^

Fix bottom 7 Å of all atoms:

.. code-block:: python

    from teros.core.fixed_atoms import get_fixed_atoms_list
    from aiida import orm

    slab = orm.load_node(<STRUCTURE_PK>)
    fixed = get_fixed_atoms_list(slab, fix_type='bottom', fix_thickness=7.0)
    print(f"Fixed atoms: {fixed}")
    # Output: [1, 2, 3, 4, 5, 6, 7, 8, 21, 22, 23, 24, 26]

Fix top 5 Å of Ag atoms only:

.. code-block:: python

    fixed = get_fixed_atoms_list(
        slab,
        fix_type='top',
        fix_thickness=5.0,
        fix_elements=['Ag']
    )

Fix 4 Å around center (2 Å above and below):

.. code-block:: python

    fixed = get_fixed_atoms_list(slab, fix_type='center', fix_thickness=4.0)

Notes
^^^^^

* Returns empty list if ``fix_type=None`` or ``fix_thickness<=0``
* Indices are 1-based (FORTRAN convention) for direct use in VASP POSCAR or CP2K input
* Z-coordinate is assumed to be perpendicular to the slab surface
* For ``fix_type='center'``, the region is ``fix_thickness/2`` above and below the slab center


add_fixed_atoms_to_cp2k_parameters()
--------------------------------------

.. autofunction:: teros.core.fixed_atoms.add_fixed_atoms_to_cp2k_parameters

Description
^^^^^^^^^^^

Add FIXED_ATOMS constraint section to CP2K input parameters.

This function modifies CP2K input parameters dictionary to include atomic constraints. It creates the necessary ``MOTION/CONSTRAINT/FIXED_ATOMS`` section with the specified atom list and components.

Parameters
^^^^^^^^^^

:param parameters: CP2K input parameters dictionary
:type parameters: dict

:param fixed_atoms_list: List of 1-based atom indices to fix
:type fixed_atoms_list: list[int]

:param fix_components: Components to fix - 'XYZ', 'XY', or 'Z'
:type fix_components: str, default='XYZ'

Returns
^^^^^^^

Modified parameters dictionary with constraint section

:rtype: dict

Example
^^^^^^^

.. code-block:: python

    from teros.core.fixed_atoms import (
        get_fixed_atoms_list,
        add_fixed_atoms_to_cp2k_parameters,
    )
    from teros.core.builders.aimd_builder_cp2k import get_aimd_defaults_cp2k

    # Get base parameters
    params = get_aimd_defaults_cp2k(cutoff=400, timestep=1.0)

    # Identify fixed atoms
    slab = orm.load_node(<STRUCTURE_PK>)
    fixed_list = get_fixed_atoms_list(slab, fix_type='bottom', fix_thickness=7.0)

    # Add constraints
    params = add_fixed_atoms_to_cp2k_parameters(
        params,
        fixed_list,
        fix_components='XYZ'  # Fully rigid
    )

    # Result: params now contains:
    # params['MOTION']['CONSTRAINT']['FIXED_ATOMS'] = {
    #     'COMPONENTS_TO_FIX': 'XYZ',
    #     'LIST': [1, 2, 3, 4, 5, ...]
    # }

CP2K Input Format
^^^^^^^^^^^^^^^^^

The generated CP2K input section looks like:

.. code-block:: text

    &MOTION
      &CONSTRAINT
        &FIXED_ATOMS
          COMPONENTS_TO_FIX XYZ
          LIST 1 2 3 4 5 6 7 8 21 22 23 24 26
        &END FIXED_ATOMS
      &END CONSTRAINT
      &MD
        ...
      &END MD
    &END MOTION

Notes
^^^^^

* The function modifies the input dictionary in-place and returns it
* Creates the ``MOTION/CONSTRAINT`` hierarchy if it doesn't exist
* Overwrites any existing ``FIXED_ATOMS`` section
* Components can be:
  
  * ``'XYZ'``: Fix all three Cartesian components (fully rigid)
  * ``'XY'``: Fix in-plane motion only
  * ``'Z'``: Fix perpendicular motion only


add_fixed_atoms_to_vasp_parameters()
-------------------------------------

.. autofunction:: teros.core.fixed_atoms.add_fixed_atoms_to_vasp_parameters

Description
^^^^^^^^^^^

Add selective dynamics to VASP POSCAR structure with ASE constraints.

This function creates a new StructureData with selective dynamics enabled. Fixed atoms are constrained using ASE's ``FixAtoms`` constraint, which translates to ``F F F`` flags in VASP POSCAR.

Parameters
^^^^^^^^^^

:param structure: Input slab structure
:type structure: orm.StructureData

:param fixed_atoms_list: List of 1-based atom indices to fix
:type fixed_atoms_list: list[int]

:param parameters: VASP INCAR parameters dictionary
:type parameters: dict

Returns
^^^^^^^

Tuple of (modified_parameters, constrained_structure):

* Modified INCAR parameters with ``IBRION``, ``POTIM``, ``NSW`` adjusted for constraints
* New StructureData with ASE FixAtoms constraint applied

:rtype: tuple[dict, orm.StructureData]

Example
^^^^^^^

.. code-block:: python

    from teros.core.fixed_atoms import (
        get_fixed_atoms_list,
        add_fixed_atoms_to_vasp_parameters,
    )

    # Get base parameters
    vasp_params = {
        'PREC': 'Normal',
        'ENCUT': 400,
        'IBRION': 0,
        'NSW': 100,
        # ... other INCAR tags ...
    }

    # Identify fixed atoms
    slab = orm.load_node(<STRUCTURE_PK>)
    fixed_list = get_fixed_atoms_list(slab, fix_type='bottom', fix_thickness=7.0)

    # Add constraints
    modified_params, constrained_slab = add_fixed_atoms_to_vasp_parameters(
        slab,
        fixed_list,
        vasp_params
    )

    # Use constrained_slab in VASP calculation
    # POSCAR will show selective dynamics:
    # Selective dynamics
    # Cartesian
    # 0.0 0.0 0.0  F F F  # Fixed
    # 1.5 0.0 0.0  T T T  # Free

VASP POSCAR Format
^^^^^^^^^^^^^^^^^^

The resulting POSCAR will have selective dynamics:

.. code-block:: text

    Material
     1.0
       10.0  0.0  0.0
        0.0 10.0  0.0
        0.0  0.0 25.0
    Ag O
      12  8
    Selective dynamics
    Cartesian
      0.000  0.000  0.000  F F F  # Fixed atom 1
      1.500  0.000  0.000  F F F  # Fixed atom 2
      ...
      0.000  0.000 20.000  T T T  # Free atom N-1
      1.500  0.000 20.000  T T T  # Free atom N

Notes
^^^^^

* Creates a new StructureData node; original structure is unchanged
* ASE's ``FixAtoms`` constraint is used internally
* Constraint information is stored in the ASE Atoms object
* VASP automatically reads selective dynamics flags from POSCAR
* Fixed atoms have zero forces in VASP output


Usage in Workflows
==================

Automatic Fixed Atoms (Recommended)
------------------------------------

Let PS-TEROS calculate fixed atoms automatically:

.. code-block:: python

    from teros.core.workgraph import build_core_workgraph

    wg = build_core_workgraph(
        workflow_preset='aimd_only',
        calculator='cp2k',  # or 'vasp'

        # ... structure and AIMD configuration ...

        # Fixed atoms (automatic calculation)
        fix_atoms=True,
        fix_type='bottom',
        fix_thickness=7.0,
        fix_elements=None,      # All elements
        fix_components='XYZ',   # Fully rigid

        name='AIMD_Fixed',
    )

PS-TEROS will:

1. Call ``get_fixed_atoms_list()`` for each slab
2. Call ``add_fixed_atoms_to_cp2k_parameters()`` or ``add_fixed_atoms_to_vasp_parameters()``
3. Generate proper input files with constraints

Manual Fixed Atoms (Advanced)
------------------------------

For custom control, pre-compute fixed atoms:

.. code-block:: python

    from teros.core.fixed_atoms import get_fixed_atoms_list
    from aiida import orm

    # Load slabs
    slabs = {
        'slab_111': orm.load_node(<STRUCTURE_PK_1>),
        'slab_110': orm.load_node(<STRUCTURE_PK_2>),
    }

    # Calculate fixed atoms for each
    fixed_atoms_lists = {}
    for label, slab in slabs.items():
        fixed = get_fixed_atoms_list(
            slab,
            fix_type='bottom',
            fix_thickness=7.0,
        )
        fixed_atoms_lists[label] = fixed
        print(f"{label}: {len(fixed)} atoms fixed")

    # Pass to AIMD function directly
    from teros.core.aimd_cp2k import aimd_single_stage_scatter_cp2k

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
        fixed_atoms_lists=fixed_atoms_lists,  # Pre-computed
        fix_components='XYZ',
    )

Troubleshooting
===============

Issue: No atoms fixed
---------------------

**Check:**

1. ``fix_type`` is not None
2. ``fix_thickness > 0``
3. Slab thickness exceeds ``fix_thickness``
4. Element symbols match if ``fix_elements`` is specified

.. code-block:: python

    # Debug: Print slab info
    from ase import Atoms
    atoms = slab.get_ase()
    positions = atoms.get_positions()
    z_coords = positions[:, 2]
    print(f"Slab Z range: {z_coords.min():.2f} to {z_coords.max():.2f} Å")
    print(f"Slab height: {z_coords.max() - z_coords.min():.2f} Å")
    print(f"Fix thickness: {fix_thickness:.2f} Å")

    # If fix_thickness > slab height, reduce it
    fix_thickness = (z_coords.max() - z_coords.min()) * 0.3  # Fix bottom 30%

Issue: Too many atoms fixed
---------------------------

**Solution:** Reduce ``fix_thickness`` or use element-specific fixing:

.. code-block:: python

    # Fix only oxygen in bottom 5 Å
    fixed = get_fixed_atoms_list(
        slab,
        fix_type='bottom',
        fix_thickness=5.0,
        fix_elements=['O']
    )

Issue: Constraints not applied
------------------------------

**For CP2K:**

Check the input file:

.. code-block:: bash

    verdi calcjob inputcat <CP2K_PK> aiida.inp | grep -A 5 CONSTRAINT

Should show ``&FIXED_ATOMS`` section.

**For VASP:**

Check POSCAR:

.. code-block:: bash

    verdi calcjob inputcat <VASP_PK> POSCAR | head -20

Should show ``Selective dynamics`` line and ``F/T`` flags.

See Also
========

* :doc:`/workflows/aimd-molecular-dynamics` - Full AIMD workflow with fixed atoms
* :doc:`/how-to/aimd-stages` - Multi-stage AIMD with constraints
* :doc:`aimd_cp2k` - CP2K AIMD implementation
* :doc:`aimd` - VASP AIMD implementation

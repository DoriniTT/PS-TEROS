.. _campaign-workgraph:

=========================================
Run references and slabs in one WorkGraph
=========================================

``build_vasp_campaign_workgraph`` puts the reference calculations and the slab
terminations of a phase diagram into one AiiDA WorkGraph, with one recipe for all
of them. Its outputs are nested by group, label and block, so that one graph PK
gives the analysis everything it needs. Use it when the references and the slabs
of a diagram should be submitted, and read back, as one unit.

The other builders keep their signatures, graphs and outputs. The
:doc:`reference builder <reference-thermochemistry>` runs the references alone, and
``build_surface_workgraph`` runs one calculation per structure. The signatures are
listed in the :doc:`API reference <api>`.

Build the graph
---------------

The example below is the SnO\ :sub:`2`\ (110) campaign of the VASP examples in the
repository: O\ :sub:`2`, rutile SnO\ :sub:`2` and α-Sn as references, each relaxed
and computed statically, with vibrations for the references, and three symmetric
terminations of SnO\ :sub:`2`\ (110). The runnable version with a command line is
``examples/vasp_campaign/campaign.py``.

.. code-block:: python

   import math

   import psteros

   TWO_PI = 2.0 * math.pi
   # Relaxed rutile SnO2 lattice in Å, from an earlier run. The slabs are cut from it;
   # see "Current limitation" below.
   RELAXED_A, RELAXED_C = 4.8301, 3.2434

   recipe = psteros.SurfaceWorkflowConfig(
       backend="vasp",
       calculation=psteros.VaspCalculationConfig(
           code_label="VASP-6.4@cluster",
           incar={"encut": 520, "prec": "Accurate", "ediff": 1e-7, "ismear": 0, "sigma": 0.05,
                  "lreal": False, "lasph": True, "ibrion": 2, "nsw": 100, "ediffg": -0.005},
           potential_mapping={"Sn": "Sn_d", "O": "O"},
           kpoints_spacing=0.03 * TWO_PI,  # Å^-1 with the 2*pi, as VASP's KSPACING
       ),
       execution=psteros.ExecutionPolicy(
           computer="cluster", queue="standard", max_concurrent_jobs=1,
           resources={"num_machines": 1, "num_mpiprocs_per_machine": 32},
           max_wallclock_seconds=12 * 3600,
       ),
       name="sno2_110",
   )

   cell_relax = psteros.CalculationOverride(parameters={"INCAR": {"isif": 3}})  # relax the cell too
   references = {
       "o2": psteros.ReferenceSystem(
           psteros.triplet_o2_cell(cell_length=12.0), "gas", symmetry_number=2, spin=1.0,
           override=psteros.CalculationOverride(
               parameters={"INCAR": {"ispin": 2, "nupdown": 2}},
               kpoints_distance=5.0 * TWO_PI,  # Gamma only in the 12 Å box
           ),
       ),
       "sno2": psteros.ReferenceSystem(
           psteros.rutile_sno2_bulk(), "solid", supercell=(2, 2, 3),
           block_overrides={
               "relax": cell_relax,
               "vibrations": psteros.CalculationOverride(
                   kpoints_distance=0.06 * TWO_PI, metadata={"max_wallclock_seconds": 24 * 3600}
               ),
           },
       ),
       "sn": psteros.ReferenceSystem(
           psteros.alpha_sn_bulk(), "solid", supercell=(2, 2, 2), block_overrides={"relax": cell_relax},
       ),
   }

   # Finished structures, cut from the relaxed lattice; ISIF = 2 keeps the cell fixed.
   slab_override = psteros.CalculationOverride(parameters={"INCAR": {"isif": 2, "ediffg": -0.02}})
   slabs = {
       f"slab_{termination}": psteros.SlabSystem(
           psteros.sno2_110_slab(
               termination=termination, triple_layers=3, vacuum_angstrom=15.0, a=RELAXED_A, c=RELAXED_C
           )[0],  # the structure; [1] is its SlabIdentity
           override=slab_override,
       )
       for termination in ("o", "sno", "sn2o")
   }

   graph = psteros.build_vasp_campaign_workgraph(
       references,
       slabs,
       recipe,
       reference_blocks=(psteros.Relax(), psteros.Static(), psteros.Vibrations(incar={"ediff": 1e-8})),
       slab_blocks=(psteros.Relax(), psteros.Static()),
       submit=False,  # build only; with submit=True, graph.pk is the PK of the campaign
   )

The labels of ``references`` and ``slabs`` name the tasks and the outputs; they must
follow the rules in :ref:`campaign-checks`. Each block takes the structure of the
block before it, or of the block named in ``structure_from``. The default blocks are
relax → static → vibrations for the references and relax → static for the slabs.
Pass ``reference_blocks`` or ``slab_blocks`` to change them, as this example does for
the vibrations.

A ``SlabSystem`` holds a finished structure, a pymatgen ``Structure`` or ``Slab``, an
AiiDA ``StructureData`` or a node PK, with the changes for its calculations. Its
``override`` applies to every block of the slab, and its ``block_overrides`` to the
block of the given name. A bare structure in ``slabs`` is taken as
``SlabSystem(structure)``. A slab is computed as a solid, in its own cell. Slabs are
finished structures in this version; see :ref:`campaign-limitation`.

The INCAR of each block is built as in the reference builder: the recipe INCAR, the
block defaults, the block's own ``incar``, the override of the reference or slab, its
override for that block, and finally the tags the block requires.

With ``submit=False`` the graph is built and checked, and nothing is submitted. The
runnable example prints every VASP task of the graph before it submits anything.

The outputs
-----------

The graph outputs are nested. Under ``references`` and ``slabs``, the label comes
first, then the block, then the port:

.. code-block:: text

   references
   ├── o2
   │   ├── relax          energy, structure, misc, remote, retrieved
   │   ├── static         energy, misc, remote, retrieved
   │   └── vibrations     frequencies, misc, remote, retrieved
   ├── sno2               relax, static and vibrations, as o2
   └── sn                 relax, static and vibrations, as o2
   slabs
   ├── slab_o
   │   ├── relax          energy, structure, misc, remote, retrieved
   │   └── static         energy, misc, remote, retrieved
   ├── slab_sno           relax and static, as slab_o
   └── slab_sn2o          relax and static, as slab_o

The ports are:

``energy``
   Electronic energy in eV of a relaxation or a static block, extrapolated to σ → 0.
``frequencies``
   Vibrational frequencies in cm\ :sup:`-1` of a vibrations block. Imaginary modes are
   negative.
``structure``
   The relaxed structure of a relaxation. A static block has no ``structure`` output;
   its input is the relaxed structure of the block before it.
``misc``, ``remote``, ``retrieved``
   The parsed VASP output (an AiiDA ``Dict``), the remote working folder and the
   retrieved files. Every block has them.

Read the results
----------------

Once the graph has run, its PK is all you need. The readers below take that PK. The
nested outputs are attributes of the graph node. As for the reference graphs,
aiida-workgraph attaches them only when every task has finished successfully.
``campaign_results`` and ``reference_results`` work while the graph runs and after a
failure. ``campaign_terminations`` and ``reference_thermochemistry`` need the blocks
they read to have finished.

.. code-block:: python

   from aiida import load_profile, orm

   load_profile()
   graph = orm.load_node(PK)
   graph.outputs.slabs.slab_o.static.energy.value                    # eV
   graph.outputs.references.sno2.vibrations.frequencies.get_list()   # cm^-1

``campaign_results(pk)``
   ``{"references": {label: {block: result}}, "slabs": {label: {block: result}}}``, with
   the same values as the outputs, as plain Python. Each result has ``kind``,
   ``state`` (for example ``finished`` or ``running``), ``pk`` of the VASP work chain,
   ``energy`` (eV) or ``frequencies`` (cm\ :sup:`-1`), ``structure`` (an AiiDA
   ``StructureData``: the relaxed structure of a relaxation, else the input of the
   block), ``misc`` (a dictionary), ``remote`` and ``retrieved``. A value is ``None``
   until its block has produced it.

   .. code-block:: python

      results = psteros.campaign_results(PK)
      results["references"]["sno2"]["static"]["energy"]    # eV, None until the static block has finished
      results["slabs"]["slab_sno"]["relax"]["state"]

``campaign_terminations(pk, *, energy_block=None, surfaces=2)``
   One ``SlabTermination`` per slab, ready for :doc:`the phase diagram <phase-diagram>`.
   The energy is that of the static block of each slab, or of the block named by
   ``energy_block``. The composition and the face area are read from the structure of
   that block. While the block has not finished, the function raises ``ValueError``,
   naming the slab and the state of the block.

   .. code-block:: python

      refs = psteros.campaign_results(PK)["references"]
      bulk = refs["sno2"]["static"]["structure"].get_pymatgen_structure()   # the cell of the energy
      metal = refs["sn"]["static"]["structure"].get_pymatgen_structure()
      oxide = psteros.BinaryOxideReferences(
          bulk_energy_ev=refs["sno2"]["static"]["energy"],
          bulk_composition=bulk.composition,
          oxygen_molecule_energy_ev=refs["o2"]["static"]["energy"],
          metal_energy_per_atom_ev=refs["sn"]["static"]["energy"] / len(metal),
      )
      diagram = psteros.surface_phase_diagram(psteros.campaign_terminations(PK), oxide)
      diagram.plot("sno2_110_phase_diagram.png")
      diagram.to_csv("sno2_110_phase_diagram.csv")

``reference_thermochemistry(pk)`` and ``reference_results(pk)``
   The same readers as for a reference graph, applied to the references of the
   campaign graph. ``reference_thermochemistry`` builds the objects of
   :doc:`reference-thermochemistry` from the electronic energy of the static block and
   the frequencies of the vibrations block. They give Δμ\ :sub:`O` as a temperature and
   an O\ :sub:`2` pressure.

   .. code-block:: python

      systems = psteros.reference_thermochemistry(PK)   # {"o2": IdealGasMolecule, "sno2": HarmonicSolid, ...}
      psteros.delta_mu_oxygen_ev(systems["o2"], 600.0, 0.2)   # Δμ_O (eV) of O2 at 600 K and 0.2 bar

.. _campaign-checks:

What is checked before anything is built
----------------------------------------

Each check runs before the graph creates any AiiDA node. Each error names the
offending label, block or field.

* ``config`` must be a VASP recipe, with ``backend="vasp"`` and a
  ``VaspCalculationConfig`` (``TypeError``).
* At least one reference or one slab is required (``ValueError``). The O\ :sub:`2`
  reference alone is a valid smoke test.
* ``config.role_overrides`` must be empty (``ValueError``). Put each change on the
  ``ReferenceSystem`` or ``SlabSystem`` it belongs to.
* A label starts with a letter and contains letters and digits, with single
  underscores between them. ``o2``, ``sno2`` and ``slab_sn2o`` are valid; ``2sno``,
  ``_sn``, ``sn_``, ``a__b``, ``sn-2`` and Python keywords such as ``class`` are not.
  The error lists the invalid labels (``ValueError``). Block names follow the same rule.
* A label is either a reference or a slab, never both (``ValueError``, naming the
  label).
* Each reference is a ``ReferenceSystem`` (``TypeError``).
* Overrides are ``CalculationOverride`` objects (``TypeError``). VASP overrides put
  their tags under ``"INCAR"`` only (``ValueError``), for references and slabs alike.
* ``block_overrides`` may name only blocks of their own group (``ValueError``, naming
  the label and the unknown blocks).
* A relaxation whose merged INCAR has ``NSW <= 0`` is an error, named by label and
  block.
* ``slab_blocks`` cannot contain a ``Vibrations`` block (``ValueError``).
* The block lists are checked as in ``psteros.blocks.check_blocks``: at least one
  block, unique names, and a ``structure_from`` that names an earlier block.

.. _campaign-limitation:

Current limitation: the slabs are given, not rebuilt
----------------------------------------------------

The slabs are finished structures, cut beforehand from a bulk lattice that was
relaxed beforehand, as the Quickstart of the repository README recommends. The graph
relaxes the references, but it does not rebuild the slabs from the bulk it relaxes.
Each slab keeps the lattice it was cut from, and the builder does not check that this
lattice is the one the bulk relaxes to. If the two differ, the slab and bulk energies
refer to different cells, and nothing in the graph warns about it.

Rebuilding the slabs from the relaxed bulk inside the graph is not part of this
version. To use the graph, relax the bulk first, in an earlier graph or run, and cut
the slabs from its relaxed lattice with ``sno2_110_slab(a=..., c=...)``, where ``a``
and ``c`` are the lattice constants in Å. The example does this with the lattice of an
earlier run.

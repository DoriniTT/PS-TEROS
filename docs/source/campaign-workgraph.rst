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

      oxide = psteros.campaign_references(PK, host="sno2")   # BinaryOxideReferences, from the same graph
      diagram = psteros.surface_phase_diagram(psteros.campaign_terminations(PK), oxide)
      diagram.plot("sno2_110_phase_diagram.png")
      diagram.to_csv("sno2_110_phase_diagram.csv")

   ``campaign_references`` builds the oxide reference of the same diagram for any binary or
   ternary oxide of the graph; see :ref:`campaign-any-material`.

``reference_thermochemistry(pk)`` and ``reference_results(pk)``
   The same readers as for a reference graph, applied to the references of the
   campaign graph. ``reference_thermochemistry`` builds the objects of
   :doc:`reference-thermochemistry` from the electronic energy of the static block and
   the frequencies of the vibrations block. They give Δμ\ :sub:`O` as a temperature and
   an O\ :sub:`2` pressure.

   .. code-block:: python

      systems = psteros.reference_thermochemistry(PK)   # {"o2": IdealGasMolecule, "sno2": HarmonicSolid, ...}
      psteros.delta_mu_oxygen_ev(systems["o2"], 600.0, 0.2)   # Δμ_O (eV) of O2 at 600 K and 0.2 bar

.. _campaign-any-material:

Any material: unary, binary, ternary
------------------------------------

The graph has the same layout for every material: references and slabs, each label with its
blocks and ports, as above. The layout does not count elements, so an elemental metal, a
binary oxide, a ternary oxide and an intermetallic are built and read in the same way. What
differs is the list of labels, and the readers turn the results into analysis input for each
material.

What the graph records
~~~~~~~~~~~~~~~~~~~~~~

A campaign graph submitted with ``submit=True`` stores a description in the node extra
``psteros_campaign`` (version 1). For every label of each group, it records the phase and the
composition of the input structure of that label: element counts of the cell whose energies
the blocks compute.

.. code-block:: text

   psteros_campaign                  node extra, version 1
     references.systems
       o2        {"phase": "gas",   "composition": {"O": 2}}
       sno2      {"phase": "solid", "composition": {"Sn": 2, "O": 4}}
     slabs.systems
       slab_o    {"phase": "solid", "composition": {"Sn": 6, "O": 12}}

The readers take the composition from the structure of the energy block when that block has
one, and from this record otherwise, so they work while the graph is still running.

``campaign_entries(pk)`` turns the graph into one ``CampaignEntry`` per label, references first
and in graph order:

.. code-block:: python

   for entry in psteros.campaign_entries(PK):
       print(entry.label, entry.phase, dict(entry.composition), entry.energy_ev, entry.energy_per_atom_ev)

``label``, ``group``
   The label, and its group: ``"references"`` or ``"slabs"``.
``phase``
   ``"gas"`` or ``"solid"``. A slab is always ``"solid"``.
``composition``
   Positive integer count of each element in the cell of ``energy_ev``, for example
   ``{"Sn": 6, "O": 12}``. The property ``atoms`` is their sum.
``energy_ev``
   Total electronic energy of that cell, in eV. ``None`` until the block has delivered it. The
   property ``energy_per_atom_ev`` (eV per atom) is ``None`` as long as ``energy_ev`` is.
``block``, ``state``
   The block that gives the energy, and its state. By default it is the last static block of the
   group, else its last relaxation; ``reference_energy_block`` and ``slab_energy_block`` name
   another, and are checked even when their group has no labels. The state is ``finished`` once
   the block has delivered its energy. Until then it says why: for example ``running`` or
   ``not started``. When the VASP work chain has finished but its energy task has not delivered,
   the state names that task, for example ``finished, slab_o_static_energy not started``.
``surface_area_angstrom2``
   For a slab, the area of one face in Å², the plane of its first two lattice vectors. ``None``
   for a reference.

The readers
~~~~~~~~~~~

The three general readers take a campaign PK, or a list of ``CampaignEntry`` objects (see
below), as their first argument. They share one rule for oxygen: when the campaign has an
O\ :sub:`2` gas reference, it is the oxygen reservoir, and an O atom or another single-element
O reference is not a candidate for it. Several candidates for one element are an error naming
them, unless ``reservoirs`` picks one, as in ``reservoirs={"O": "o2"}``.

``campaign_chemical_potentials(source, *, reservoirs=None, elements=None)``
   ``{element: eV per atom}`` of the elements that have a single-element reference, at the
   energy per atom of that reference: the elemental limits. For oxygen it is E(O\ :sub:`2`)/2.
   ``elements`` restricts the result, and its checks, to those elements, so that an unrelated
   reference does not matter; an element of ``elements`` without a reference is an error.
``campaign_references(source, *, host, reservoirs=None, exclude=(), independent=None)``
   The reference object of the phase diagram of the oxide ``host``: a ``BinaryOxideReferences``
   or a ``TernaryOxideReferences``, as the examples below show.
``campaign_surface_energies(source, *, chemical_potentials_ev=None, reservoirs=None, surfaces=2)``
   ``{slab label: γ in eV/Å²}``, with the formula below. ``chemical_potentials_ev`` (eV per atom)
   defaults to the elemental limits; without it, only the references of the elements of the slabs
   are read.

.. math::

   \gamma = \frac{E_\mathrm{slab} - \sum_i N_i\,\mu_i}{\mathrm{surfaces}\cdot A}

Here A is the area of one face, N\ :sub:`i` the count of element i in the slab, and μ\ :sub:`i` its
chemical potential in eV per atom. Multiply by ``psteros.EV_PER_ANGSTROM2_TO_J_PER_M2`` for J/m².

The default chemical potentials are for unary slabs only. For a slab of several elements,
``campaign_surface_energies`` asks for explicit ``chemical_potentials_ev`` and names the slabs. The
reason is that a compound is in equilibrium with its elements only where their chemical potentials
add up to the energy of its formula unit. The elemental limits, with each element at its own
reference, generally do not satisfy this, and surface energies taken there can come out negative.

Unary: Au
~~~~~~~~~

An elemental metal is a unary system. Its only reference is the metal, and its surface energies
need no phase diagram. The default chemical potential is exact here, because the slabs are unary.
This example cuts an Au(111) and an Au(100) slab from the bulk lattice, which should be the relaxed
one (see :ref:`campaign-limitation`):

.. code-block:: python

   import psteros
   from pymatgen.core import Lattice, Structure
   from pymatgen.core.surface import SlabGenerator

   bulk = Structure.from_spacegroup("Fm-3m", Lattice.cubic(A_AU), ["Au"], [[0, 0, 0]])  # A_AU: relaxed, in Å

   def cut(miller):
       generator = SlabGenerator(bulk, miller, min_slab_size=12.0, min_vacuum_size=15.0, center_slab=True)
       return generator.get_slabs(symmetrize=True)[0]

   graph = psteros.build_vasp_campaign_workgraph(
       {"au": psteros.ReferenceSystem(bulk, "solid")},
       {"au_111": cut((1, 1, 1)), "au_100": cut((1, 0, 0))},
       recipe,  # the VASP recipe of the example above
       reference_blocks=(psteros.Relax(), psteros.Static()),
       submit=True,
   )

   # once the graph has finished; graph.pk is its PK
   gamma = psteros.campaign_surface_energies(graph.pk)   # {"au_111": eV/Å², "au_100": eV/Å²}

Both faces of each slab count (``surfaces=2``). The chemical potential of gold is its bulk energy
per atom, so nothing else is needed. Multiply by ``psteros.EV_PER_ANGSTROM2_TO_J_PER_M2`` for J/m².

Binary: SnO2
~~~~~~~~~~~~

The SnO\ :sub:`2` campaign above is a binary oxide. ``campaign_references`` takes the host ``sno2``,
a solid of tin and oxygen, and finds the references that the model needs: the O\ :sub:`2` gas
reference for oxygen, and the single-element reference ``sn``, which sets the O-poor limit. It
returns the ``BinaryOxideReferences`` of :doc:`phase-diagram`:

.. code-block:: python

   oxide = psteros.campaign_references(PK, host="sno2")   # BinaryOxideReferences
   diagram = psteros.surface_phase_diagram(psteros.campaign_terminations(PK), oxide)
   diagram.plot("sno2_110_phase_diagram.png")

Without a single-element reference of the metal, the O-poor limit is unknown, and the diagram needs
``delta_mu_range``. The SnO\ :sub:`2` slabs have two elements, so their surface energies are not read
at the default chemical potentials: ``campaign_surface_energies`` asks for ``chemical_potentials_ev``,
and the phase diagram is the reader that gives them as a function of Δμ\ :sub:`O`.

Ternary: SrTiO3
~~~~~~~~~~~~~~~

The graph holds the host ``srtio3``, the single-element references ``sr`` and ``ti``, the O\ :sub:`2`
gas ``o2``, and two more oxides of the same elements, ``sro`` and ``tio2``. The competing phases of the
diagram are every other solid reference made of the elements of the host that is not a chosen
reservoir, so ``sro`` and ``tio2`` are added automatically. ``sr`` and ``ti`` are the chosen
reservoirs, the only candidates of their elements, so they are not competing phases. An unchosen
elemental phase is a competing phase too: if the graph also held a second Sr allotrope,
``reservoirs={"Sr": "sr"}`` would choose ``sr`` and make the other allotrope a competing phase.

.. code-block:: python

   oxide = psteros.campaign_references(PK, host="srtio3")   # TernaryOxideReferences
   diagram = psteros.ternary_surface_phase_diagram(psteros.campaign_terminations(PK), oxide)
   diagram.plot("srtio3_001_phase_diagram.png")

The element on the horizontal axis is Sr, the alphabetically first one by default;
``independent="Ti"`` puts Ti there instead. The choice changes the axes, not the physics.

Intermetallic: PdIn
~~~~~~~~~~~~~~~~~~~

PdIn has no phase-diagram model yet. ``campaign_references(PK, host="pdin")`` raises ``ValueError``,
so the intermetallic is read through the chemical potentials of its elements and through its surface
energies at chosen chemical potentials. The graph holds the single-element references ``pd`` and
``indium`` (``in`` is a Python keyword and cannot be a label), the compound ``pdin``, and slabs such
as ``pdin_110``.

Bulk PdIn is in equilibrium with its elements only where the chemical potentials of Pd and In add up
to the energy of one formula unit. The elemental values, with each element at its own reference, do
not satisfy this. The two ends of that line are the limits to use: the Pd-rich limit, with Pd at its
elemental value, and the In-rich limit, with In at its elemental value.

.. code-block:: python

   mu = psteros.campaign_chemical_potentials(PK, elements=("Pd", "In"))   # {"Pd": ..., "In": ...}, eV per atom
   pdin = next(entry for entry in psteros.campaign_entries(PK) if entry.label == "pdin")
   formula_unit_ev = 2 * pdin.energy_per_atom_ev                  # two atoms per PdIn formula unit

   pd_rich = {"Pd": mu["Pd"], "In": formula_unit_ev - mu["Pd"]}
   in_rich = {"Pd": formula_unit_ev - mu["In"], "In": mu["In"]}
   gamma_pd_rich = psteros.campaign_surface_energies(PK, chemical_potentials_ev=pd_rich)   # eV/Å²
   gamma_in_rich = psteros.campaign_surface_energies(PK, chemical_potentials_ev=in_rich)

Without ``chemical_potentials_ev``, ``campaign_surface_energies`` refuses the PdIn slab, because it
has two elements. A stoichiometric slab, with as many Pd as In atoms, gives the same γ at both
limits, because only the sum of the two chemical potentials enters; the limits differ for slabs that
are not stoichiometric.

When a reader refuses
~~~~~~~~~~~~~~~~~~~~~

The readers do not guess. Each refusal is a ``ValueError`` that names the labels, the block, the
element or the option involved:

* ``campaign_references`` needs a solid host with oxygen and one metal (a binary oxide) or two metals
  (a ternary oxide), and an O\ :sub:`2` gas reference. Any other host, such as an intermetallic or a
  nitride, has no phase-diagram model yet.
* Every other reference may contain only the elements of the host. A reference of another element,
  such as a solid of Si for SnO2, is an error, and so is a gas with more than one element, such as
  ``h2o``: neither is a reservoir of this model. ``exclude=("h2o",)`` leaves the named labels out.
* For a binary host, any other solid reference is an error: a solid oxygen, a compound such as SnO,
  or an elemental phase other than the chosen reservoir. For a ternary host the same solids are
  competing phases.
* Several candidates for one element are an error naming them, unless ``reservoirs`` chooses one. A
  label that is not a single-element reference of that element is an error too. For oxygen the
  candidates are the O\ :sub:`2` gas references, as described above.
* A host in ``exclude``, or a label of ``reservoirs`` in ``exclude``, is an error. So is an unknown
  label in ``exclude`` or ``reservoirs``, and ``independent`` for a binary host.
* An error from the reference model itself, such as an oxide that is unstable against its elements,
  starts with the host label, for example ``'sno2': ...``. The message for a missing O\ :sub:`2`
  reference lists the gas references of the campaign.
* ``campaign_surface_energies`` needs explicit ``chemical_potentials_ev`` for a slab of several
  elements, and ``reservoirs`` cannot be combined with ``chemical_potentials_ev``. A slab without a
  surface area, or with an element that has no chemical potential, is an error naming the slab (and
  the element).
* A reader needs the energies it uses. A block that has not delivered its energy is an error naming
  the label and the state of that block.
* Energies and areas must be finite numbers: NaN and infinity are rejected when a ``CampaignEntry``
  is made.

Hand-typed entries
~~~~~~~~~~~~~~~~~~

The readers also take a list of ``CampaignEntry`` objects in place of the PK, so energies from any
source work without AiiDA. Each label must appear once.

.. code-block:: python

   entries = [
       psteros.CampaignEntry(label="au", group="references", phase="solid", composition={"Au": 4}, energy_ev=e_au),
       psteros.CampaignEntry(label="au_111", group="slabs", phase="solid", composition={"Au": 48},
                             energy_ev=e_au_111, surface_area_angstrom2=area_111),
   ]
   gamma = psteros.campaign_surface_energies(entries)   # {"au_111": eV/Å²}

Here ``e_au`` is the energy in eV of the four-atom bulk cell, ``e_au_111`` the energy in eV of the
48-atom slab and ``area_111`` the area in Å² of one face. ``campaign_references`` and
``campaign_chemical_potentials`` take such a list in the same way.

.. _campaign-checks:

What is checked before anything is built
----------------------------------------

Each check runs while the graph is built, before anything is submitted. Each error names
the offending label, block or field.

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
* Two label/block pairs that give the same task names are an error naming both pairs. The
  reference ``sn`` with block ``o_relax`` and the slab ``sn_o`` with block ``relax`` both give
  ``sn_o_relax``, so one of them must be renamed (``ValueError``).
* A PK that does not belong to a structure, for example the PK of a ``Dict`` node, is a
  ``TypeError`` naming the label.

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

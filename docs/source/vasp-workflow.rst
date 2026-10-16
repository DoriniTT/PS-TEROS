.. _vasp-workflow:

==========================================
Prepare a VASP relaxation-to-static graph
==========================================

Use this guide after the :doc:`first tutorial <tutorial>` when you want a final
static VASP calculation on the geometry of a preceding relaxation, for the
bulk, the references and the slabs of a surface study. The guide builds the
graph only; it does not submit work.

Why use two calculation stages?
-------------------------------

A relaxation moves the atoms (and, with ``ISIF=3``, the cell) until the forces
fall below ``EDIFFG``. VASP keeps the plane-wave basis of the starting cell
during the relaxation, so the energy of a relaxation that changed the cell
carries a basis-set (Pulay) error. A static calculation on the relaxed
structure gives the energy to use in the thermodynamics. Keeping the stages
separate also makes the handoff visible in the AiiDA provenance record.

Before you start
----------------

You need an active AiiDA profile, a registered ``vasp.vasp`` code and an
uploaded POTCAR family (see :doc:`installation`). The labels and scheduler
settings below are placeholders. Replace them with values from your own
environment, and converge the numerical settings for your material before
submitting a calculation.

Choose one execution policy
---------------------------

Both stages must use the same execution policy, with the settings of your own
computer and scheduler. The builder turns its queue, resource, wall-time, and MPI choices into AiiDA task
metadata. The registered code in each recipe selects the actual AiiDA computer.

.. code-block:: python

   import psteros

   execution = psteros.ExecutionPolicy(
       computer="my-cluster",
       queue="my-queue",
       resources={"num_machines": 1, "num_mpiprocs_per_machine": 32},
       max_wallclock_seconds=86_400,
       max_concurrent_jobs=1,
   )

``max_concurrent_jobs`` limits how many calculations of the graph run at once;
1 runs them one after the other, ``None`` removes the limit. Options you leave
out, such as ``queue`` or ``max_wallclock_seconds``, are not sent to the
scheduler. The ``queue`` is passed to AiiDA as ``queue_name``, so the scheduler
plugin of your computer writes it in its own syntax (``#PBS -q``,
``#SBATCH --partition``).

Define the two recipes
----------------------

Both stages share the code, POTCAR family and mapping, cutoff and k-point
spacing; only the ionic settings differ. ``ENCUT`` and ``EDIFF`` are in eV,
``EDIFFG`` in eV/Å (negative: a force criterion) and ``kpoints_spacing`` in
Å⁻¹.

.. code-block:: python

   ELECTRONIC = {"ENCUT": 520, "PREC": "Accurate", "EDIFF": 1.0e-6,
                 "ISMEAR": 0, "SIGMA": 0.05, "LREAL": False}

   def vasp_recipe(name, ionic, overrides=None):
       return psteros.SurfaceWorkflowConfig(
           backend="vasp",
           calculation=psteros.VaspCalculationConfig(
               code_label="vasp@my-cluster",
               incar={**ELECTRONIC, **ionic},
               potential_family="PBE",
               potential_mapping={"Sn": "Sn_d"},
               kpoints_spacing=0.20,
               max_iterations=3,      # restarts, e.g. after the walltime
           ),
           execution=execution,
           name=name,
           role_overrides=overrides or {},
       )

   RELAX = {"IBRION": 2, "NSW": 200, "ISIF": 2, "EDIFFG": -0.01}
   STATIC = {"IBRION": -1, "NSW": 0}

The builder refuses a relaxation INCAR that does not move the ions
(``IBRION >= 0`` and ``NSW > 0``) and a static INCAR that does.

.. _vasp-precision:

Or choose a precision level
---------------------------

Instead of writing the cutoff, the k-point spacing and the convergence tags
yourself, name a level: ``"fast"``, ``"moderate"`` or ``"precise"``. The
level fills in what your recipe leaves out, and everything you write
explicitly wins.

.. code-block:: python

   RELAX = {"IBRION": 2, "NSW": 200, "ISIF": 2}     # only the ionic settings
   STATIC = {"IBRION": -1, "NSW": 0}

   relax = psteros.VaspCalculationConfig("vasp@my-cluster", RELAX, precision="moderate")
   static = psteros.VaspCalculationConfig("vasp@my-cluster", STATIC, precision="moderate")
   # or, for a calculation with nothing else to say: VaspCalculationConfig("vasp@my-cluster")

============  =========  ===================  =========  =====================
level         ``ENCUT``  ``kpoints_spacing``  ``EDIFF``  relaxation ``EDIFFG``
============  =========  ===================  =========  =====================
``fast``      400 eV     0.44 Å⁻¹             1e-5 eV    -0.05 eV/Å
``moderate``  550 eV     0.31 Å⁻¹             1e-6 eV    -0.03 eV/Å
``precise``   700 eV     0.19 Å⁻¹             1e-7 eV    -0.01 eV/Å
============  =========  ===================  =========  =====================

Every level also sets ``PREC = Accurate``, Gaussian smearing (``ISMEAR = 0``,
``SIGMA = 0.05``), ``LREAL = False``, ``ALGO = Normal``, ``LASPH = True``,
``NELM = 200``, ``LWAVE`` and ``LCHARG`` off and PBE (``GGA = PE``). The rules:

* **Your values win.** A tag in ``incar``, in a ``CalculationOverride`` or in
  a block overrides the level, in any upper or lower case, and an explicit
  ``kpoints_spacing`` (even 0.20) overrides the level's spacing.
* **``EDIFFG`` only for a relaxation**, that is when the final INCAR of the
  structure has ``IBRION`` 1, 2 or 3 and ``NSW > 0``; a static or a
  finite-difference calculation gets none. A static INCAR that an override
  turns into a relaxation gets it too.
* **Spin is off unless the structure is magnetic.** ``ISPIN = 1``, except when
  an element has a magnetic moment in the aiida-vasp protocol (Fe, Co, Ni, Mn,
  Cr, V, ...): then ``ISPIN = 2`` and ``MAGMOM`` with those moments. Give
  ``ISPIN`` or ``MAGMOM`` yourself for anything else, such as a triplet molecule.
* **POTCARs.** An element that ``potential_mapping`` does not name gets the
  POTCAR of the aiida-vasp protocol (``Sn_d``, ``Sr_sv``, ``Ti_pv``, ``O``, ...);
  the POTCAR *family* is still yours (``potential_family``).
* **A recipe without a level is unchanged.** ``precision`` is ``None`` unless
  you set it, or leave ``incar`` out altogether (then it is ``"moderate"``).
  The level applies to the surface builders and to the reference blocks.

aiida-vasp has one VASP protocol, ``balanced``, and no ``fast`` or
``stringent`` one, so PS-TEROS defines the three levels. ``moderate`` is the
``balanced`` protocol (cutoff, k-point spacing, smearing, electronic settings,
POTCARs, read from the installed aiida-vasp and checked by a test); ``fast`` and
``precise`` follow the ends of the scan of aiida-vasp's convergence work chain.
Three things differ from ``balanced`` on purpose: the functional is PBE
(``balanced`` uses PBEsol, ``GGA = PS``), the calculation is not spin-polarised
unless magnetic, and ``EDIFF`` is one number per cell instead of per atom, so
that the energies of cells of different sizes are comparable. The levels are
starting points: converge the property you compute before you trust it, and
tighten ``EDIFF`` (``1e-7`` or less) for the displaced calculations of
:doc:`vibrations`.

Per-structure settings
----------------------

``CalculationOverride`` changes one calculation of a recipe, keyed by its
label. Typical overrides of a surface study:

.. code-block:: python

   slab, _ = psteros.sno2_110_slab(termination="o", triple_layers=9, vacuum_angstrom=20.0)

   cell_relaxation = psteros.CalculationOverride(parameters={"INCAR": {"ISIF": 3}})
   triplet_o2 = psteros.CalculationOverride(
       parameters={"INCAR": {"ISPIN": 2, "MAGMOM": [1.0, 1.0]}},
       kpoints_distance=10.0,          # Gamma only for the molecule
   )
   frozen_centre = psteros.CalculationOverride(fixed_sites=psteros.central_sites(slab, half_width=1.5))

   relax = vasp_recipe("sno2_110", RELAX, {
       "bulk": cell_relaxation, "metal": cell_relaxation,
       "o2": triplet_o2, "sno2_110_o": frozen_centre,
   })
   static = vasp_recipe("sno2_110_static", STATIC, {"o2": triplet_o2})

* ``parameters={"INCAR": {...}}`` updates the INCAR of that structure.
* ``kpoints_distance`` replaces the k-point spacing; a large value gives a
  Γ-only mesh for a molecule in a box.
* ``fixed_sites`` keeps the listed sites fixed with selective dynamics;
  ``psteros.central_sites(slab, half_width)`` selects the sites near the
  mid-plane of a slab.
* A slab with two different faces needs a dipole correction:
  ``{"LDIPOL": True, "IDIPOL": 3, "DIPOL": [0.5, 0.5, 0.5]}``.

``ChargeNeutralSurfaceStudy.vasp_overrides(stage)`` and
``PolarSurfaceStudy.vasp_overrides(stage)`` write all of these for you, for
``stage="relax"`` and ``stage="static"``.

Connect the stages
------------------

Pass the labelled starting structures and the two recipes to the builder. Each
static task receives the relaxed ``StructureData`` of its relaxation.

.. code-block:: python

   structures = {
       "bulk": psteros.rutile_sno2_bulk(),
       "metal": psteros.alpha_sn_bulk(),
       "o2": psteros.triplet_o2_cell(cell_length=12.0),
       "sno2_110_o": slab,
   }
   graph = psteros.build_relax_static_workgraph(structures, relax, static)

   print(graph.name)
   print(graph.max_number_jobs)

You should see:

.. code-block:: text

   sno2_110_relax_static
   1

The graph holds two tasks per label, ``<label>_relax_vasp`` and
``<label>_static_vasp``, and exposes ``<label>_relaxed_structure``,
``<label>_relax_misc`` and ``<label>_static_misc`` (the energies) for each.
This remains an unsubmitted graph because ``submit`` was not set to ``True``.

Read the results
----------------

When the graph has finished, ``read_vasp_results`` returns the static energies
(``energy_extrapolated``, eV) and the relaxed structures by label:

.. code-block:: python

   from aiida import orm

   energies, relaxed = psteros.read_vasp_results(orm.load_node(graph.pk), structures)

They go straight into :doc:`phase-diagram`, or into the ``analyse`` method of
``ChargeNeutralSurfaceStudy`` and ``PolarSurfaceStudy``.

Before submission
-----------------

Check the items that the library cannot choose for you:

* the starting structure and any atoms you intend to constrain;
* convergence of ``ENCUT``, k-point sampling, vacuum, and slab thickness, and
  an ``ENCUT`` about 30 % above the largest ``ENMAX`` for cell relaxations;
* that every calculation you plan to compare uses the same POTCARs, ``ENCUT``,
  ``PREC`` and functional; and
* the queue, resources, wall-time, and MPI settings accepted by your AiiDA
  computer and scheduler.

The :doc:`SnO2 surface-energy model <examples>` explains how compatible bulk,
slab, and oxygen-reference energies are combined after the calculations finish,
and :doc:`phase-diagram` turns them into a surface phase diagram. The same graph
can be built with Quantum ESPRESSO, see :doc:`qe-workflow`.

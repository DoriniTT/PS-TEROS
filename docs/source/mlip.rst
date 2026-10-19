.. _mlip:

=====================================================
Relax terminations with a fine-tuned MACE model
=====================================================

Relaxing every termination with DFT limits how many terminations a study can
afford. The MLIP route of psteros replaces those relaxations: a few DFT
*static* calculations teach a MACE foundation model your system, the model
relaxes every slab on a GPU, and one more DFT static at each relaxed geometry
gives the energies of the surface phase diagram. In one line: **MACE finds the
relaxed geometry, DFT gives its energy.**

It is the route psteros recommends when you have a GPU: the slab relaxations,
usually the bottleneck of a study, cost minutes on the GPU instead of DFT
hours, so a study can afford many more terminations. Without a GPU, relax the
slabs with DFT as before (:doc:`vasp-workflow`, :doc:`qe-workflow`,
:doc:`cp2k-workflow`); nothing in those workflows changes.

.. figure:: ../images/mlip_route.png
   :width: 100%
   :alt: Flowchart of the MLIP route

   The MLIP route: DFT statics, the electron-count check, a replay fine-tune
   and a MACE relaxation on a GPU, one active-learning round, the final DFT
   statics that enter the diagram and, optionally, the model's vibrations.

.. warning::

   Converge the DFT recipe for your system and check a new system against a
   few DFT relaxations before relying on the route. Use the final DFT
   energies (``"mlip_dft"``), not the model's own (``"mlip"``).

Metals, semiconductors and insulators
-------------------------------------

The route itself does not care whether the compound is a metal. The slabs do:
in a metal any cut of the bulk has a sensible electron count, but a slab of a
semiconductor or an insulator must be charge-neutral, or passivated, which
pymatgen's ``SlabGenerator`` alone does not ensure. So psteros asks for one of
two things:

* **Slabs cut by psteros** with electron counting:
  ``MlipRoute(ChargeNeutralSurfaceStudy(...))`` (symmetric, charge-neutral
  slabs of an element, a binary or a ternary compound) or
  ``MlipRoute(PolarSurfaceStudy(...))`` (polar faces on a pseudo-hydrogen
  passivated bottom); see :ref:`mlip-route-studies`. These run for any
  compound, with no metal check.
* **Your own slabs** (``MlipSurfaceStudy``, ``build_mlip_surface_workgraph``):
  the graph checks that the bulk is a metal (partly filled states at the
  Fermi level, ``require_metallic``) and stops otherwise. With
  ``allow_non_metallic=True`` it runs anyway and warns you that the electron
  count of the slabs is then your responsibility.

The recipe
----------

1. **Relax the references with DFT, once**: the bulk, the metal and the gas
   molecule in a box, with the recipe of the study. Cut the slabs at the
   relaxed bulk cell and do not relax them.
2. **One DFT static of every structure** (bulk, references, slabs as built):
   energies and forces, the training labels.
3. **Check the electron count**: slabs cut by psteros pass; your own slabs
   need a metallic bulk (or ``allow_non_metallic=True``).
4. **Replay fine-tune of MACE-MH-1** on a GPU: a new head learns the DFT
   labels while the original head keeps training on the foundation model's
   own labels of ~150 small structures of the same elements, so the model does
   not forget general chemistry.
5. **MACE relaxes every structure**, cell fixed, slab centre held in place.
6. **One active-learning round**: one DFT static at each relaxed geometry
   joins the labels, the model is fine-tuned again from MH-1 on all of them
   and relaxes the slabs again from the start. The first model has only seen
   unrelaxed slabs; these statics add labels near the relaxed geometries.
7. **One final DFT static at each relaxed geometry**; these energies, for the
   slabs and the references alike, enter the diagram.

That is three DFT statics per termination instead of a DFT relaxation.
Optionally, the final model also gives the harmonic vibrations of every
structure (see :ref:`mlip-vibrations`).

Run it: ``MlipSurfaceStudy``
----------------------------

.. code-block:: python

   import psteros
   from aiida import orm

   study = psteros.MlipSurfaceStudy(bulk, {"Ru": ru_metal, "O": o2_box}, slabs)
   print(study)      # labels, fixed sites, number of DFT statics and GPU jobs

   static = psteros.SurfaceWorkflowConfig(backend="cp2k", name="ruo2", calculation=cp2k_recipe,
                                          execution=cpu_execution, role_overrides={"ref_O": triplet_o2})
   fine_tune = psteros.MaceFineTuneConfig("/work/me/models/mace-mh-1.model",
                                          replay_file="replay-data-mh-1-omat-pbe.xyz")
   job = psteros.MlipJobConfig("mace-python@gpu-cluster", psteros.ExecutionPolicy(
       queue="gpu", with_mpi=False, max_wallclock_seconds=6 * 3600,
       resources={"num_machines": 1, "num_mpiprocs_per_machine": 8, "ngpus_per_machine": 1}))
   graph = study.build(static, fine_tune, job=job, submit=True)

   result = study.analyse(orm.load_node(graph.pk))   # the final DFT energies only
   print(result)
   result.diagram.plot("ruo2_110_mlip.png")

``MlipSurfaceStudy(bulk, references, slabs, fixed_half_width_angstrom=1.5,
rattled=0)`` labels the calculations ``bulk``, ``ref_<element>`` and the slab
labels; it keeps the sites within ``fixed_half_width_angstrom`` of each slab
mid-plane in place during the MACE relaxations and adds ``rattled`` rattled
copies of each slab to the training statics. ``analyse`` also takes a mapping
of energies typed by hand. A full script is
``examples/mlip_ruo2/ruo2_mlip_graph.py``; the guide is
``docs/MLIP_SURFACES.md``.

Instead of the active-learning round, ``rattled=5`` adds five rattled copies
of each slab to the first training set; with ``active_learning_rounds=0`` one
fine-tune then runs, with more DFT statics up front.

.. _mlip-route-studies:

Run it: ``MlipRoute`` with the studies of psteros
-------------------------------------------------

``MlipRoute(study)`` runs a ``ChargeNeutralSurfaceStudy`` or a
``PolarSurfaceStudy`` on the MLIP route. The study cuts the slabs (one or
more orientations, any number of elements) and keeps its labels, its settings
and its analysis; the route replaces the DFT relaxations:

.. code-block:: python

   study = psteros.ChargeNeutralSurfaceStudy(
       sno2_bulk, [(1, 1, 0), (1, 0, 0)], references={"Sn": tin, "O": psteros.triplet_o2_cell()})
   route = psteros.MlipRoute(study)        # the centre of each slab held in place, as above
   print(route)                            # slabs, fixed sites, DFT statics and GPU jobs

   static = psteros.SurfaceWorkflowConfig(
       backend="vasp", name="sno2", execution=cpu_execution,
       calculation=psteros.VaspCalculationConfig(
           "vasp@cpu-cluster", incar={"ENCUT": 520, "EDIFF": 1e-6, "ISMEAR": 0, "SIGMA": 0.05,
                                      "IBRION": -1, "NSW": 0, "LREAL": False}),
       role_overrides=study.vasp_overrides("static"))
   graph = route.build(static, fine_tune, job=job, submit=True)

   result = route.analyse(orm.load_node(graph.pk))   # the study's own analysis, on the final DFT energies
   result.diagram.plot("sno2_mlip.png")

``route.read(graph)`` gives the energies and the MACE-relaxed structures, the
same pair as ``read_vasp_results``. A ternary study gives a
``TernarySurfacePhaseDiagram``, an element its surface energies.

**Polar surfaces.** ``MlipRoute(PolarSurfaceStudy(...))`` runs the polar
faces with VASP (the pseudo-hydrogen POTCARs); use
``study.potential_mapping(...)`` and ``study.vasp_overrides("static")`` in the
static recipe. Every atom of a polar slab relaxes, as in the polar method, and
``analyse`` returns the ``PolarStudyResult`` with the bottom check of the
MACE-relaxed slabs. For MACE the pseudo-hydrogens are hydrogen atoms; the
fine-tune learns them from their DFT labels like any other atom.

.. _mlip-vibrations:

Vibrations on the GPU
---------------------

With ``vibrations=psteros.VibrationsConfig(...)``, ``build`` adds one GPU job
in which the deployed model of the last round gives the harmonic modes of
every structure at its relaxed geometry, by the same Γ-point finite
differences as ``build_vibrations_workgraph`` (:doc:`vibrations`): the
sites held in place in the relaxation are not displaced (their region is
counted as bulk, so it must have the bulk stoichiometry, which ``build``
checks), a bulk is displaced in its ``supercells`` entry, and the gas
references are molecules. No DFT calculation is added, so a large supercell
costs little. ``analyse(graph, temperature_k=...)`` then uses free energies:

.. code-block:: python

   graph = route.build(static, fine_tune, job=job, submit=True,
                       vibrations=psteros.VibrationsConfig(supercells={"bulk": (3, 3, 3), "ref_Sn": (3, 3, 3)}))
   result = route.analyse(orm.load_node(graph.pk), temperature_k=800.0)

The energies are still the final DFT statics; only the vibrational terms come
from the model. These are not mixed sources in the sense below: a constant
offset per element has no curvature, so it does not enter the modes. The
modes of the model have not yet been compared with DFT finite differences for
a series of systems, so check a few: ``route.dft_vibrations(graph, static)``
builds the finite-difference graph of ``build_vibrations_workgraph`` on the
MACE geometries with the same fixed sites (QE, VASP, and here also CP2K with
``RUN_TYPE ENERGY_FORCE``), and
``route.analyse(graph, temperature_k=T, vibrations=psteros.read_vibrations(check, labels))``
uses those modes instead. ``read_vibrations`` reads both graphs.
Polar slabs have no vibrations yet (their pseudo chemical potentials are 0 K
energies).

One energy source per diagram
-----------------------------

With

.. math::

   \gamma = \frac{1}{nA}\Big[E_\mathrm{slab} - \sum_i N_i\,\mu_i\Big],
   \qquad \mu_i = E_i^\mathrm{ref} + \Delta\mu_i ,

adding a constant :math:`\delta_Z` per atom of element :math:`Z` to every energy
changes :math:`\gamma` by :math:`\sum_Z N_Z \delta_Z - \sum_Z N_Z \delta_Z = 0`, but
only if the slab and every reference carry the same :math:`\delta_Z`: they must
come from the same model, or the same DFT. A model's slab energy with a DFT
bulk energy brings the offsets straight into :math:`\gamma`.
``check_one_energy_source`` raises a ``ValueError`` that names the labels of
each source (``"dft"``, ``"mlip"``, ``"mlip_dft"``). The same identity makes the
E0s of the fine-tuned head irrelevant for :math:`\gamma`; ``fit_element_shifts``
only makes the fit start close.

The building blocks (job folders, no AiiDA)
--------------------------------------------

1. **Labels.** DFT static calculations (for example with
   :doc:`cp2k-workflow`) of the bulk, the elemental references, the gas
   molecule and every termination as built, each a ``LabelledStructure``
   (energy in eV, forces in eV/Å). With smearing, use the free energy, whose
   derivatives are the forces.
2. **Fine-tune.** ``MaceFineTuneConfig`` holds the recipe (default: MACE-MH-1,
   head ``omat_pbe``, float64, batch 1, constant learning rate 1e-3, EMA
   0.999, universal loss with energy and force weights 10, 100 passes, the
   final EMA model). With ``replay=True`` it is a multi-head replay fine-tune:
   ``select_replay_structures`` picks from the foundation model's replay data
   every distinct structure made only of the training elements;
   ``prepare_finetune_job`` writes a job folder (``job.json``, the labels,
   the replay geometries and the script ``finetune.py``) to run on the GPU
   computer, one model per seed. The job labels the replay structures with
   the foundation model, fits the E0 shifts, checks that the untrained model
   reproduces the foundation model, trains and evaluates both heads.
3. **Relax.** ``prepare_relax_job`` writes a folder in which every model
   relaxes every structure with ASE (``relax.py``; positions only, fixed sites
   kept, ``MlipRelaxConfig``: FIRE to 0.01 eV/Å by default).
   ``relax_with_calculator`` does the same in Python with any ASE calculator.
4. **Diagram.** Type the energies into ``SlabTermination`` and the reference
   classes of :doc:`phase-diagram`, all from one source.

.. code-block:: python

   import psteros

   config = psteros.MaceFineTuneConfig(
       "/work/me/models/mace-mh-1.model",
       foundation_model_sha256="a522eb7f59c7879963d41586528f4980baf33e086c94aa92e3eafdeccad3be47",
       replay_file="replay-data-mh-1-omat-pbe.xyz",
   )
   train, valid, record = psteros.select_replay_structures(config.replay_file, elements=["Ru", "O"])
   folder = psteros.prepare_finetune_job("jobs/ruo2_m2", frames, config, name="ruo2_m2",
                                         replay_train=train, replay_valid=valid)
   # run "python finetune.py" in that folder on a GPU node, then:
   models = {f"m2_{seed}": (f"jobs/ruo2_m2/models/ruo2_m2_seed{seed}.model", "Default") for seed in config.seeds}
   psteros.prepare_relax_job("jobs/ruo2_relax", slabs, models,
                             fixed_sites={label: psteros.central_sites(slab) for label, slab in slabs.items()})

``rutile_110_slab`` builds the six RuO\ :sub:`2`\ (110) terminations of the
example (``RUTILE_110_TERMINATIONS``) for any rutile MO\ :sub:`2`.


The whole route as one WorkGraph
--------------------------------

``build_mlip_surface_workgraph`` wires the route for QE, VASP or
CP2K statics and a GPU computer with a MACE environment:

1. one DFT static of every structure (slabs and references, the references
   already relaxed with the same recipe) and of the optional
   ``training_structures`` (for example ``rattled_copies`` of the slabs);
2. a replay fine-tune (``MaceFineTuneConfig``) on those labels, after checking
   that ``metallic_bulk`` is a metal;
3. the model of the first seed relaxes every structure (cell fixed,
   ``fixed_sites`` kept in place); the other seeds relax them too, as a
   committee in ``mlip_round<r>_relax_results``;
4. ``active_learning_rounds`` times (1 by default): one DFT
   static at each relaxed geometry joins the labels, a new fine-tune from the
   foundation model trains on all of them, and the new model relaxes the
   structures again from the start;
5. with ``dft_check`` (default), one DFT static at each final geometry;
6. with ``vibrations``, the harmonic modes of every structure from the
   deployed model (``<label>_mlip_vibrations``).

``allow_non_metallic=True`` skips the metal check of step 2, with a warning
about the electron count of the slabs. Every task is wired when the graph is built; calcfunctions only turn finished
calculations into the next job's inputs (the training folder, the relaxed
structure). The fine-tune and relaxation jobs are aiida-shell jobs that run
``finetune.py`` and ``relax.py`` with the Python of ``MlipJobConfig.code_label``;
a relaxation links the remote folder of its fine-tune, so both run on that
computer and the folder, which holds the models, is never cleaned. A
relaxation that does not reach ``fmax`` stops the graph instead of sending an
unrelaxed geometry to DFT.

.. code-block:: python

   import psteros

   static = psteros.SurfaceWorkflowConfig(backend="cp2k", name="ruo2", calculation=cp2k_recipe,
                                          execution=cpu_execution, role_overrides={"o2": triplet})
   fine_tune = psteros.MaceFineTuneConfig("/work/me/models/mace-mh-1.model",
                                          replay_file="replay-data-mh-1-omat-pbe.xyz")
   job = psteros.MlipJobConfig("mace-python@gpu-cluster", psteros.ExecutionPolicy(
       queue="gpu", with_mpi=False, max_wallclock_seconds=6 * 3600,
       resources={"num_machines": 1, "num_mpiprocs_per_machine": 8, "ngpus_per_machine": 1}))
   graph = psteros.build_mlip_surface_workgraph(
       structures, static, fine_tune, job=job, metallic_bulk="bulk",
       fixed_sites={label: psteros.central_sites(structures[label]) for label in slabs}, submit=True)
   # when it has finished (or while it runs, for the parts that have):
   energies = psteros.read_mlip_energies(graph_node, list(structures), energy_source="mlip_dft")
   relaxations = psteros.read_mlip_results(graph_node, list(structures))

The energies of ``read_mlip_energies`` all come from one source, so the slabs
and the references go into one diagram. The route needs three DFT statics per
termination (start, active learning, check), however many steps a DFT
relaxation of that termination would take.

.. _mlip:

=====================================================
Relax terminations with a fine-tuned MACE model
=====================================================

Relaxing every termination with DFT limits how many terminations a study can
afford. The opt-in MLIP route of psteros replaces those relaxations: DFT
*static* calculations label a fine-tune of a MACE foundation model, the
fine-tuned model relaxes the terminations on a GPU, and the surface phase
diagram uses its energies, or one DFT static on each relaxed geometry. Nothing
in the existing workflows changes when you do not use it.

.. warning::

   The route is for **metallic systems only** (metals, alloys, metallic
   oxides such as RuO\ :sub:`2`). ``require_metallic`` stops it when the bulk
   has no partly filled states at the Fermi level. Semiconductors and
   insulators need a separate treatment. The route is being validated on
   RuO\ :sub:`2`\ (110) against full DFT relaxations and Reuter and Scheffler,
   Phys. Rev. B 65, 035406 (2002); see ``docs/plans/MACE_SURFACE_PLAN.md``.

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

The steps
---------

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
validation (``RUTILE_110_TERMINATIONS``) for any rutile MO\ :sub:`2`.

.. _qe-first-workflow:

=================
QE-first workflow
=================

Backend policy
--------------

psteros uses Quantum ESPRESSO through ``aiida-quantumespresso`` as its
primary calculation backend.  A QE recipe consists of a code label,
an ``aiida-pseudo`` family, explicit QE namelists, reciprocal-space spacing,
and an execution policy.  VASP remains implemented through
``VaspCalculationConfig`` for an existing VASP project, but is not the
default or the reference path for new psteros work.

The Bohr reference policy is deliberately narrow: one machine, one MPI rank,
the ``gpu_a100`` queue, and ``max_concurrent_jobs=1``.  It is enforced by
``ExecutionPolicy`` and by the containing project controller.  A psteros
graph is provenance-preserving, but it does not invent a second monitoring
service; an Odyssey Auto project owns polling, recovery, and archival.

Relaxation followed by a static calculation
-------------------------------------------

Use the two-stage builder when the final energy must use a relaxed structure.
Both stages are serial, including when the mapping contains several slabs.

.. code-block:: python

   import psteros

   policy = psteros.ExecutionPolicy(max_wallclock_seconds=86400)
   relax = psteros.SurfaceWorkflowConfig(
       backend="qe",
       calculation=psteros.QeCalculationConfig(
           code_label="QE-7.6-PW-GPU-A100@bohr",
           pseudo_family="SSSP/1.3/PBE/precision",
           parameters={
               "CONTROL": {
                   "calculation": "relax", "tstress": True, "tprnfor": True,
                   "forc_conv_thr": 1.0e-3, "nstep": 200,
               },
               "SYSTEM": {"ecutwfc": 75.0, "ecutrho": 600.0},
               "ELECTRONS": {"conv_thr": 1.0e-10},
               "IONS": {"ion_dynamics": "bfgs"},
           },
       ),
       execution=policy,
       name="sno2_110_relax",
   )
   static = psteros.SurfaceWorkflowConfig(
       backend="qe",
       calculation=psteros.QeCalculationConfig(
           code_label="QE-7.6-PW-GPU-A100@bohr",
           pseudo_family="SSSP/1.3/PBE/precision",
           parameters={
               "CONTROL": {"calculation": "scf", "tstress": True, "tprnfor": True},
               "SYSTEM": {"ecutwfc": 75.0, "ecutrho": 600.0},
               "ELECTRONS": {"conv_thr": 1.0e-10},
           },
       ),
       execution=policy,
       name="sno2_110_static",
   )
   o_slab, _ = psteros.sno2_110_slab(termination="o", triple_layers=9)
   graph = psteros.build_qe_relax_static_workgraph(
       {"sno2_110_o": o_slab}, relax, static
   )

The final SCF task receives the relaxed ``StructureData`` output directly,
while pseudopotentials are mapped from the immutable input structure.  This
avoids rebuilding an intermediate structure in client code.

For a constrained QE relaxation, use
``psteros.qe_fixed_coordinate_flags`` rather than constructing raw booleans.
In aiida-quantumespresso, ``True`` means *fixed* and is written as QE's
``0`` coordinate flag; ``False`` means free and is written as ``1``.  For
example, to fix only a selected central set of sites:

.. code-block:: python

   override = psteros.CalculationOverride(
       settings={"FIXED_COORDS": psteros.qe_fixed_coordinate_flags(
           len(o_slab), fixed_site_indices={12, 13, 14, 15, 16, 17},
       )}
   )

SnO2(110) analysis convention
------------------------------

For a symmetric rutile SnO2(110) slab in equilibrium with bulk SnO2,
``psteros.surface_energy_oxide_equilibrium`` evaluates

.. math::

   \gamma(\Delta\mu_O) =
   \frac{E_{slab} - N_{Sn}E_{SnO_2} - (N_O - 2N_{Sn})
   [E_{O_2}/2 + \Delta\mu_O]}{2A}.

The chemically allowed oxygen interval must be established from the same
functional, pseudopotentials, and converged numerical settings.  In the
maintained reference campaign it is constrained by alpha-Sn and litharge SnO
in addition to the O-rich limit ``Delta mu_O <= 0``.  Do not combine energy
references from another backend or pseudopotential family.

Reproduced SnO2(110) case study
--------------------------------

The maintained end-to-end psteros example is a real calculation of the clean,
symmetric rutile SnO2(110) O, SnO, and Sn2O terminations.  The frozen models
have nine triple layers, 20 Angstrom of vacuum, and formulas Sn18O36,
Sn18O34, and Sn18O32.  psteros created the AiiDA QE relax/static workflows;
each calculation used QE 7.6 with ``SSSP/1.3/PBE/precision``, PBE, 90/720 Ry
wavefunction/charge-density cutoffs, 0.30 Angstrom :sup:`-1` k-point spacing,
six fixed central atoms, one Bohr A100, one MPI rank, and a global limit of
one active job.

The SnO termination must use a metallic electronic treatment.  The source
study reports a SnO-surface band crossing the Fermi level, while O and Sn2O
are insulating; see `Kuefner et al., Physical Review B 86, 075320 (2012)
<https://link.aps.org/doi/10.1103/PhysRevB.86.075320>`_.  Accordingly, the
campaign qualified ``occupations='smearing'``, Marzari-Vanderbilt smearing,
``degauss=0.01 Ry``, and ``nbnd=280`` for SnO and alpha-Sn.  The 0.01- and
0.02-Ry probes differed by only ``0.000118779 eV Angstrom^-2`` in the common
surface-energy normalization, and a 300-band requalification differed by
``1.639725e-10 eV Angstrom^-2``.  This follows QE's documented metallic
occupation guidance; see `INPUT_PW <https://www.quantum-espresso.org/Doc/INPUT_PW.html>`_.

The resulting PBE/SSSP chemical-potential domain is
``-2.2906400569 <= Delta mu_O <= 0.0 eV``.  Within this explicitly limited
O/SnO/Sn2O clean-termination set, Sn2O is stable from ``-2.2906400569`` to
``-1.8221539181 eV`` and O is stable from that crossover to the O-rich limit.
The converged metallic SnO line is never lowest in this domain.  At the
O-rich limit the surface energies are 1.28661, 3.01428, and 4.02227 J m
:sup:`-2` for O, SnO, and Sn2O, respectively.

This is a 0 K electronic-energy result for the stated PBE/SSSP model and
candidate set, not a global experimental phase diagram: vibrational, thermal,
pressure, adsorbate, and reconstruction contributions were intentionally not
included.  The psteros/Tessera campaign records the exact AiiDA provenance and
the plot-ready analysis receipt; VASP remains a supported secondary psteros
backend, but was not used in this QE-first reference workflow.

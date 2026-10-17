=========
Changelog
=========

Complete development history and feature releases for PS-TEROS.

Recent Updates
==============

2.0.0 (unreleased)
------------------

**PS-TEROS v2**

* Public API ``import psteros``: VASP as the central engine, Quantum ESPRESSO
  as the second one, typed recipes, ``build_surface_workgraph`` and
  ``build_relax_static_workgraph``
* The VASP adapter passes the INCAR in aiida-vasp's ``incar`` namespace and
  fixes atoms with selective dynamics (``CalculationOverride(fixed_sites=...)``)
* No machine-specific defaults: ``ExecutionPolicy`` assumes no computer, queue
  or scheduler directive, and every recipe takes the user's own code label
* Surface phase diagrams of any binary or ternary compound
  (``BinaryReferences``, ``TernaryReferences``); the oxide classes are kept
* Charge-neutral symmetric terminations for semiconductors and insulators
  (``find_charge_neutral_terminations``) and their calculation set in the v2
  workflow (``ChargeNeutralSurfaceStudy``)
* Absolute surface energies of polar surfaces with pseudo-hydrogen
  passivation (``find_polar_terminations``, ``PolarSurfaceStudy``)
* Precision levels for VASP recipes (``VaspCalculationConfig(..., precision="fast" | "moderate" | "precise")``):
  the level fills in the cutoff, the k-point spacing, the convergence tags and the POTCARs that the recipe leaves out,
  and explicit values win; ``incar`` is now optional (leaving it out selects ``"moderate"``). aiida-vasp has a single
  VASP protocol, ``balanced``, which is the ``moderate`` level; ``fast`` and ``precise`` are defined by PS-TEROS, with
  PBE and no spin polarisation unless the structure is magnetic. Recipes that do not use it are unchanged
  (see :ref:`vasp-precision`)
* Opt-in vibrational free energies (harmonic, Γ-point finite differences):
  ``build_vibrations_workgraph``, ``read_vibrations``, ``HarmonicVibrations``,
  ``solid_free_energy_ev`` and ``molecule_reference_energy_ev``; slabs with a
  frozen centre use a partial Hessian, bulks a Γ-only supercell. Without them
  every result is unchanged (see :doc:`/vibrations`)
* Notes in :doc:`/vibrations` on imaginary modes, other compounds (ternary oxides) and running a long graph (queue limits, one daemon worker, recovery)
* A mode task of ``build_vibrations_workgraph`` that starts with some retrieved folders missing now raises a ``ValueError`` that
  names them, instead of a bare ``KeyError``
* CP2K backend (``pip install '.[cp2k]'``): ``Cp2kCalculationConfig``, ``build_cp2k_surface_workgraph``,
  ``build_cp2k_relax_static_workgraph`` and the readers ``read_cp2k_results`` (eV), ``read_cp2k_forces`` (eV/Å) and
  ``read_cp2k_dos``, on aiida-cp2k's ``Cp2kBaseWorkChain`` (see :doc:`/cp2k-workflow`). The QE/VASP builders and
  ``build_vibrations_workgraph`` now raise ``ValueError`` for a CP2K recipe; QE and VASP graphs are unchanged
* See :doc:`/phase-diagram` and ``CHANGE.md`` for details

Earlier Unreleased Features
---------------------------

**Relaxed VASP structures keep their kind names**

* ``build_relax_static_workgraph`` failed at the static step of every structure with a custom kind, such as the
  pseudo-hydrogens ``H0p75`` and ``H1p25`` of polar surfaces (``KeyError: 'H'``): aiida-vasp returns the relaxed
  structure with element names only, so the static calculation had no POTCAR for ``H``.
  A ``restore_vasp_kinds`` step now gives each relaxed atom the kind of the starting atom it came from. It is added
  only for structures with custom kinds; graphs of all other structures are unchanged. No results change.
* aiida-vasp also restarts a work chain (after a walltime or a failed line search) from that kind-less relaxed
  structure and failed with ``KeyError: 'H'``. The potential mapping now also carries the element name of a custom kind
  whose POTCAR is unambiguous (one pseudo-hydrogen kind per structure); with ``H0p75`` and ``H1p25`` in one slab a
  restart still stops with the error. No results change.

**VASP k-point spacing now matches its documented unit**

* ``VaspCalculationConfig.kpoints_spacing`` is documented in Å⁻¹ with the 2π of VASP's ``KSPACING``, but was passed
  unchanged to aiida-vasp, which multiplies it by 2π. A value such as 0.2 Å⁻¹ therefore gave a single Γ k-point
  for any cell of about 5 Å or more (GaAs, ``a = 5.75`` Å: ``1×1×1`` instead of ``6×6×6``) and no error.
  The adapter now converts the value. VASP calculations run with an earlier version used a much coarser mesh
  than intended and should be repeated.

**VASP k-point spacing now matches its documented unit**

* ``VaspCalculationConfig.kpoints_spacing`` is documented in Å⁻¹ with the 2π of VASP's ``KSPACING``, but was passed
  unchanged to aiida-vasp, which multiplies it by 2π. A value such as 0.2 Å⁻¹ therefore gave a single Γ k-point
  for any cell of about 5 Å or more (GaAs, ``a = 5.75`` Å: ``1×1×1`` instead of ``6×6×6``) and no error.
  The adapter now converts the value. VASP calculations run with an earlier version used a much coarser mesh
  than intended and should be repeated.

**Free energies of the reference systems (VASP, additive)**

* New ``Relax``, ``Static`` and ``Vibrations`` blocks, ``ReferenceSystem`` and ``build_vasp_reference_workgraph`` run
  relax → static → vibrations for bulks and molecules; ``reference_results`` and ``reference_thermochemistry`` read
  them back, and ``IdealGasMolecule``, ``HarmonicSolid``, ``FreeEnergy``, ``free_energies``, ``delta_mu_oxygen_ev`` and
  ``oxygen_pressure_bar`` turn them into free energies and read Δμ\ :sub:`O` as T and p(O\ :sub:`2`). Nothing existing
  changes; see :ref:`reference-thermochemistry` and ``examples/vasp_sno2_phase_diagram_vibrations/``.
* ``build_surface_workgraph`` with ``backend="vasp"`` now passes a flat INCAR to aiida-vasp under its ``incar``
  namespace and merges override INCAR tags, ``kpoints_distance`` and ``settings`` into it (results change for the VASP
  path only: before, aiida-vasp 5 rejected the flat INCAR at run time, and override ``kpoints_distance`` and ``settings``
  were ignored; OUTCAR and vasprun.xml are now kept in ``<label>_retrieved``). ``ExecutionPolicy`` has the optional
  ``extra_options``.
* Run on Lovelace with VASP 6.5.1, the new blocks needed three fixes, all in this release: the ``Vibrations`` block runs
  without band parallelisation (``ISYM = 0`` for a gas, ``NCORE = 1`` for a solid; VASP otherwise stops with "requested a
  change of the k-point set"), and the supercell of a solid lists its atoms grouped by element (VASP otherwise found no
  symmetry: 432 instead of 8 displacements for the 72-atom SnO\ :sub:`2` cell). The reference blocks also convert
  ``kpoints_spacing`` to aiida-vasp's unit like the surface builder.

**CP2K Calculator Support for AIMD**

* CP2K integration for ab initio molecular dynamics simulations
* Efficient Born-Oppenheimer MD with GPW/GAPW methods
* Automatic BASIS_MOLOPT and GTH_POTENTIALS file generation
* Seamless workflow: VASP for bulk/slab → CP2K for AIMD
* Sequential stage support with restart chaining
* See :doc:`/workflows/aimd-molecular-dynamics` for usage

**Fixed Atoms Constraints**

* Atomic constraint support for both VASP and CP2K calculators
* Flexible constraint specification: bottom, top, or center regions
* Element-specific and component-specific fixing (XYZ, XY, Z)
* Automatic constraint calculation for auto-generated slabs
* Manual control for advanced use cases
* See :doc:`/api/fixed_atoms` for API details

**Electronic Properties Module (DOS & Band Structure)**

* Comprehensive electronic structure calculations for bulk and slabs
* Material-agnostic builders for DOS and band structure
* Integration with vasp.v2.bands workchain and seekpath
* See `CHANGE.md <https://github.com/your-repo/PS-TEROS/blob/main/CHANGE.md#unreleased---electronic-properties-module-dos--band-structure>`_ for details

**AIMD Module (Ab Initio Molecular Dynamics)**

* Sequential AIMD stages with automatic restart chaining
* Parallel MD on multiple slab terminations
* Temperature ramping and multi-phase protocols
* See `CHANGE.md <https://github.com/your-repo/PS-TEROS/blob/main/CHANGE.md#unreleased---aimd-module-ab-initio-molecular-dynamics>`_ for details

**Relaxation Energy Calculation**

* Optional E_relaxed - E_unrelaxed calculation for slabs
* Quantifies energetic stabilization from surface relaxation
* See `CHANGE.md <https://github.com/your-repo/PS-TEROS/blob/main/CHANGE.md#unreleased---relaxation-energy-module>`_ for details

Version History
===============

The complete changelog with detailed API changes, implementation notes, and usage examples is maintained in ``CHANGE.md`` at the repository root.

**View the full changelog**: `CHANGE.md on GitHub <https://github.com/your-repo/PS-TEROS/blob/main/CHANGE.md>`_

Major Releases
--------------

**v1.0.0** - AiiDA/WorkGraph Modernization

* Updated to latest AiiDA and AiiDA-WorkGraph
* Added cleavage energy module
* Added restart functionality
* Added default builders
* Manual termination input support

**v0.2.0** - User-Provided Slab Structures

* Support for custom slab structures via ``input_slabs``
* Bypass automatic termination generation
* Full backward compatibility

**v0.1.2** - Bug Fixes

* Fixed enthalpy of formation calculation
* Corrected oxygen chemical potential limits

**v0.1.1** - Defect Analysis

* Added DEFECT_TYPES file generation
* Stoichiometric deviation analysis

Migration Guides
================

**Upgrading from v0.x to v1.0**

v1.0.0 requires updated AiiDA and AiiDA-WorkGraph versions. See the `migration guide <https://github.com/your-repo/PS-TEROS/blob/main/docs/migrations/v1.0.md>`_ for details.

**Adding New Features to Existing Workflows**

New features (electronic properties, AIMD, relaxation energy) are fully backward compatible. Enable them by setting the appropriate boolean flags. See :doc:`workflows/intermediate-with-features` for examples.

Contributing
============

Found a bug or have a feature request? Please open an issue on `GitHub <https://github.com/your-repo/PS-TEROS/issues>`_.

See Also
========

* :doc:`contributing` - Contributing guidelines
* :doc:`workflows/index` - Working examples showcasing features

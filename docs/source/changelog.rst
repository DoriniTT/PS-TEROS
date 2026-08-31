=========
Changelog
=========

Complete development history and feature releases for PS-TEROS.

Recent Updates
==============

Unreleased Features
-------------------

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

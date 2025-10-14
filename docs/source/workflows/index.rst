=========
Workflows
=========

PS-TEROS workflows are complete, production-ready examples organized by complexity. Each workflow builds on previous concepts and demonstrates realistic use cases.

Workflow Overview
=================

Choose your workflow based on your experience level and project requirements:

Beginner: Surface Energy Calculation
-------------------------------------

**Best for**: First-time users, basic surface energy calculations

**What you'll learn**:

* Setting up a complete PS-TEROS workflow
* Bulk and reference structure relaxation
* Formation enthalpy calculation
* Slab generation and relaxation
* Surface energy calculation for binary oxides

**Material**: Ag₂O (binary oxide)

**Runtime**: ~2-4 hours on typical cluster

:doc:`Start the Beginner Workflow → <beginner-surface-energy>`

Intermediate: Adding Features
------------------------------

**Best for**: Users comfortable with basics, ready for advanced features

**Builds on**: Beginner workflow

**New features**:

* Relaxation energy calculation (E_relaxed - E_unrelaxed)
* Cleavage energy for complementary terminations
* Electronic properties (DOS and band structure) for bulk

**Shows**: How to enable optional features with boolean flags

**Material**: Ag₂O with extended features

**Runtime**: ~4-6 hours

:doc:`Start the Intermediate Workflow → <intermediate-with-features>`

Advanced: Complete Workflow
----------------------------

**Best for**: Production calculations, publication-quality results

**Includes**: All PS-TEROS features

**Capabilities**:

* Ternary oxide support (Ag₃PO₄)
* Surface thermodynamics with chemical potential sampling
* Electronic properties for both bulk and selected slabs
* Per-slab configuration for heterogeneous terminations
* Complete provenance tracking

**Material**: Ag₃PO₄ (ternary oxide)

**Runtime**: ~8-12 hours

:doc:`Start the Advanced Workflow → <advanced-complete>`

AIMD: Molecular Dynamics
-------------------------

**Best for**: Finite-temperature studies, dynamic surface behavior

**Prerequisites**: Understanding of basic workflow

**Capabilities**:

* Sequential AIMD stages (equilibration → production)
* Automatic restart chaining between stages
* Temperature ramping protocols
* Parallel MD on multiple slab terminations
* Trajectory analysis

**Materials**: Ag₂O or Ag₃PO₄ slabs

**Runtime**: Variable (depends on MD length)

:doc:`Start the AIMD Workflow → <aimd-molecular-dynamics>`

Choosing Your Workflow
=======================

**If you're brand new**: Start with :doc:`Beginner <beginner-surface-energy>`

**If you know the basics**: Try :doc:`Intermediate <intermediate-with-features>`

**If you're publishing**: Use :doc:`Advanced <advanced-complete>`

**If you need MD**: Go to :doc:`AIMD <aimd-molecular-dynamics>`

All workflows are **fully working examples** that you can adapt for your own systems.

Need Help?
==========

* Stuck? Check the :doc:`How-To Guides </how-to/index>`
* Need API details? See the :doc:`API Reference </api/index>`
* Questions about theory? Read the :doc:`Theory </theory>` page

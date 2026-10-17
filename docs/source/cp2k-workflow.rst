.. _cp2k-workflow:

============================================
Prepare a CP2K relaxation-to-static graph
============================================

psteros also runs CP2K (Quickstep, GPW) through aiida-cp2k, installed with
``pip install '.[cp2k]'``. The CP2K builders are separate functions,
``psteros.build_cp2k_relax_static_workgraph`` and
``psteros.build_cp2k_surface_workgraph``; the QE and VASP builders refuse a
CP2K recipe. This guide builds a graph for a metallic oxide, RuO\ :sub:`2`, and
does not submit work.

Before you start
----------------

You need an AiiDA profile with a registered ``cp2k`` code. The CP2K graphs
use only aiida-cp2k's ``Cp2kBaseWorkChain``, so the AiiDA daemon needs
aiida-cp2k but not psteros. The labels, basis sets and numerical settings
below are placeholders: converge the cutoff, the basis and the k-points for
your material before using the energies.
In the snippets, ``bulk`` and ``slab`` are pymatgen structures of RuO\ :sub:`2`
and of a symmetric RuO\ :sub:`2`\ (110) slab.

Write the recipes
-----------------

``Cp2kCalculationConfig.parameters`` is the CP2K input as the nested mapping
of aiida-cp2k, with upper-case section names. A metal needs diagonalisation
with smearing and added orbitals (the OT minimiser does not support smearing).
psteros adds what depends on the structure: one ``KIND`` per kind with the
names of ``basis_sets`` and ``potentials``, the k-point mesh, the fixed atoms
of a relaxation, and, for ``ENERGY_FORCE``, the force and density-of-states
files that the readers use.

.. code-block:: python

   import psteros

   execution = psteros.ExecutionPolicy(
       queue="my-queue",
       resources={"num_machines": 1, "num_mpiprocs_per_machine": 128},
       max_wallclock_seconds=4 * 3600,
       max_concurrent_jobs=3,
   )

   dft = {
       "QS": {"METHOD": "GPW", "EPS_DEFAULT": 1.0e-12},
       "MGRID": {"NGRIDS": 4, "CUTOFF": 600, "REL_CUTOFF": 60},
       "XC": {"XC_FUNCTIONAL": {"_": "PBE"}},
       "POISSON": {"PERIODIC": "XYZ", "PSOLVER": "PERIODIC"},
       "SCF": {
           "EPS_SCF": 1.0e-7,
           "MAX_SCF": 300,
           "ADDED_MOS": 40,
           "DIAGONALIZATION": {"ALGORITHM": "STANDARD"},
           "MIXING": {"METHOD": "BROYDEN_MIXING", "ALPHA": 0.2, "NBUFFER": 8},
           "SMEAR": {"METHOD": "FERMI_DIRAC", "ELECTRONIC_TEMPERATURE": 300},
       },
   }

   def cp2k_recipe(name, run_type, motion=None, overrides=None):
       parameters = {"GLOBAL": {"RUN_TYPE": run_type}, "FORCE_EVAL": {"METHOD": "Quickstep", "DFT": dft}}
       if motion:
           parameters["MOTION"] = motion
       return psteros.SurfaceWorkflowConfig(
           backend="cp2k",
           name=name,
           execution=execution,
           role_overrides=overrides or {},
           calculation=psteros.Cp2kCalculationConfig(
               code_label="cp2k@my-cluster",
               parameters=parameters,
               basis_sets={"Ru": "DZVP-MOLOPT-PBE-GTH-q16", "O": "DZVP-MOLOPT-PBE-GTH-q6"},
               potentials={"Ru": "GTH-PBE-q16", "O": "GTH-PBE-q6"},
               basis_set_file="BASIS_MOLOPT_UZH",
               potential_file="POTENTIAL_UZH",
               kpoints_spacing=0.20,
           ),
       )

The O\ :sub:`2` reference is spin-polarised and runs at Γ; a slab keeps its
central layers fixed during the relaxation (``FIXED_ATOMS``, one-based in the
CP2K input, zero-based in ``fixed_sites``). ``CalculationOverride.parameters``
updates the recipe section by section:

.. code-block:: python

   o2 = psteros.CalculationOverride(
       kpoints_distance=100.0,
       parameters={"FORCE_EVAL": {"DFT": {"UKS": True, "MULTIPLICITY": 3}}},
   )
   fixed = psteros.CalculationOverride(fixed_sites=tuple(psteros.central_sites(slab)))
   # 0.01 eV/A is 1.94e-4 Ha/bohr; the bulk also relaxes its cell
   cell = psteros.CalculationOverride(
       parameters={"GLOBAL": {"RUN_TYPE": "CELL_OPT"}, "MOTION": {"CELL_OPT": {"MAX_FORCE": 1.94e-4}}})

   relax = cp2k_recipe("ruo2_110", "GEO_OPT", {"GEO_OPT": {"OPTIMIZER": "BFGS", "MAX_FORCE": 1.94e-4}},
                       overrides={"o2": o2, "slab": fixed, "bulk": cell})
   static = cp2k_recipe("ruo2_110_static", "ENERGY_FORCE", overrides={"o2": o2})

``kpoints_spacing`` is in Å\ :sup:`-1` with the 2π, as for VASP: the
Monkhorst-Pack mesh has ``ceil(|b_i| / kpoints_spacing)`` points along each
reciprocal vector. A 1×1×1 mesh (or ``kpoints_spacing=None``) runs at Γ only,
without a ``KPOINTS`` section; give a molecule a large
``CalculationOverride(kpoints_distance=...)``. ``basis_set_file`` and
``potential_file`` are read from CP2K's data directory, or from
``input_files`` (file name -> local path), which copies them into every job
and keeps them in the provenance.

The builder checks every label before it builds anything: each relaxation
needs ``GLOBAL/RUN_TYPE`` ``GEO_OPT`` or ``CELL_OPT`` (recipe plus override),
each static calculation ``ENERGY_FORCE`` or ``ENERGY``. Kind names must be the
element followed by digits, if any, which is how aiida-cp2k writes them.

Build and read the graph
------------------------

.. code-block:: python

   graph = psteros.build_cp2k_relax_static_workgraph(
       {"bulk": bulk, "slab": slab, "o2": psteros.triplet_o2_cell()},
       relax,
       static,
   )

Each label gets ``<label>_relax_cp2k`` and ``<label>_static_cp2k`` tasks and
the graph outputs ``<label>_relaxed_structure``,
``<label>_relax_parameters``, ``<label>_static_parameters`` and the
``_retrieved`` folders of both stages. After the run:

.. code-block:: python

   from aiida import orm

   node = orm.load_node(graph_pk)
   energies, structures = psteros.read_cp2k_results(node, ["bulk", "slab", "o2"])
   forces = psteros.read_cp2k_forces(node, ["slab"])  # eV/A, one row per atom
   dos = psteros.read_cp2k_dos(node, ["bulk"])["bulk"]
   print(dos.partially_occupied_bins())  # > 0: partly filled states, a metal

Energies are in eV (CP2K's total ``FORCE_EVAL`` energy times
27.211386245988). With smearing that energy is the electronic free energy,
whose derivatives are the forces. The force file's unit is read from its
header. ``psteros.build_cp2k_surface_workgraph`` runs one calculation per label
(the recipe's run type, ``ENERGY_FORCE`` by default) with the outputs
``<label>_parameters``, ``<label>_structure`` and ``<label>_retrieved``, which
the same readers accept.

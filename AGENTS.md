# PS-TEROS

Ab initio surface thermodynamics with AiiDA WorkGraphs (Quantum ESPRESSO and VASP backends).
See README.md for the user workflow (structures -> WorkGraph -> phase diagram).

- Public API: everything exported in `psteros/__init__.py` (`__all__`). Typed recipes in `psteros/config.py`,
  WorkGraph builders in `psteros/workflow.py`, engine adapters in `psteros/backends/`, pure-Python analysis in
  `psteros/thermodynamics.py`, `psteros/phase_diagram.py`, `psteros/phase_diagram_ternary.py`.
- `psteros/core/` and `psteros/experimental/` are legacy code: do not build new features on them; reuse ideas, not imports.
- Tests: `python -m pytest tests/unit` (pure Python, no AiiDA profile); WorkGraph tests build graphs with `submit=False`.

## New features only add; nothing that works today changes

Every new feature (vibrational contributions, new calculation types, new analyses, ...) builds on top of the
current program. A script, notebook or finished AiiDA graph that works today must give the same result after the
feature lands.

- **No breaking changes to the public API.** Do not rename, remove or reorder public functions, classes,
  arguments, dataclass fields, WorkGraph output names (`{label}_static_parameters`, `{label}_relaxed_structure`, ...),
  CSV columns or figure defaults. Do not change a default value or a unit.
- **New behaviour is opt-in.** Add it as a new keyword argument with a default that reproduces today's behaviour
  exactly (e.g. `temperature=None` means "0 K total energies, as before"), a new optional dataclass field placed
  after the existing ones, or a new function/class. Leaving the new option out must give bit-identical numbers.
- **Prove it with a test.** Each feature adds a test showing that the old call path still returns the old result,
  next to the tests of the new behaviour. Existing tests are not edited to make them pass.
- **Leave working builders alone.** The existing builders in `psteros/workflow.py` (`build_surface_workgraph`,
  the relax -> static builders, ...) keep their signature, graph and outputs. A more general builder is a new
  function next to them, not a rewrite of them.
- **Mixed inputs are an error, not a guess.** If an option must apply to every structure to be consistent
  (e.g. vibrations given for some terminations but not others), raise a `ValueError` that names the offending
  labels, in line with the strict recipe harmony of the rest of psteros.

## Building blocks (one consistent way to add a calculation)

New calculation steps are composable blocks, so users learn one pattern (idea taken from kiln's stages/tiles):

- **One block = one step**, a frozen dataclass holding a `SurfaceWorkflowConfig` plus its own typed options
  (e.g. `Relax`, `Static`, `Vibrations` with a `VibrationsConfig`). Validate everything in `__post_init__` or
  before the graph is built; error messages name the block and the structure label.
- **Blocks are linked by name**: a block takes the structure (and, if useful, the restart folder) of an earlier
  block. A remote folder that a later block restarts from is never cleaned.
- **Calcfunctions only for run-time values** (displaced structures of a relaxed slab, frequencies from forces);
  anything known at build time is computed directly. Wire all tasks at build time; no dynamic sub-graphs.
- **Outputs are namespaced and additive**: a block adds new outputs (`{label}_vibrations`, ...) and never changes
  existing ones. Provide a reader (like `read_vibrations(pk)`) that turns them into pure-Python objects and works
  while a graph runs or after it fails.
- **Physics stays pure Python.** Each feature has an AiiDA-free layer (dataclasses + functions in `psteros/`)
  that can be tested and used with energies typed by hand; the AiiDA block only produces its inputs.
- **Backend parity**: write the backend-neutral version first (works with every adapter in `psteros/backends/`);
  engine-specific shortcuts (VASP `IBRION=5`, QE `ph.x`) come later as options of the same block.

## User-facing consistency

- Same style as the existing API: typed frozen dataclasses, keyword-only arguments, explicit units in names
  (`_ev`, `_angstrom2`, `_cm1`, `temperature` in K), `submit=False` by default.
- Every new public name goes into `psteros/__init__.py` and `__all__`, a docstring with the formula and its
  convention, a section in `docs/source/`, and a short example in `examples/`.
- Record the user-visible addition in `docs/source/changelog.rst`.

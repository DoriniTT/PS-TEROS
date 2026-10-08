# PS-TEROS: guidelines for contributors and agents

PS-TEROS builds AiiDA WorkGraphs for oxide surface thermodynamics and analyses
their energies in pure Python. People run real campaigns with it, so every new
feature **builds on** the program: it never breaks, renames or silently changes
what already works.

## The public API is a contract

The public API is everything exported from `psteros/__init__.py` (`__all__`)
and the documented signatures in `docs/source/api.rst`.

- **Add, don't change.** A new feature goes in new modules, new functions and
  new optional keyword arguments whose defaults reproduce today's behaviour.
  Never remove or rename a public name, reorder positional arguments, change a
  default, or change the meaning or units of an existing field, return value
  or WorkGraph output name.
- **The same inputs give the same graph and the same numbers.** An existing
  script must build the same WorkGraph (task names, inputs, outputs) and an
  existing analysis must return the same values after your change. New physics
  (corrections, extra terms) is opt-in, never switched on by default.
- **Bug fixes that change results are not silent.** If existing behaviour is
  wrong, fix it in its own commit, note it in `CHANGE.md`, and say in the
  docstring what changed and from which version.
- **Deprecate before removing.** Keep the old name working, emit a
  `DeprecationWarning` that names the replacement, and document it.
- `psteros/core/`, `psteros/compat.py` and `psteros/experimental/` are the
  legacy (pre-1.0) code. Don't build new features on them and don't change
  them, except to fix a bug that blocks their users.

## How a feature fits in

Follow the shape of the existing modules so that every feature is used the
same way:

1. **Typed recipe:** frozen dataclasses with no AiiDA nodes, validated in
   `__post_init__` with an error that names the bad field and the allowed
   values (see `psteros/config.py`). A recipe can be built, printed and
   tested without an AiiDA profile.
2. **Graph builder:** `build_<something>_workgraph(structures, recipe, *, submit=False)`
   returns the WorkGraph and only submits it when asked. Import AiiDA,
   aiida-workgraph and plugins inside the function, never at module level of a
   pure-Python module, so `import psteros` works without a profile.
   Name tasks and outputs with stable, documented labels (`<label>_<stage>`)
   so results can be found by name.
3. **Blocks for multi-step calculations:** a calculation step (relax, static,
   vibrations, ...) is a block, like the tiles in kiln: one module that
   validates its recipe and adds its tasks to the graph, returning its output
   sockets, the structure for the next block and its remote folder. Blocks are
   linked by name (`structure_from="relax"`), so a new step is a new block,
   not a new branch in an existing builder.
4. **Pure-Python analysis:** functions or frozen dataclasses that take plain
   numbers (eV, Å, K, bar) and return plain numbers or a result object with
   `to_csv`/`plot` where useful (see `psteros/phase_diagram.py`). Analysis
   never needs AiiDA; it also accepts values read from any other source.
5. **Units in names:** `_ev`, `_angstrom2`, `_k`, `_bar`, `_cm1` suffixes on
   fields and arguments; state the convention (for example
   `mu_O = E(O2)/2 + Delta mu_O`) in the docstring.
6. **Show the parts:** when a quantity is a sum of physical terms (energies,
   free energies, corrections), return the breakdown as well as the total, so
   every number is auditable.

Keep the public entry points few and obvious: export new user-facing objects
from `psteros/__init__.py` and add them to `__all__`; keep helpers private
(`_name`) or in their module.

## Backends

Quantum ESPRESSO (`aiida-quantumespresso`) is the primary backend and VASP
(`aiida-vasp`, `vasp.v2.vasp`) the secondary one. A feature may support one
backend first; it must then reject the other with a clear error, never fall
back silently. Backend adapters live in `psteros/backends/`.

## Tests and docs

- Every change runs `python -m pytest tests/unit` green; add tests for every
  new public object, including its validation errors.
- Pure-Python physics is checked against known values (textbook or literature
  numbers, or an independent implementation such as ASE) with the source cited
  in the test.
- WorkGraph builders get construction tests that build the graph with
  `submit=False` on the throwaway profile from `tests/conftest.py`; nothing is
  ever submitted from the tests.
- Document new public objects in `docs/source/api.rst`, add a how-to or
  example under `examples/` for a new workflow, and add an entry to
  `CHANGE.md`.

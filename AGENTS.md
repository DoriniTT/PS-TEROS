# psteros contributor guidance

The public package is `psteros`; do not add compatibility imports for the
retired package name or reintroduce a generic composition API.

Quantum ESPRESSO through `aiida-quantumespresso` is the primary backend.  VASP
is a maintained secondary backend.  Keep backend-specific code in
`psteros/backends/`, keep configuration free of AiiDA nodes, and keep
thermodynamic analysis pure Python where possible.

The public entry point is `psteros.build_surface_workgraph`.  It must preserve
the explicit execution policy and never create an unbounded fan-out.  The
maintained SnO2(110) campaign runs on Bohr `gpu_a100` with one active job.

Run focused tests after a change, then the complete suite when the affected
surface includes legacy modules.  Real calculations are controlled by the
external Tessera Odyssey Auto project, not by an in-package watcher.


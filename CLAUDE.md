# CLAUDE.md

This file is a thin mirror of `/home/trevizam/git/PS-TEROS/AGENTS.md`.

For this repository:

- Read `/home/trevizam/git/PS-TEROS/AGENTS.md` first. It is the canonical instruction file.
- The main entry point is `teros.core.workgraph.build_core_workgraph()`.
- Restart the AiiDA daemon after Python code changes.
- Run the narrowest relevant tiered tests before broader suites.
- If you touch shared lego logic, check whether `quantum-lego` also needs the same update.

For the full workflow architecture, presets, lego module guidance, and testing details, use:

- `/home/trevizam/git/PS-TEROS/AGENTS.md`

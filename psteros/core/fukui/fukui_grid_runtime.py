"""Helpers to locate the external FukuiGrid runtime."""

from __future__ import annotations

import os
from pathlib import Path


_MODULE_DIR = Path(__file__).resolve().parent
DEFAULT_FUKUI_GRID_DIR = _MODULE_DIR.parent.parent / 'external' / 'FukuiGrid'


def _normalize_import_root(candidate: Path) -> Path | None:
    """Return a sys.path entry that makes `import FukuiGrid` work."""
    candidate = candidate.expanduser()

    if candidate.is_file():
        if candidate.name == 'FukuiGrid.py':
            return candidate.parent
        return None

    if (candidate / 'FukuiGrid.py').is_file():
        return candidate

    nested_module = candidate / 'FukuiGrid' / 'FukuiGrid.py'
    if nested_module.is_file():
        return nested_module.parent

    return None


def get_fukui_grid_search_paths() -> list[Path]:
    """Return candidate locations in priority order."""
    candidates: list[Path] = []

    env_path = os.environ.get('FUKUI_GRID_PATH')
    if env_path:
        candidates.append(Path(env_path))

    candidates.extend(
        [
            DEFAULT_FUKUI_GRID_DIR,
            Path.home() / 'programs',
            Path.home() / 'programs' / 'FukuiGrid',
            Path.home() / 'programs' / 'FukuiGrid.py',
        ]
    )

    unique: list[Path] = []
    seen: set[Path] = set()
    for candidate in candidates:
        expanded = candidate.expanduser()
        if expanded in seen:
            continue
        seen.add(expanded)
        unique.append(expanded)

    return unique


def resolve_fukui_grid_import_root() -> Path | None:
    """Find a sys.path entry that exposes the FukuiGrid module."""
    for candidate in get_fukui_grid_search_paths():
        import_root = _normalize_import_root(candidate)
        if import_root is not None:
            return import_root
    return None


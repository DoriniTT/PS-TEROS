#!/home/thiagotd/envs/aiida/bin/python
r"""
Unified surface modifier for ASE slabs: vacancies, H decoration, and combined modes.

Modes
-----
1) vacancies : remove surface O atoms in all non-empty combinations.
2) hydrogen  : add H on surface O atoms in all non-empty combinations.
3) combine   : for each subset V of surface O (to remove), add H to an independent
               subset H of the remaining surface O (S \ V); enumerate all pairs (V, H).
               By default, skip variants where both V and H are empty.

Surfaces
--------
- top, bottom, or both. If 'both', operations are generated independently per face
  and combined via cartesian products.

Assumptions
-----------
- Slab normal is aligned with +z axis.
- H is placed at O ± (0,0,oh_dist) depending on surface (top:+, bottom:-).

Examples
--------
# 1) Vacancies on top
python surface_modes.py POSCAR --mode vacancies --z-window 0.5 --outdir vac_top --fmt vasp

# 2) Hydrogenation (OH) on bottom with 1.0 Å O–H
python surface_modes.py POSCAR --mode hydrogen --which-surface bottom --oh-dist 1.00 --outdir H_bottom

# 3) Combined (vacancies + H on remaining sites) on both surfaces
python surface_modes.py POSCAR --mode combine --which-surface both --outdir combo_both

# 4) Complete mode - runs all three modes
python surface_modes.py POSCAR --mode complete --which-surface top

# Include unmodified baseline variant as well
python surface_modes.py POSCAR --mode hydrogen --include-empty --outdir H_with_ref
"""
from __future__ import annotations
import argparse
import json
from dataclasses import dataclass
from itertools import combinations, product
from pathlib import Path
from typing import List, Sequence, Tuple, Dict, Any

import numpy as np
from ase.io import read, write
from ase import Atoms, Atom

# ----------------------------- supercell -----------------------------

def make_supercell(atoms: Atoms, nx: int, ny: int, nz: int) -> Atoms:
    """Create supercell by repeating unit cell.

    Args:
        atoms: Input structure
        nx, ny, nz: Repetitions along each axis

    Returns:
        Supercell structure
    """
    return atoms.repeat((nx, ny, nz))

def calculate_surface_area(atoms: Atoms) -> float:
    """Calculate xy-plane surface area of the cell.

    Args:
        atoms: Structure with cell vectors

    Returns:
        Surface area in nm²
    """
    cell = atoms.get_cell()
    # Cross product of cell vectors [0] and [1] gives area vector
    # Magnitude of this vector is the xy-plane area
    cross = np.cross(cell[0], cell[1])
    area_angstrom2 = np.linalg.norm(cross)
    # Convert Ų to nm² (1 nm² = 100 Ų)
    area_nm2 = area_angstrom2 / 100.0
    return area_nm2

def calculate_coverage(n_species: int, area_nm2: float, precision: int = 4) -> float:
    """Calculate surface coverage in species per nm².

    Args:
        n_species: Number of species (vacancies or H atoms)
        area_nm2: Surface area in nm²
        precision: Decimal places for rounding

    Returns:
        Coverage in species/nm²
    """
    if area_nm2 == 0:
        return 0.0
    coverage = n_species / area_nm2
    return round(coverage, precision)

def generate_variant_name(mode: str, index: int, coverage: float | Tuple[float, float]) -> str:
    """Generate simplified variant filename based on mode and coverage.

    Format:
        - vacancies: vac_{index:03d}_{coverage}
        - hydrogen: oh_{index:03d}_{coverage}
        - combine: combo_{index:03d}_{vac_cov}_{oh_cov}

    Args:
        mode: One of 'vacancies', 'hydrogen', 'combine'
        index: Sequential variant index
        coverage: Single coverage value or tuple of (vac_coverage, oh_coverage) for combine

    Returns:
        Simplified filename (without extension)
    """
    mode_map = {
        'vacancies': 'vac',
        'hydrogen': 'oh',
        'combine': 'combo'
    }

    mode_prefix = mode_map.get(mode, mode)

    if mode == 'combine' and isinstance(coverage, tuple):
        # Dual coverage: combo_001_1.5_0.75
        vac_cov, oh_cov = coverage
        return f"{mode_prefix}_{index:03d}_{vac_cov}_{oh_cov}"
    else:
        # Single coverage: vac_001_0.25 or oh_005_1.75
        return f"{mode_prefix}_{index:03d}_{coverage}"

def deduplicate_by_coverage(
    variants_list: List[Tuple[Dict[str, Any], Atoms, Any]],
    mode: str
) -> Tuple[List[Tuple[Dict[str, Any], Atoms, Any]], Dict[str, int]]:
    """Deduplicate variants by coverage, keeping first occurrence.

    Args:
        variants_list: List of (variant_data, atoms_object, coverage_key) tuples
        mode: Mode name for reporting

    Returns:
        (deduplicated_list, stats_dict)
    """
    seen_coverages = {}
    deduplicated = []

    for variant in variants_list:
        variant_data, atoms_obj, coverage_key = variant
        if coverage_key not in seen_coverages:
            seen_coverages[coverage_key] = True
            deduplicated.append(variant)

    stats = {
        "total_generated": len(variants_list),
        "unique_coverages": len(seen_coverages),
        "kept": len(deduplicated)
    }

    return deduplicated, stats


def sample_by_coverage_bins(
    variants_list: List[Tuple[Dict[str, Any], Atoms, Any]],
    n_bins: int | None,
    mode: str
) -> List[Tuple[Dict[str, Any], Atoms, Any]]:
    """Sample variants by dividing coverage range into bins.

    Args:
        variants_list: List of (variant_data, atoms_object, coverage_key) tuples
        n_bins: Number of bins to divide coverage range into (None = no binning)
        mode: Mode name ('vacancies', 'hydrogen', 'combine')

    Returns:
        Sampled variants list
    """
    if n_bins is None or n_bins <= 0:
        return variants_list

    if len(variants_list) == 0:
        return variants_list

    # Extract coverages
    coverages = [variant[2] for variant in variants_list]

    # Check if this is 2D coverage (tuple)
    if isinstance(coverages[0], tuple):
        # 2D binning for combined mode
        vac_coverages = [c[0] for c in coverages]
        oh_coverages = [c[1] for c in coverages]

        vac_min, vac_max = min(vac_coverages), max(vac_coverages)
        oh_min, oh_max = min(oh_coverages), max(oh_coverages)

        # Handle edge cases
        if vac_max == vac_min:
            vac_width = 1.0
        else:
            vac_width = (vac_max - vac_min) / n_bins

        if oh_max == oh_min:
            oh_width = 1.0
        else:
            oh_width = (oh_max - oh_min) / n_bins

        # Create 2D grid of bin centers
        sampled = []
        sampled_indices = set()

        for i in range(n_bins):
            for j in range(n_bins):
                vac_center = vac_min + (i + 0.5) * vac_width
                oh_center = oh_min + (j + 0.5) * oh_width

                # Find variant closest to this grid cell center (Euclidean distance)
                min_dist = float('inf')
                closest_idx = -1

                for idx, variant in enumerate(variants_list):
                    if idx in sampled_indices:
                        continue
                    vac_cov, oh_cov = variant[2]
                    dist = ((vac_cov - vac_center)**2 + (oh_cov - oh_center)**2)**0.5
                    if dist < min_dist:
                        min_dist = dist
                        closest_idx = idx

                if closest_idx >= 0:
                    sampled.append(variants_list[closest_idx])
                    sampled_indices.add(closest_idx)

        return sampled

    # 1D binning for single coverage
    min_cov = min(coverages)
    max_cov = max(coverages)

    # If all same coverage or only one variant, return as-is
    if max_cov == min_cov or len(variants_list) <= n_bins:
        return variants_list

    bin_width = (max_cov - min_cov) / n_bins
    bin_centers = [min_cov + (i + 0.5) * bin_width for i in range(n_bins)]

    # For each bin, select structure closest to center
    sampled = []
    sampled_indices = set()

    for center in bin_centers:
        # Find variant with coverage closest to this bin center
        min_dist = float('inf')
        closest_idx = -1

        for idx, variant in enumerate(variants_list):
            if idx in sampled_indices:
                continue
            coverage = variant[2]
            dist = abs(coverage - center)
            if dist < min_dist:
                min_dist = dist
                closest_idx = idx

        if closest_idx >= 0:
            sampled.append(variants_list[closest_idx])
            sampled_indices.add(closest_idx)

    return sampled


def add_bin_metadata(
    variants_list: List[Tuple[Dict[str, Any], Atoms, Any]],
    n_bins: int,
    mode: str
) -> List[Tuple[Dict[str, Any], Atoms, Any]]:
    """Add bin_id and bin_center metadata to each variant.

    Args:
        variants_list: List of sampled variants
        n_bins: Number of bins used
        mode: Mode name

    Returns:
        Variants with metadata added
    """
    if len(variants_list) == 0:
        return variants_list

    coverages = [v[2] for v in variants_list]

    # Check if 2D coverage
    if isinstance(coverages[0], tuple):
        # 2D metadata
        vac_coverages = [c[0] for c in coverages]
        oh_coverages = [c[1] for c in coverages]

        vac_min, vac_max = min(vac_coverages), max(vac_coverages)
        oh_min, oh_max = min(oh_coverages), max(oh_coverages)

        vac_width = (vac_max - vac_min) / n_bins if vac_max > vac_min else 1.0
        oh_width = (oh_max - oh_min) / n_bins if oh_max > oh_min else 1.0

        result = []
        for variant_data, atoms_obj, coverage_key in variants_list:
            vac_cov, oh_cov = coverage_key

            # Determine which bin this belongs to
            vac_bin = int((vac_cov - vac_min) / vac_width) if vac_width > 0 else 0
            vac_bin = min(vac_bin, n_bins - 1)
            oh_bin = int((oh_cov - oh_min) / oh_width) if oh_width > 0 else 0
            oh_bin = min(oh_bin, n_bins - 1)

            vac_center = vac_min + (vac_bin + 0.5) * vac_width
            oh_center = oh_min + (oh_bin + 0.5) * oh_width

            variant_data['bin_id_vac'] = vac_bin
            variant_data['bin_id_oh'] = oh_bin
            variant_data['bin_center_vac'] = round(vac_center, 4)
            variant_data['bin_center_oh'] = round(oh_center, 4)

            result.append((variant_data, atoms_obj, coverage_key))

        return result
    else:
        # 1D metadata
        min_cov = min(coverages)
        max_cov = max(coverages)
        bin_width = (max_cov - min_cov) / n_bins if max_cov > min_cov else 1.0

        result = []
        for variant_data, atoms_obj, coverage_key in variants_list:
            # Determine which bin this belongs to
            bin_id = int((coverage_key - min_cov) / bin_width) if bin_width > 0 else 0
            bin_id = min(bin_id, n_bins - 1)

            bin_center = min_cov + (bin_id + 0.5) * bin_width

            variant_data['bin_id'] = bin_id
            variant_data['bin_center'] = round(bin_center, 4)

            result.append((variant_data, atoms_obj, coverage_key))

        return result


# ----------------------------- helpers -----------------------------

@dataclass
class SurfaceSelection:
    species: str
    indices: List[int]
    z_coords: List[float]
    which: str  # "top" or "bottom"

def select_surface_atoms(
    atoms: Atoms,
    species: str = "O",
    z_window: float = 0.5,
    which_surface: str = "top",
) -> SurfaceSelection:
    """Return indices of `species` on the chosen surface by z-window criterion."""
    assert which_surface in {"top", "bottom"}
    idx_species = [i for i, a in enumerate(atoms) if atoms[i].symbol == species]
    if not idx_species:
        raise ValueError(f"No atoms of species '{species}' found.")

    z = atoms.get_positions()[:, 2]
    z_species = np.array([z[i] for i in idx_species])

    if which_surface == "top":
        ref = z_species.max()
        mask = z_species >= (ref - z_window)
        which = "top"
    else:
        ref = z_species.min()
        mask = z_species <= (ref + z_window)
        which = "bottom"

    sel_idx_species = np.array(idx_species)[mask].tolist()
    sel_z = z_species[mask].tolist()

    if not sel_idx_species:
        raise ValueError(
            f"No '{species}' atoms found within z_window={z_window} Å on {which_surface} surface."
        )

    return SurfaceSelection(
        species=species,
        indices=sel_idx_species,
        z_coords=sel_z,
        which=which,
    )

def all_subsets(items: Sequence[int], allow_empty: bool) -> List[Tuple[int, ...]]:
    """Power set in increasing size order; include empty if allow_empty."""
    out: List[Tuple[int, ...]] = []
    n = len(items)
    start = 0 if allow_empty else 1
    for r in range(start, n + 1):
        out.extend(combinations(items, r))
    return out

def delete_indices_copy(atoms: Atoms, idx_to_delete: Sequence[int]) -> Atoms:
    new = atoms.copy()
    for i in sorted(idx_to_delete, reverse=True):
        del new[i]
    return new

def add_H_on_indices(atoms: Atoms, indices: Sequence[int], which_surface: str, oh_dist: float) -> Atoms:
    assert which_surface in {"top", "bottom"}
    sign = +1.0 if which_surface == "top" else -1.0
    delta = np.array([0.0, 0.0, sign * oh_dist])
    new = atoms.copy()
    pos = new.get_positions()
    for i in indices:
        h_xyz = pos[i] + delta
        new.append(Atom('H', position=h_xyz))
    return new

def add_H_at_positions(atoms: Atoms, positions: Sequence[Tuple[float, float, float]], which_surface: str, oh_dist: float) -> Atoms:
    """Add H atoms at specified O positions (uses absolute coords, not indices)."""
    assert which_surface in {"top", "bottom"}
    sign = +1.0 if which_surface == "top" else -1.0
    delta = np.array([0.0, 0.0, sign * oh_dist])
    new = atoms.copy()
    for o_xyz in positions:
        h_xyz = np.array(o_xyz) + delta
        new.append(Atom('H', position=h_xyz))
    return new

def cartesian_with_empty(a_sets: List[Tuple[int, ...]], b_sets: List[Tuple[int, ...]], skip_both_empty=True):
    a_with_empty = [tuple()] + a_sets
    b_with_empty = [tuple()] + b_sets
    for a, b in product(a_with_empty, b_with_empty):
        if skip_both_empty and (not a and not b):
            continue
        yield a, b

# ----------------------------- core class -----------------------------

class SurfaceModifier:
    def __init__(
        self,
        atoms: Atoms,
        species: str = "O",
        z_window: float = 0.5,
        which_surface: str = "top",
        oh_dist: float = 0.98,
        include_empty: bool = False,
        prefix: str = "surf",
        outdir: Path | str = "outputs",
        fmt: str = "vasp",
        tag_original: bool = False,
        supercell: Tuple[int, int, int] | None = None,
        deduplicate_by_coverage: bool = False,
        coverage_precision: int = 4,
        coverage_bins: int | None = None,
        max_vacancies: int | None = None,
        max_OH: int | None = None,
    ):
        assert which_surface in {"top", "bottom", "both"}
        assert fmt in {"vasp", "cif", "xyz"}

        # Create supercell if requested
        if supercell is not None:
            atoms = make_supercell(atoms, *supercell)

        self.atoms = atoms
        self.species = species
        self.z_window = z_window
        self.which_surface = which_surface
        self.oh_dist = oh_dist
        self.include_empty = include_empty
        self.prefix = prefix
        self.outdir = Path(outdir)
        self.fmt = fmt
        self.tag_original = tag_original
        self.deduplicate_by_coverage = deduplicate_by_coverage
        self.coverage_precision = coverage_precision
        self.coverage_bins = coverage_bins
        self.max_vacancies = max_vacancies
        self.max_OH = max_OH

        # Calculate and store surface area for coverage calculations
        self.surface_area_nm2 = calculate_surface_area(atoms)

        self.outdir.mkdir(parents=True, exist_ok=True)
        if tag_original:
            base = self.outdir / f"000_original.{ 'vasp' if fmt=='vasp' else fmt }"
            write(base.as_posix(), atoms, direct=True) if fmt == "vasp" else write(base.as_posix(), atoms)

    # ---------- selection ----------
    def _select(self):
        if self.which_surface in {"top", "bottom"}:
            sel = select_surface_atoms(self.atoms, self.species, self.z_window, self.which_surface)
            return {"single": sel}
        else:
            top = select_surface_atoms(self.atoms, self.species, self.z_window, "top")
            bot = select_surface_atoms(self.atoms, self.species, self.z_window, "bottom")
            return {"top": top, "bottom": bot}

    # ---------- writing ----------
    def _write(self, atoms: Atoms, name: str):
        path = self.outdir / f"{name}.{ 'vasp' if self.fmt=='vasp' else self.fmt }"
        if self.fmt == "vasp":
            write(path.as_posix(), atoms, direct=True)
        else:
            write(path.as_posix(), atoms)
        return str(path)

    # ---------- filtering ----------
    def _filter_by_max_modifications(
        self,
        variants: List[Tuple[Dict[str, Any], Atoms, Any]],
        mode: str
    ) -> Tuple[List[Tuple[Dict[str, Any], Atoms, Any]], int]:
        """Filter variants by maximum number of modifications.

        Args:
            variants: List of (variant_data, atoms_object, coverage_key) tuples
            mode: Mode name ('vacancies', 'hydrogen', 'combine')

        Returns:
            (filtered_variants, count_filtered)
        """
        # Check if filtering is needed
        if mode == 'vacancies' and self.max_vacancies is None:
            return variants, 0
        if mode == 'hydrogen' and self.max_OH is None:
            return variants, 0
        if mode == 'combine' and self.max_vacancies is None and self.max_OH is None:
            return variants, 0

        original_count = len(variants)
        filtered = []

        for v in variants:
            meta = v[0]

            if mode == 'vacancies':
                if meta['count_deleted'] <= self.max_vacancies:
                    filtered.append(v)

            elif mode == 'hydrogen':
                if meta['count_H'] <= self.max_OH:
                    filtered.append(v)

            elif mode == 'combine':
                vac_ok = self.max_vacancies is None or meta['count_deleted'] <= self.max_vacancies
                oh_ok = self.max_OH is None or meta['count_H'] <= self.max_OH
                if vac_ok and oh_ok:
                    filtered.append(v)

        filtered_count = original_count - len(filtered)
        return filtered, filtered_count

    # ---------- modes ----------
    def run_vacancies(self) -> Dict[str, Any]:
        manifest = {"mode": "vacancies", "variants": []}

        # Generate all variants in memory first
        variants_to_process = []

        if self.which_surface in {"top", "bottom"}:
            sel = self._select()["single"]
            sets = all_subsets(sel.indices, allow_empty=self.include_empty)
            sets = sorted(sets, key=lambda t: (len(t), t))
            for subset in sets:
                if not subset and not self.include_empty:
                    continue
                var = delete_indices_copy(self.atoms, subset)

                # Calculate coverage
                coverage = calculate_coverage(len(subset), self.surface_area_nm2, self.coverage_precision)

                variant_data = {
                    "deleted_indices": list(subset),
                    "surface": sel.which,
                    "count_deleted": len(subset),
                    "surface_area_nm2": self.surface_area_nm2,
                    "vacancy_coverage": coverage
                }
                variants_to_process.append((variant_data, var, coverage))
        else:
            top = self._select()["top"]
            bot = self._select()["bottom"]
            top_sets = all_subsets(top.indices, allow_empty=self.include_empty)
            bot_sets = all_subsets(bot.indices, allow_empty=self.include_empty)
            for tset, bset in cartesian_with_empty(top_sets[1:] if not self.include_empty else top_sets,
                                                   bot_sets[1:] if not self.include_empty else bot_sets,
                                                   skip_both_empty=not self.include_empty):
                delete_union = sorted(set(tset).union(bset))
                var = delete_indices_copy(self.atoms, delete_union)

                # Calculate coverage
                coverage = calculate_coverage(len(delete_union), self.surface_area_nm2, self.coverage_precision)

                variant_data = {
                    "deleted_top": list(tset),
                    "deleted_bottom": list(bset),
                    "count_deleted": len(delete_union),
                    "surface": "both",
                    "surface_area_nm2": self.surface_area_nm2,
                    "vacancy_coverage": coverage
                }
                variants_to_process.append((variant_data, var, coverage))

        # Filter by max_vacancies if specified (BEFORE deduplication)
        total_before_filter = len(variants_to_process)
        variants_to_process, filtered_count = self._filter_by_max_modifications(
            variants_to_process, 'vacancies'
        )
        if filtered_count > 0:
            manifest["filter_stats"] = {
                "max_vacancies": self.max_vacancies,
                "max_OH": self.max_OH,
                "total_generated": total_before_filter,
                "filtered_out": filtered_count,
                "kept": len(variants_to_process)
            }

        # Deduplicate if requested
        if self.deduplicate_by_coverage and len(variants_to_process) > 0:
            variants_to_process, stats = deduplicate_by_coverage(variants_to_process, "vacancies")
            manifest["deduplication_stats"] = stats

        # Apply binning if requested
        if self.coverage_bins is not None and self.deduplicate_by_coverage:
            before_binning = len(variants_to_process)
            variants_to_process = sample_by_coverage_bins(
                variants_to_process,
                self.coverage_bins,
                'vacancies'
            )
            # Add bin metadata
            variants_to_process = add_bin_metadata(
                variants_to_process,
                self.coverage_bins,
                'vacancies'
            )
            manifest["binning_stats"] = {
                "bins_requested": self.coverage_bins,
                "before_binning": before_binning,
                "after_binning": len(variants_to_process)
            }

        # Write files for kept variants with simplified names
        for idx, (variant_data, var, coverage) in enumerate(variants_to_process):
            # Generate simplified filename based on coverage
            simple_name = generate_variant_name('vacancies', idx, coverage)
            path = self._write(var, simple_name)
            variant_data["name"] = simple_name
            variant_data["file"] = path
            manifest["variants"].append(variant_data)

        return manifest

    def run_hydrogen(self) -> Dict[str, Any]:
        manifest = {"mode": "hydrogen", "variants": []}

        # Generate all variants in memory first
        variants_to_process = []

        if self.which_surface in {"top", "bottom"}:
            sel = self._select()["single"]
            sets = all_subsets(sel.indices, allow_empty=self.include_empty)
            sets = sorted(sets, key=lambda t: (len(t), t))
            for subset in sets:
                if not subset and not self.include_empty:
                    continue
                var = self.atoms.copy() if not subset else add_H_on_indices(self.atoms, subset, sel.which, self.oh_dist)

                # Calculate coverage
                coverage = calculate_coverage(len(subset), self.surface_area_nm2, self.coverage_precision)

                variant_data = {
                    "decorated_indices": list(subset),
                    "surface": sel.which,
                    "count_H": len(subset),
                    "surface_area_nm2": self.surface_area_nm2,
                    "OH_coverage": coverage
                }
                variants_to_process.append((variant_data, var, coverage))
        else:
            top = self._select()["top"]
            bot = self._select()["bottom"]
            top_sets = all_subsets(top.indices, allow_empty=self.include_empty)
            bot_sets = all_subsets(bot.indices, allow_empty=self.include_empty)
            for tset, bset in cartesian_with_empty(top_sets[1:] if not self.include_empty else top_sets,
                                                   bot_sets[1:] if not self.include_empty else bot_sets,
                                                   skip_both_empty=not self.include_empty):
                var = self.atoms.copy()
                if tset:
                    var = add_H_on_indices(var, tset, "top", self.oh_dist)
                if bset:
                    var = add_H_on_indices(var, bset, "bottom", self.oh_dist)

                # Calculate coverage
                coverage = calculate_coverage(len(tset) + len(bset), self.surface_area_nm2, self.coverage_precision)

                variant_data = {
                    "decorated_top": list(tset),
                    "decorated_bottom": list(bset),
                    "count_H": len(tset) + len(bset),
                    "surface": "both",
                    "surface_area_nm2": self.surface_area_nm2,
                    "OH_coverage": coverage
                }
                variants_to_process.append((variant_data, var, coverage))

        # Filter by max_OH if specified (BEFORE deduplication)
        total_before_filter = len(variants_to_process)
        variants_to_process, filtered_count = self._filter_by_max_modifications(
            variants_to_process, 'hydrogen'
        )
        if filtered_count > 0:
            manifest["filter_stats"] = {
                "max_vacancies": self.max_vacancies,
                "max_OH": self.max_OH,
                "total_generated": total_before_filter,
                "filtered_out": filtered_count,
                "kept": len(variants_to_process)
            }

        # Deduplicate if requested
        if self.deduplicate_by_coverage and len(variants_to_process) > 0:
            variants_to_process, stats = deduplicate_by_coverage(variants_to_process, "hydrogen")
            manifest["deduplication_stats"] = stats

        # Apply binning if requested
        if self.coverage_bins is not None and self.deduplicate_by_coverage:
            before_binning = len(variants_to_process)
            variants_to_process = sample_by_coverage_bins(
                variants_to_process,
                self.coverage_bins,
                'hydrogen'
            )
            # Add bin metadata
            variants_to_process = add_bin_metadata(
                variants_to_process,
                self.coverage_bins,
                'hydrogen'
            )
            manifest["binning_stats"] = {
                "bins_requested": self.coverage_bins,
                "before_binning": before_binning,
                "after_binning": len(variants_to_process)
            }

        # Write files for kept variants with simplified names
        for idx, (variant_data, var, coverage) in enumerate(variants_to_process):
            # Generate simplified filename based on coverage
            simple_name = generate_variant_name('hydrogen', idx, coverage)
            path = self._write(var, simple_name)
            variant_data["name"] = simple_name
            variant_data["file"] = path
            manifest["variants"].append(variant_data)

        return manifest

    def run_combine(self) -> Dict[str, Any]:
        """
        Combined enumeration:
        For each subset V (vacancies) of surface O, choose an independent subset H of (S \\ V) to decorate.
        """
        manifest = {"mode": "combine", "variants": []}

        # Generate all variants in memory first
        variants_to_process = []

        def enumerate_single(face_sel: SurfaceSelection):
            S = face_sel.indices
            V_sets = all_subsets(S, allow_empty=self.include_empty)
            for V in sorted(V_sets, key=lambda t: (len(t), t)):
                remaining = [i for i in S if i not in V]
                H_sets = all_subsets(remaining, allow_empty=self.include_empty)
                for H in sorted(H_sets, key=lambda t: (len(t), t)):
                    if not (V or H) and not self.include_empty:
                        continue
                    # Capturar posições dos O que receberão H ANTES de deletar
                    h_positions = [tuple(self.atoms.positions[i]) for i in H] if H else []
                    # Deletar vacâncias
                    var = delete_indices_copy(self.atoms, V) if V else self.atoms.copy()
                    # Adicionar H usando as posições capturadas
                    if h_positions:
                        var = add_H_at_positions(var, h_positions, face_sel.which, self.oh_dist)

                    # Calculate coverages
                    vac_coverage = calculate_coverage(len(V), self.surface_area_nm2, self.coverage_precision)
                    oh_coverage = calculate_coverage(len(H), self.surface_area_nm2, self.coverage_precision)
                    coverage_key = (vac_coverage, oh_coverage)

                    variant_data = {
                        "surface": face_sel.which,
                        "deleted_indices": list(V), "decorated_indices": list(H),
                        "count_deleted": len(V), "count_H": len(H),
                        "surface_area_nm2": self.surface_area_nm2,
                        "vacancy_coverage": vac_coverage,
                        "OH_coverage": oh_coverage
                    }
                    variants_to_process.append((variant_data, var, coverage_key))

        if self.which_surface in {"top", "bottom"}:
            sel = self._select()["single"]
            enumerate_single(sel)
        else:
            top = self._select()["top"]
            bot = self._select()["bottom"]

            # Build cartesian products of (V_top, H_top) x (V_bot, H_bot)
            # Generate all pairs on each face, then combine.
            def generate_face_pairs(sel):
                S = sel.indices
                V_sets = all_subsets(S, allow_empty=self.include_empty)
                pairs = []
                for V in V_sets:
                    remaining = [i for i in S if i not in V]
                    H_sets = all_subsets(remaining, allow_empty=self.include_empty)
                    for H in H_sets:
                        if not (V or H) and not self.include_empty:
                            continue
                        pairs.append((V, H))
                # deterministic order
                pairs = sorted(pairs, key=lambda vh: (len(vh[0]), vh[0], len(vh[1]), vh[1]))
                return pairs

            top_pairs = generate_face_pairs(top)
            bot_pairs = generate_face_pairs(bot)

            for (Vt, Ht), (Vb, Hb) in product(top_pairs, bot_pairs):
                if not (Vt or Ht or Vb or Hb) and not self.include_empty:
                    continue
                # Capturar posições dos O que receberão H ANTES de deletar
                ht_positions = [tuple(self.atoms.positions[i]) for i in Ht] if Ht else []
                hb_positions = [tuple(self.atoms.positions[i]) for i in Hb] if Hb else []
                # Deletar vacâncias
                var = self.atoms.copy()
                if Vt or Vb:
                    delete_union = sorted(set(Vt).union(Vb))
                    var = delete_indices_copy(var, delete_union)
                # Adicionar H usando as posições capturadas
                if ht_positions:
                    var = add_H_at_positions(var, ht_positions, "top", self.oh_dist)
                if hb_positions:
                    var = add_H_at_positions(var, hb_positions, "bottom", self.oh_dist)

                # Calculate coverages for both surfaces
                total_vac = len(Vt) + len(Vb)
                total_h = len(Ht) + len(Hb)
                vac_coverage = calculate_coverage(total_vac, self.surface_area_nm2, self.coverage_precision)
                oh_coverage = calculate_coverage(total_h, self.surface_area_nm2, self.coverage_precision)
                coverage_key = (vac_coverage, oh_coverage)

                variant_data = {
                    "surface": "both",
                    "deleted_top": list(Vt), "decorated_top": list(Ht),
                    "deleted_bottom": list(Vb), "decorated_bottom": list(Hb),
                    "count_deleted": total_vac, "count_H": total_h,
                    "surface_area_nm2": self.surface_area_nm2,
                    "vacancy_coverage": vac_coverage,
                    "OH_coverage": oh_coverage
                }
                variants_to_process.append((variant_data, var, coverage_key))

        # Filter by max_vacancies and max_OH if specified (BEFORE deduplication)
        total_before_filter = len(variants_to_process)
        variants_to_process, filtered_count = self._filter_by_max_modifications(
            variants_to_process, 'combine'
        )
        if filtered_count > 0:
            manifest["filter_stats"] = {
                "max_vacancies": self.max_vacancies,
                "max_OH": self.max_OH,
                "total_generated": total_before_filter,
                "filtered_out": filtered_count,
                "kept": len(variants_to_process)
            }

        # Deduplicate if requested
        if self.deduplicate_by_coverage and len(variants_to_process) > 0:
            variants_to_process, stats = deduplicate_by_coverage(variants_to_process, "combine")
            manifest["deduplication_stats"] = stats

        # Apply binning if requested (2D binning for dual coverage)
        if self.coverage_bins is not None and self.deduplicate_by_coverage:
            before_binning = len(variants_to_process)
            variants_to_process = sample_by_coverage_bins(
                variants_to_process,
                self.coverage_bins,
                'combine'
            )
            # Add bin metadata
            variants_to_process = add_bin_metadata(
                variants_to_process,
                self.coverage_bins,
                'combine'
            )
            manifest["binning_stats"] = {
                "bins_requested": self.coverage_bins,
                "before_binning": before_binning,
                "after_binning": len(variants_to_process),
                "binning_type": "2D_grid"
            }

        # Write files for kept variants with simplified names
        for idx, (variant_data, var, coverage) in enumerate(variants_to_process):
            # Generate simplified filename based on dual coverage
            simple_name = generate_variant_name('combine', idx, coverage)
            path = self._write(var, simple_name)
            variant_data["name"] = simple_name
            variant_data["file"] = path
            manifest["variants"].append(variant_data)

        return manifest

    def run_complete(self) -> Dict[str, Any]:
        """
        Execute all three modes (vacancies, hydrogen, combine) with all files in the same directory.
        """
        # All modes write to the same directory
        # Run vacancies mode
        manifest_vac = self.run_vacancies()
        with open(self.outdir / "manifest_vacancies.json", "w") as f:
            json.dump(manifest_vac, f, indent=2)
        
        # Run hydrogen mode
        manifest_h = self.run_hydrogen()
        with open(self.outdir / "manifest_hydrogen.json", "w") as f:
            json.dump(manifest_h, f, indent=2)
        
        # Run combine mode
        manifest_combo = self.run_combine()
        with open(self.outdir / "manifest_combine.json", "w") as f:
            json.dump(manifest_combo, f, indent=2)
        
        # Create combined manifest
        combined_manifest = {
            "mode": "complete",
            "total_variants": len(manifest_vac["variants"]) + len(manifest_h["variants"]) + len(manifest_combo["variants"]),
            "modes": {
                "vacancies": {
                    "count": len(manifest_vac["variants"]),
                    "manifest_file": "manifest_vacancies.json"
                },
                "hydrogen": {
                    "count": len(manifest_h["variants"]),
                    "manifest_file": "manifest_hydrogen.json"
                },
                "combine": {
                    "count": len(manifest_combo["variants"]),
                    "manifest_file": "manifest_combine.json"
                }
            }
        }

        # Add deduplication stats if present
        if self.deduplicate_by_coverage:
            combined_manifest["deduplication_summary"] = {
                "vacancies": manifest_vac.get("deduplication_stats", {}),
                "hydrogen": manifest_h.get("deduplication_stats", {}),
                "combine": manifest_combo.get("deduplication_stats", {})
            }

        return combined_manifest

# ----------------------------- CLI -----------------------------

def main():
    p = argparse.ArgumentParser(description="Unified surface modifier: vacancies, hydrogenation, combined (ASE).")
    p.add_argument("input", help="Input structure (ASE-readable: POSCAR/CONTCAR, CIF, etc.)")
    p.add_argument("--mode", choices=["vacancies", "hydrogen", "combine", "complete"], required=True,
                   help="Which modification mode to run. 'complete' runs all three modes.")
    p.add_argument("--species", default="O", help="Target species (default: O)")
    p.add_argument("--z-window", type=float, default=0.5, help="Å window from extreme z to include as surface (default: 0.5)")
    p.add_argument("--which-surface", choices=["top", "bottom", "both"], default="top",
                   help="Which surface to act on (default: top)")
    p.add_argument("--oh-dist", type=float, default=0.98, help="O–H distance in Å placed along ±z (default: 0.98)")
    p.add_argument("--outdir", default="surf_variants", help="Output directory (default: surf_variants)")
    p.add_argument("--fmt", default="vasp", choices=["vasp", "cif", "xyz"], help="Output format (default: vasp)")
    p.add_argument("--tag-original", action="store_true", help="Write a copy of the original structure as 000_original.*")
    p.add_argument("--prefix", default="surf", help="Filename prefix (default: surf)")
    p.add_argument("--include-empty", action="store_true",
                   help="Include baseline (no-op) variants where applicable.")
    p.add_argument("--no-manifest", action="store_true", help="Do not write JSON manifest")

    # New arguments for supercell and coverage
    p.add_argument("--supercell", type=int, nargs=3, metavar=("NX", "NY", "NZ"),
                   help="Create NX×NY×NZ supercell before modifications (e.g., --supercell 2 2 1)")
    p.add_argument("--deduplicate-by-coverage", action="store_true",
                   help="Keep only one variant per unique coverage value")
    p.add_argument("--coverage-precision", type=int, default=4,
                   help="Decimal places for coverage rounding (default: 4)")
    p.add_argument("--coverage-bins", type=int, default=None,
                   help="Number of bins for coverage sampling (requires --deduplicate-by-coverage)")
    p.add_argument("--max-vacancies", type=int, default=None,
                   help="Maximum number of O vacancies per structure (None = unlimited)")
    p.add_argument("--max-OH", type=int, default=None,
                   help="Maximum number of OH groups per structure (None = unlimited)")

    args = p.parse_args()

    # Validate coverage-bins argument
    if args.coverage_bins is not None:
        if not args.deduplicate_by_coverage:
            p.error("--coverage-bins requires --deduplicate-by-coverage to be enabled")
        if args.coverage_bins <= 0:
            p.error("--coverage-bins must be a positive integer")

    # Adjust default outdir for complete mode
    if args.mode == "complete" and args.outdir == "surf_variants":
        args.outdir = "all_possibilities"
    
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    atoms: Atoms = read(args.input)

    sm = SurfaceModifier(
        atoms=atoms,
        species=args.species,
        z_window=args.z_window,
        which_surface=args.which_surface,
        oh_dist=args.oh_dist,
        include_empty=args.include_empty,
        prefix=args.prefix,
        outdir=outdir,
        fmt=args.fmt,
        tag_original=args.tag_original,
        supercell=tuple(args.supercell) if args.supercell else None,
        deduplicate_by_coverage=args.deduplicate_by_coverage,
        coverage_precision=args.coverage_precision,
        coverage_bins=args.coverage_bins,
        max_vacancies=args.max_vacancies,
        max_OH=args.max_OH,
    )

    if args.mode == "vacancies":
        manifest = sm.run_vacancies()
    elif args.mode == "hydrogen":
        manifest = sm.run_hydrogen()
    elif args.mode == "combine":
        manifest = sm.run_combine()
    else:  # complete
        manifest = sm.run_complete()

    if not args.no_manifest:
        with open(outdir / "manifest.json", "w") as f:
            json.dump(manifest, f, indent=2)

    if args.mode == "complete":
        print(f"Mode: {args.mode}. Wrote {manifest['total_variants']} variants total to '{outdir}'.")
        print(f"  - Vacancies: {manifest['modes']['vacancies']['count']} variants")
        print(f"  - Hydrogen: {manifest['modes']['hydrogen']['count']} variants")
        print(f"  - Combine: {manifest['modes']['combine']['count']} variants")

        # Print deduplication stats for complete mode if present
        if args.deduplicate_by_coverage:
            # Load individual manifests to check for deduplication stats
            for mode_name in ["vacancies", "hydrogen", "combine"]:
                manifest_file = outdir / f"manifest_{mode_name}.json"
                if manifest_file.exists():
                    with open(manifest_file, "r") as f:
                        mode_manifest = json.load(f)
                    if "deduplication_stats" in mode_manifest:
                        stats = mode_manifest["deduplication_stats"]
                        print(f"  - {mode_name.capitalize()} deduplication: {stats['total_generated']} -> {stats['kept']} variants")
                    if "binning_stats" in mode_manifest:
                        stats = mode_manifest["binning_stats"]
                        print(f"  - {mode_name.capitalize()} binning: {stats['before_binning']} -> {stats['after_binning']} ({stats['bins_requested']} bins)")
    else:
        print(f"Mode: {args.mode}. Wrote {len(manifest['variants'])} variants to '{outdir}'.")

        # Print deduplication stats for single mode if present
        if "deduplication_stats" in manifest:
            stats = manifest["deduplication_stats"]
            removed = stats['total_generated'] - stats['kept']
            print(f"Deduplication: {stats['total_generated']} generated -> {stats['kept']} kept (removed {removed} duplicates)")

        # Print binning stats if present
        if "binning_stats" in manifest:
            stats = manifest["binning_stats"]
            print(f"Binning: {stats['before_binning']} unique -> {stats['after_binning']} sampled ({stats['bins_requested']} bins)")

    if not args.no_manifest:
        print(f"Manifest: {outdir / 'manifest.json'}")

if __name__ == "__main__":
    main()
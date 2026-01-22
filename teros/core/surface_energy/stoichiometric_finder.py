"""
Stoichiometric and Symmetric Surface Finder.

This module provides experimental functionality to find surfaces that are BOTH
stoichiometric (preserve bulk composition) AND symmetric (equivalent top/bottom).

For simple metals like Au, Ag, Cu, nearly all surfaces meet these criteria.
For intermetallics and oxides, finding such surfaces is challenging but valuable
because they allow using the simple surface energy formula:

    γ = (E_slab - N·E_bulk/atom) / (2A)

instead of requiring chemical potential-dependent thermodynamics.

Strategies:
    1. filter_first: Generate all slabs with symmetrize=False, filter post-hoc
    2. symmetrize_check: Generate with symmetrize=True, verify stoichiometry preserved
    3. thickness_scan: Scan thickness range (10-30Å) to find valid terminations
    4. bond_preservation: Preserve chemical units (PO4, SO4, etc.) using bonds parameter

When no valid surface is found, a NoStoichiometricSymmetricSurfaceError is raised
with a recommendation to use the thermodynamics approach instead.
"""

from __future__ import annotations

import typing as t
import warnings
from dataclasses import dataclass, field

import numpy as np
from aiida import orm
from aiida_workgraph import task

if t.TYPE_CHECKING:
    from pymatgen.core.surface import Slab
    from pymatgen.core import Structure


# Default strategies in order of execution
DEFAULT_STRATEGIES = [
    'filter_first',      # Fastest: generate all, filter post-hoc
    'symmetrize_check',  # Quick: symmetrize=True, check if stoich preserved
    'thickness_scan',    # Thorough: scan thickness range
    'bond_preservation', # Complex materials: preserve chemical units
]


class NoStoichiometricSymmetricSurfaceError(ValueError):
    """
    Raised when no stoichiometric and symmetric surface can be found.

    This error indicates that the simple surface energy formula cannot be used
    for this Miller index, and ab initio atomistic thermodynamics with chemical
    potential dependencies should be used instead.

    Attributes:
        miller_index: The Miller index that failed
        strategies_tried: List of strategies that were attempted
        recommendation: Suggested alternative approach
    """

    def __init__(
        self,
        miller_index: tuple[int, int, int],
        strategies_tried: list[str],
        n_terminations_checked: int = 0,
        n_stoichiometric: int = 0,
        n_symmetric: int = 0,
    ):
        self.miller_index = miller_index
        self.strategies_tried = strategies_tried
        self.n_terminations_checked = n_terminations_checked
        self.n_stoichiometric = n_stoichiometric
        self.n_symmetric = n_symmetric

        message = (
            f"No stoichiometric AND symmetric surface found for Miller index {miller_index}.\n"
            f"  Strategies tried: {', '.join(strategies_tried)}\n"
            f"  Terminations checked: {n_terminations_checked}\n"
            f"  Stoichiometric only: {n_stoichiometric}\n"
            f"  Symmetric only: {n_symmetric}\n"
            f"\n"
            f"RECOMMENDATION: Use the thermodynamics module with chemical potential\n"
            f"dependencies instead:\n"
            f"    from teros.core.thermodynamics import build_thermodynamics_workgraph\n"
            f"\n"
            f"This approach calculates γ(Δμ_O, Δμ_M) for non-stoichiometric surfaces."
        )
        super().__init__(message)


@dataclass
class SlabSearchResult:
    """
    Result from a stoichiometric+symmetric slab search.

    Attributes:
        slab: The PyMatGen Slab object (None if no valid slab found)
        is_stoichiometric: Whether slab composition matches bulk
        is_symmetric: Whether top and bottom surfaces are equivalent
        strategy_used: Which search strategy found this slab
        thickness_angstrom: Slab thickness in Angstroms
        termination_index: Index of the termination (0, 1, 2, ...)
        miller_index: Miller index tuple
        bonds_broken: Number of bonds broken during slab creation
        warnings: Any warnings generated during search
    """
    slab: t.Optional['Slab']
    is_stoichiometric: bool
    is_symmetric: bool
    strategy_used: str
    thickness_angstrom: float
    termination_index: int
    miller_index: tuple[int, int, int] = (0, 0, 1)
    bonds_broken: int = 0
    warnings: list[str] = field(default_factory=list)

    @property
    def is_valid(self) -> bool:
        """Check if this result represents a valid stoichiometric+symmetric slab."""
        return self.slab is not None and self.is_stoichiometric and self.is_symmetric

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            'is_stoichiometric': self.is_stoichiometric,
            'is_symmetric': self.is_symmetric,
            'strategy_used': self.strategy_used,
            'thickness_angstrom': self.thickness_angstrom,
            'termination_index': self.termination_index,
            'miller_index': list(self.miller_index),
            'bonds_broken': self.bonds_broken,
            'warnings': self.warnings,
            'is_valid': self.is_valid,
        }


@dataclass
class MillerFeasibilityReport:
    """
    Feasibility report for finding stoichiometric+symmetric surfaces.

    This report helps users understand whether a given Miller index is
    suitable for simple surface energy calculations or requires thermodynamics.

    Attributes:
        miller_index: The Miller index analyzed
        has_valid_surfaces: Whether any valid surfaces were found
        n_terminations_checked: Total number of terminations examined
        n_stoichiometric: Number of stoichiometric terminations
        n_symmetric: Number of symmetric terminations
        n_both: Number of terminations that are BOTH
        best_thickness: Recommended thickness if valid surface found
        recommended_strategy: Which strategy worked best
        notes: Additional information and recommendations
    """
    miller_index: tuple[int, int, int]
    has_valid_surfaces: bool
    n_terminations_checked: int
    n_stoichiometric: int
    n_symmetric: int
    n_both: int
    best_thickness: t.Optional[float] = None
    recommended_strategy: str = ''
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            'miller_index': list(self.miller_index),
            'has_valid_surfaces': self.has_valid_surfaces,
            'n_terminations_checked': self.n_terminations_checked,
            'n_stoichiometric': self.n_stoichiometric,
            'n_symmetric': self.n_symmetric,
            'n_both': self.n_both,
            'best_thickness': self.best_thickness,
            'recommended_strategy': self.recommended_strategy,
            'notes': self.notes,
        }

    def __str__(self) -> str:
        status = "OK" if self.has_valid_surfaces else "NEEDS THERMODYNAMICS"
        return (
            f"Miller {self.miller_index}: {status}\n"
            f"  Checked: {self.n_terminations_checked} | "
            f"Stoich: {self.n_stoichiometric} | "
            f"Sym: {self.n_symmetric} | "
            f"Both: {self.n_both}\n"
            f"  {'Best thickness: ' + f'{self.best_thickness:.1f}Å' if self.best_thickness else ''}"
            f"  {', '.join(self.notes) if self.notes else ''}"
        )


def get_bulk_composition(structure: t.Union[orm.StructureData, 'Structure']) -> dict:
    """
    Get the reduced composition of a bulk structure.

    Args:
        structure: AiiDA StructureData or PyMatGen Structure

    Returns:
        Dictionary mapping element symbol to count (reduced formula)
    """
    from pymatgen.io.ase import AseAtomsAdaptor

    if isinstance(structure, orm.StructureData):
        ase_atoms = structure.get_ase()
        pmg_structure = AseAtomsAdaptor.get_structure(ase_atoms)
    else:
        pmg_structure = structure

    # Get reduced composition
    composition = pmg_structure.composition.reduced_composition
    return {str(el): int(amt) for el, amt in composition.items()}


def check_slab_properties(
    slab: 'Slab',
    bulk_composition: dict,
) -> tuple[bool, bool, dict]:
    """
    Check if a slab is stoichiometric and symmetric.

    Args:
        slab: PyMatGen Slab object
        bulk_composition: Bulk composition as {element: count} dict

    Returns:
        Tuple of (is_stoichiometric, is_symmetric, details_dict)
    """
    # Check symmetry using PyMatGen's built-in method
    is_symmetric = slab.is_symmetric()

    # Check stoichiometry
    slab_composition = slab.composition.reduced_composition
    slab_comp_dict = {str(el): int(amt) for el, amt in slab_composition.items()}

    # Compare compositions (they should have the same ratios)
    is_stoichiometric = (slab_comp_dict == bulk_composition)

    # If direct comparison fails, check if ratios match
    if not is_stoichiometric:
        # Normalize to smallest integer ratios
        bulk_gcd = np.gcd.reduce(list(bulk_composition.values()))
        slab_gcd = np.gcd.reduce(list(slab_comp_dict.values()))

        bulk_norm = {k: v // bulk_gcd for k, v in bulk_composition.items()}
        slab_norm = {k: v // slab_gcd for k, v in slab_comp_dict.items()}

        is_stoichiometric = (bulk_norm == slab_norm)

    details = {
        'slab_composition': slab_comp_dict,
        'bulk_composition': bulk_composition,
        'is_symmetric': is_symmetric,
        'is_stoichiometric': is_stoichiometric,
        'thickness': slab.get_orthogonal_c_slab_scale(slab) * slab.lattice.c,
    }

    return is_stoichiometric, is_symmetric, details


def detect_bonds(
    structure: t.Union[orm.StructureData, 'Structure'],
    tolerance: float = 0.2,
) -> dict[tuple[str, str], float]:
    """
    Auto-detect bonds from covalent radii.

    This is useful for complex materials like oxides where preserving
    chemical units (PO4, SO4, etc.) is important.

    Args:
        structure: Structure to analyze
        tolerance: Tolerance factor for bond detection (default 0.2 = 20%)

    Returns:
        Dictionary mapping (element1, element2) to bond length
    """
    from pymatgen.io.ase import AseAtomsAdaptor
    from pymatgen.core.periodic_table import Element

    if isinstance(structure, orm.StructureData):
        ase_atoms = structure.get_ase()
        pmg_structure = AseAtomsAdaptor.get_structure(ase_atoms)
    else:
        pmg_structure = structure

    # Get unique elements
    elements = list(set(str(site.specie) for site in pmg_structure))

    bonds = {}
    for i, el1 in enumerate(elements):
        for el2 in elements[i:]:
            # Sum of covalent radii
            r1 = Element(el1).atomic_radius or 1.5
            r2 = Element(el2).atomic_radius or 1.5
            bond_length = (r1 + r2) * (1 + tolerance)

            # Store both orderings
            bonds[(el1, el2)] = bond_length
            if el1 != el2:
                bonds[(el2, el1)] = bond_length

    return bonds


def _strategy_filter_first(
    generator,
    bulk_composition: dict,
    miller_index: tuple[int, int, int],
    min_slab_thickness: float,
    min_vacuum_thickness: float,
    lll_reduce: bool,
    center_slab: bool,
    primitive: bool,
) -> list[SlabSearchResult]:
    """
    Strategy 1: Generate all slabs with symmetrize=False, filter post-hoc.

    This is the fastest strategy as it generates slabs once and filters.
    """
    results = []

    slabs = generator.get_slabs(
        symmetrize=False,  # Don't force symmetry
        repair=True,       # Repair broken bonds if possible
    )

    for i, slab in enumerate(slabs):
        is_stoich, is_sym, details = check_slab_properties(slab, bulk_composition)

        result = SlabSearchResult(
            slab=slab if (is_stoich and is_sym) else None,
            is_stoichiometric=is_stoich,
            is_symmetric=is_sym,
            strategy_used='filter_first',
            thickness_angstrom=details['thickness'],
            termination_index=i,
            miller_index=miller_index,
        )

        if is_stoich and is_sym:
            results.append(result)

    return results


def _strategy_symmetrize_check(
    generator,
    bulk_composition: dict,
    miller_index: tuple[int, int, int],
    min_slab_thickness: float,
    min_vacuum_thickness: float,
    lll_reduce: bool,
    center_slab: bool,
    primitive: bool,
) -> list[SlabSearchResult]:
    """
    Strategy 2: Generate with symmetrize=True, verify stoichiometry preserved.

    This strategy leverages PyMatGen's symmetrization but verifies that
    the symmetrization process didn't break stoichiometry.
    """
    results = []

    slabs = generator.get_slabs(
        symmetrize=True,  # Force symmetric slabs
        repair=True,
    )

    for i, slab in enumerate(slabs):
        is_stoich, is_sym, details = check_slab_properties(slab, bulk_composition)

        # Symmetric should always be True here, but verify
        if not is_sym:
            warnings.warn(
                f"Slab {i} for {miller_index} not symmetric despite symmetrize=True"
            )
            continue

        result = SlabSearchResult(
            slab=slab if is_stoich else None,
            is_stoichiometric=is_stoich,
            is_symmetric=is_sym,
            strategy_used='symmetrize_check',
            thickness_angstrom=details['thickness'],
            termination_index=i,
            miller_index=miller_index,
        )

        if is_stoich:
            results.append(result)

    return results


def _strategy_thickness_scan(
    structure: 'Structure',
    bulk_composition: dict,
    miller_index: tuple[int, int, int],
    min_vacuum_thickness: float,
    lll_reduce: bool,
    center_slab: bool,
    primitive: bool,
    min_thickness: float = 10.0,
    max_thickness: float = 30.0,
    thickness_step: float = 2.0,
) -> list[SlabSearchResult]:
    """
    Strategy 3: Scan thickness range to find valid terminations.

    Some Miller indices only have stoichiometric+symmetric slabs at
    specific thicknesses. This strategy systematically scans.
    """
    from pymatgen.core.surface import SlabGenerator

    results = []
    seen_compositions = set()  # Avoid duplicates

    for thickness in np.arange(min_thickness, max_thickness + 0.1, thickness_step):
        generator = SlabGenerator(
            structure,
            miller_index,
            min_slab_size=thickness,
            min_vacuum_size=min_vacuum_thickness,
            lll_reduce=lll_reduce,
            center_slab=center_slab,
            primitive=primitive,
            in_unit_planes=False,
        )

        slabs = generator.get_slabs(symmetrize=False, repair=True)

        for i, slab in enumerate(slabs):
            # Create unique identifier to avoid duplicates
            comp_key = str(sorted(slab.composition.items()))
            n_atoms = len(slab)
            unique_key = f"{comp_key}_{n_atoms}"

            if unique_key in seen_compositions:
                continue
            seen_compositions.add(unique_key)

            is_stoich, is_sym, details = check_slab_properties(slab, bulk_composition)

            if is_stoich and is_sym:
                result = SlabSearchResult(
                    slab=slab,
                    is_stoichiometric=is_stoich,
                    is_symmetric=is_sym,
                    strategy_used='thickness_scan',
                    thickness_angstrom=details['thickness'],
                    termination_index=i,
                    miller_index=miller_index,
                )
                results.append(result)

    return results


def _strategy_bond_preservation(
    structure: 'Structure',
    bulk_composition: dict,
    miller_index: tuple[int, int, int],
    min_slab_thickness: float,
    min_vacuum_thickness: float,
    lll_reduce: bool,
    center_slab: bool,
    primitive: bool,
    bonds: dict = None,
    auto_detect_bonds: bool = True,
) -> list[SlabSearchResult]:
    """
    Strategy 4: Preserve chemical units using bonds parameter.

    This is useful for complex materials like phosphates (PO4), sulfates (SO4),
    etc. where breaking certain bonds would create unrealistic surfaces.
    """
    from pymatgen.core.surface import SlabGenerator

    results = []

    # Detect bonds if not provided
    if bonds is None and auto_detect_bonds:
        bonds = detect_bonds(structure)

    if bonds is None:
        return results  # Can't use this strategy without bonds

    generator = SlabGenerator(
        structure,
        miller_index,
        min_slab_size=min_slab_thickness,
        min_vacuum_size=min_vacuum_thickness,
        lll_reduce=lll_reduce,
        center_slab=center_slab,
        primitive=primitive,
        in_unit_planes=False,
    )

    try:
        # Get slabs with bond preservation (max_broken_bonds=0)
        slabs = generator.get_slabs(
            symmetrize=False,
            repair=True,
            bonds=bonds,
            max_broken_bonds=0,
        )

        for i, slab in enumerate(slabs):
            is_stoich, is_sym, details = check_slab_properties(slab, bulk_composition)

            if is_stoich and is_sym:
                result = SlabSearchResult(
                    slab=slab,
                    is_stoichiometric=is_stoich,
                    is_symmetric=is_sym,
                    strategy_used='bond_preservation',
                    thickness_angstrom=details['thickness'],
                    termination_index=i,
                    miller_index=miller_index,
                    bonds_broken=0,
                )
                results.append(result)

    except Exception as e:
        warnings.warn(f"Bond preservation strategy failed for {miller_index}: {e}")

    return results


def find_stoichiometric_symmetric_slabs(
    structure: t.Union[orm.StructureData, 'Structure'],
    miller_index: t.Union[list, tuple],
    min_slab_thickness: float = 15.0,
    min_vacuum_thickness: float = 15.0,
    lll_reduce: bool = True,
    center_slab: bool = True,
    primitive: bool = True,
    strategies: list[str] = None,
    bonds: dict = None,
    auto_detect_bonds: bool = False,
    max_thickness: float = 30.0,
    raise_if_not_found: bool = True,
) -> list[SlabSearchResult]:
    """
    Find stoichiometric AND symmetric slabs for a given Miller index.

    This function tries multiple strategies to find valid surfaces:
    1. filter_first: Generate all, filter post-hoc (fastest)
    2. symmetrize_check: Use symmetrize=True, verify stoichiometry
    3. thickness_scan: Scan thickness range 10-30Å
    4. bond_preservation: Preserve chemical units

    Args:
        structure: Bulk structure (AiiDA StructureData or PyMatGen Structure)
        miller_index: Miller index as list or tuple, e.g., [1, 1, 0]
        min_slab_thickness: Minimum slab thickness in Angstroms
        min_vacuum_thickness: Minimum vacuum thickness in Angstroms
        lll_reduce: Use LLL reduction for slab cell
        center_slab: Center slab in c direction
        primitive: Find primitive cell before slab generation
        strategies: List of strategies to try (default: all)
        bonds: Dictionary of bonds for bond_preservation strategy
        auto_detect_bonds: Auto-detect bonds from covalent radii
        max_thickness: Maximum thickness for thickness_scan strategy
        raise_if_not_found: Raise error if no valid surface found

    Returns:
        List of SlabSearchResult objects, sorted by thickness

    Raises:
        NoStoichiometricSymmetricSurfaceError: If no valid surface found
            and raise_if_not_found=True
    """
    from pymatgen.core.surface import SlabGenerator
    from pymatgen.io.ase import AseAtomsAdaptor

    # Convert structure if needed
    if isinstance(structure, orm.StructureData):
        ase_atoms = structure.get_ase()
        pmg_structure = AseAtomsAdaptor.get_structure(ase_atoms)
    else:
        pmg_structure = structure

    # Normalize miller index
    miller_tuple = tuple(int(m) for m in miller_index)

    # Get bulk composition
    bulk_composition = get_bulk_composition(pmg_structure)

    # Default strategies
    if strategies is None:
        strategies = DEFAULT_STRATEGIES.copy()

    # Create slab generator for filter_first and symmetrize_check
    generator = SlabGenerator(
        pmg_structure,
        miller_tuple,
        min_slab_size=min_slab_thickness,
        min_vacuum_size=min_vacuum_thickness,
        lll_reduce=lll_reduce,
        center_slab=center_slab,
        primitive=primitive,
        in_unit_planes=False,
    )

    all_results = []
    strategies_tried = []

    # Track statistics for error message
    n_terminations_checked = 0
    n_stoichiometric = 0
    n_symmetric = 0

    for strategy in strategies:
        strategies_tried.append(strategy)

        if strategy == 'filter_first':
            results = _strategy_filter_first(
                generator, bulk_composition, miller_tuple,
                min_slab_thickness, min_vacuum_thickness,
                lll_reduce, center_slab, primitive,
            )

        elif strategy == 'symmetrize_check':
            results = _strategy_symmetrize_check(
                generator, bulk_composition, miller_tuple,
                min_slab_thickness, min_vacuum_thickness,
                lll_reduce, center_slab, primitive,
            )

        elif strategy == 'thickness_scan':
            results = _strategy_thickness_scan(
                pmg_structure, bulk_composition, miller_tuple,
                min_vacuum_thickness, lll_reduce, center_slab, primitive,
                min_thickness=10.0, max_thickness=max_thickness,
            )

        elif strategy == 'bond_preservation':
            results = _strategy_bond_preservation(
                pmg_structure, bulk_composition, miller_tuple,
                min_slab_thickness, min_vacuum_thickness,
                lll_reduce, center_slab, primitive,
                bonds=bonds, auto_detect_bonds=auto_detect_bonds,
            )

        else:
            warnings.warn(f"Unknown strategy: {strategy}")
            continue

        all_results.extend(results)

        # Early exit if we found valid results
        if results:
            break

    # Update statistics (approximate from first strategy results)
    if strategies_tried:
        # Re-check all terminations for statistics
        all_slabs = generator.get_slabs(symmetrize=False, repair=True)
        for slab in all_slabs:
            n_terminations_checked += 1
            is_stoich, is_sym, _ = check_slab_properties(slab, bulk_composition)
            if is_stoich:
                n_stoichiometric += 1
            if is_sym:
                n_symmetric += 1

    # Sort results by thickness (thinnest first for efficiency)
    all_results.sort(key=lambda r: r.thickness_angstrom)

    # Raise error if no valid surfaces found
    if not all_results and raise_if_not_found:
        raise NoStoichiometricSymmetricSurfaceError(
            miller_index=miller_tuple,
            strategies_tried=strategies_tried,
            n_terminations_checked=n_terminations_checked,
            n_stoichiometric=n_stoichiometric,
            n_symmetric=n_symmetric,
        )

    return all_results


def analyze_miller_feasibility(
    structure: t.Union[orm.StructureData, 'Structure'],
    miller_indices: list[t.Union[list, tuple]],
    min_slab_thickness: float = 15.0,
    min_vacuum_thickness: float = 15.0,
    strategies: list[str] = None,
    bonds: dict = None,
    auto_detect_bonds: bool = False,
    max_thickness: float = 30.0,
) -> dict[tuple[int, int, int], MillerFeasibilityReport]:
    """
    Analyze feasibility of finding stoichiometric+symmetric surfaces.

    This diagnostic function checks multiple Miller indices and generates
    a feasibility report for each, helping users decide which orientations
    are suitable for simple surface energy calculations vs. thermodynamics.

    Args:
        structure: Bulk structure
        miller_indices: List of Miller indices to analyze
        min_slab_thickness: Minimum slab thickness
        min_vacuum_thickness: Minimum vacuum thickness
        strategies: Search strategies to use
        bonds: Bond dictionary for bond_preservation
        auto_detect_bonds: Auto-detect bonds
        max_thickness: Maximum thickness for scanning

    Returns:
        Dictionary mapping Miller index tuple to MillerFeasibilityReport

    Example:
        >>> reports = analyze_miller_feasibility(structure, [(1,0,0), (1,1,0), (1,1,1)])
        >>> for miller, report in reports.items():
        ...     print(f"{miller}: {'OK' if report.has_valid_surfaces else 'THERMODYNAMICS'}")
    """
    from pymatgen.core.surface import SlabGenerator
    from pymatgen.io.ase import AseAtomsAdaptor

    # Convert structure if needed
    if isinstance(structure, orm.StructureData):
        ase_atoms = structure.get_ase()
        pmg_structure = AseAtomsAdaptor.get_structure(ase_atoms)
    else:
        pmg_structure = structure

    bulk_composition = get_bulk_composition(pmg_structure)
    reports = {}

    for miller in miller_indices:
        miller_tuple = tuple(int(m) for m in miller)

        # Try to find valid surfaces
        try:
            results = find_stoichiometric_symmetric_slabs(
                pmg_structure,
                miller_tuple,
                min_slab_thickness=min_slab_thickness,
                min_vacuum_thickness=min_vacuum_thickness,
                strategies=strategies,
                bonds=bonds,
                auto_detect_bonds=auto_detect_bonds,
                max_thickness=max_thickness,
                raise_if_not_found=False,  # Don't raise, we want the report
            )
        except Exception as e:
            results = []
            warnings.warn(f"Error analyzing {miller_tuple}: {e}")

        # Gather statistics
        generator = SlabGenerator(
            pmg_structure,
            miller_tuple,
            min_slab_size=min_slab_thickness,
            min_vacuum_size=min_vacuum_thickness,
            in_unit_planes=False,
        )

        all_slabs = generator.get_slabs(symmetrize=False, repair=True)
        n_stoich = 0
        n_sym = 0
        n_both = 0

        for slab in all_slabs:
            is_stoich, is_sym, _ = check_slab_properties(slab, bulk_composition)
            if is_stoich:
                n_stoich += 1
            if is_sym:
                n_sym += 1
            if is_stoich and is_sym:
                n_both += 1

        # Generate notes
        notes = []
        if n_both == 0:
            if n_stoich > 0 and n_sym > 0:
                notes.append("Stoichiometric and symmetric terminations exist but don't overlap")
            elif n_stoich == 0:
                notes.append("No stoichiometric terminations found")
            elif n_sym == 0:
                notes.append("No symmetric terminations found")
            notes.append("Consider using thermodynamics approach")
        else:
            notes.append(f"Found {n_both} valid termination(s)")

        # Determine best thickness and strategy
        best_thickness = None
        best_strategy = ''
        if results:
            best_result = results[0]
            best_thickness = best_result.thickness_angstrom
            best_strategy = best_result.strategy_used

        report = MillerFeasibilityReport(
            miller_index=miller_tuple,
            has_valid_surfaces=len(results) > 0,
            n_terminations_checked=len(all_slabs),
            n_stoichiometric=n_stoich,
            n_symmetric=n_sym,
            n_both=n_both,
            best_thickness=best_thickness,
            recommended_strategy=best_strategy,
            notes=notes,
        )

        reports[miller_tuple] = report

    return reports


@task.calcfunction
def generate_stoichiometric_symmetric_slabs(
    bulk_structure: orm.StructureData,
    miller_indices: orm.List,
    min_slab_thickness: orm.Float = None,
    min_vacuum_thickness: orm.Float = None,
    lll_reduce: orm.Bool = None,
    center_slab: orm.Bool = None,
    symmetrize: orm.Bool = None,
    primitive: orm.Bool = None,
    strategies: orm.List = None,
    max_thickness: orm.Float = None,
) -> dict:
    """
    AiiDA calcfunction to generate stoichiometric+symmetric slabs.

    This is the AiiDA-compatible wrapper around find_stoichiometric_symmetric_slabs.
    It returns a dictionary of StructureData objects for valid terminations.

    Args:
        bulk_structure: Bulk structure as StructureData
        miller_indices: List of Miller indices
        min_slab_thickness: Minimum slab thickness (default: 15.0)
        min_vacuum_thickness: Minimum vacuum thickness (default: 15.0)
        lll_reduce: Use LLL reduction (default: True)
        center_slab: Center slab (default: True)
        symmetrize: This parameter is ignored (we find symmetric slabs automatically)
        primitive: Find primitive cell (default: True)
        strategies: List of strategies (default: all)
        max_thickness: Maximum thickness for scanning (default: 30.0)

    Returns:
        Dictionary with:
        - 'slabs': Namespace of {hkl_key: {term_key: StructureData}}
        - 'metadata': Namespace of {hkl_key: {term_key: Dict with search details}}
        - 'feasibility': Dict with feasibility report for each Miller index
    """
    from pymatgen.io.ase import AseAtomsAdaptor
    import ase

    # Parse parameters
    min_thickness = 15.0 if min_slab_thickness is None else min_slab_thickness.value
    min_vacuum = 15.0 if min_vacuum_thickness is None else min_vacuum_thickness.value
    use_lll = True if lll_reduce is None else lll_reduce.value
    use_center = True if center_slab is None else center_slab.value
    use_primitive = True if primitive is None else primitive.value
    max_thick = 30.0 if max_thickness is None else max_thickness.value
    strategy_list = None if strategies is None else strategies.get_list()

    # Get Miller indices
    miller_list = miller_indices.get_list()

    # Results containers
    slabs_output = {}
    metadata_output = {}
    feasibility_output = {}

    for miller in miller_list:
        miller_tuple = tuple(int(m) for m in miller)
        hkl_key = f"hkl_{miller_tuple[0]}{miller_tuple[1]}{miller_tuple[2]}"

        try:
            results = find_stoichiometric_symmetric_slabs(
                bulk_structure,
                miller_tuple,
                min_slab_thickness=min_thickness,
                min_vacuum_thickness=min_vacuum,
                lll_reduce=use_lll,
                center_slab=use_center,
                primitive=use_primitive,
                strategies=strategy_list,
                max_thickness=max_thick,
                raise_if_not_found=True,
            )

            slabs_output[hkl_key] = {}
            metadata_output[hkl_key] = {}

            for i, result in enumerate(results):
                term_key = f"term_{i}"

                # Convert slab to StructureData
                ase_atoms = AseAtomsAdaptor.get_atoms(result.slab)
                structure_data = orm.StructureData(ase=ase_atoms)

                slabs_output[hkl_key][term_key] = structure_data
                metadata_output[hkl_key][term_key] = orm.Dict(dict=result.to_dict())

            feasibility_output[hkl_key] = {
                'has_valid_surfaces': True,
                'n_valid': len(results),
            }

        except NoStoichiometricSymmetricSurfaceError as e:
            # Store the error information
            feasibility_output[hkl_key] = {
                'has_valid_surfaces': False,
                'error': str(e),
                'n_terminations_checked': e.n_terminations_checked,
                'n_stoichiometric': e.n_stoichiometric,
                'n_symmetric': e.n_symmetric,
            }
            # Re-raise to fail the workflow
            raise

    return {
        'slabs': slabs_output,
        'metadata': metadata_output,
        'feasibility': orm.Dict(dict=feasibility_output),
    }


def get_feasibility_summary(
    reports: t.Union[dict, orm.Dict],
) -> str:
    """
    Generate a human-readable summary of feasibility reports.

    Args:
        reports: Dictionary of MillerFeasibilityReport objects or Dict node

    Returns:
        Formatted summary string
    """
    if isinstance(reports, orm.Dict):
        reports_dict = reports.get_dict()
    else:
        reports_dict = reports

    lines = []
    lines.append("=" * 60)
    lines.append("STOICHIOMETRIC+SYMMETRIC SURFACE FEASIBILITY ANALYSIS")
    lines.append("=" * 60)

    n_ok = 0
    n_fail = 0

    for miller_str, report_data in sorted(reports_dict.items()):
        if isinstance(report_data, MillerFeasibilityReport):
            report = report_data
        else:
            # Reconstruct from dict
            report = MillerFeasibilityReport(**report_data)

        status = "OK" if report.has_valid_surfaces else "NEEDS THERMODYNAMICS"
        if report.has_valid_surfaces:
            n_ok += 1
        else:
            n_fail += 1

        lines.append(f"\n{miller_str}: {status}")
        lines.append(f"  Terminations: {report.n_terminations_checked}")
        lines.append(f"  Stoichiometric: {report.n_stoichiometric}")
        lines.append(f"  Symmetric: {report.n_symmetric}")
        lines.append(f"  Both: {report.n_both}")

        if report.best_thickness:
            lines.append(f"  Best thickness: {report.best_thickness:.1f} Å")
        if report.recommended_strategy:
            lines.append(f"  Strategy: {report.recommended_strategy}")
        if report.notes:
            for note in report.notes:
                lines.append(f"  Note: {note}")

    lines.append("\n" + "=" * 60)
    lines.append(f"SUMMARY: {n_ok} OK, {n_fail} need thermodynamics approach")
    lines.append("=" * 60)

    return '\n'.join(lines)

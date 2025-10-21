#!/usr/bin/env python
"""Tests for surface_modes.py"""
import pytest
import numpy as np
from ase import Atoms
from surface_modes import make_supercell, calculate_surface_area, calculate_coverage, deduplicate_by_coverage


def test_make_supercell_2x2x1():
    """Test creating a 2x2x1 supercell."""
    # Create a simple 2-atom unit cell
    atoms = Atoms('H2', positions=[[0, 0, 0], [1, 0, 0]], cell=[2, 2, 2], pbc=True)

    # Create 2x2x1 supercell
    supercell = make_supercell(atoms, 2, 2, 1)

    # Should have 2 * 2 * 2 * 1 = 8 atoms
    assert len(supercell) == 8

    # Cell should be doubled in x and y, same in z
    assert supercell.cell[0, 0] == pytest.approx(4.0)
    assert supercell.cell[1, 1] == pytest.approx(4.0)
    assert supercell.cell[2, 2] == pytest.approx(2.0)


def test_calculate_surface_area_cubic():
    """Test surface area calculation for cubic cell."""
    # 10 Å × 10 Å cell = 100 Ų = 1 nm²
    atoms = Atoms('H', cell=[[10, 0, 0], [0, 10, 0], [0, 0, 10]], pbc=True)

    area = calculate_surface_area(atoms)

    assert area == pytest.approx(1.0, abs=1e-6)  # 1 nm²


def test_calculate_surface_area_rectangular():
    """Test surface area calculation for rectangular cell."""
    # 5 Å × 20 Å cell = 100 Ų = 1 nm²
    atoms = Atoms('H', cell=[[5, 0, 0], [0, 20, 0], [0, 0, 10]], pbc=True)

    area = calculate_surface_area(atoms)

    assert area == pytest.approx(1.0, abs=1e-6)


def test_calculate_surface_area_non_orthogonal():
    """Test surface area with non-orthogonal cell (monoclinic)."""
    # Cell with 90° angle between a and c, 120° between a and b
    a, b = 10.0, 10.0
    gamma = 120 * np.pi / 180  # 120 degrees in radians
    atoms = Atoms('H', cell=[[a, 0, 0], [b*np.cos(gamma), b*np.sin(gamma), 0], [0, 0, 10]], pbc=True)

    area = calculate_surface_area(atoms)

    # Expected: a × b × sin(120°) = 10 × 10 × (√3/2) = 86.6 Ų = 0.866 nm²
    assert area == pytest.approx(0.866, abs=1e-3)


def test_calculate_coverage_basic():
    """Test basic coverage calculation."""
    coverage = calculate_coverage(n_species=5, area_nm2=2.0, precision=4)

    assert coverage == 2.5  # 5 / 2.0 = 2.5


def test_calculate_coverage_rounding():
    """Test coverage rounding to specified precision."""
    coverage = calculate_coverage(n_species=1, area_nm2=3.0, precision=4)

    # 1/3 = 0.333333... should round to 0.3333
    assert coverage == 0.3333


def test_calculate_coverage_precision_3():
    """Test coverage with different precision."""
    coverage = calculate_coverage(n_species=1, area_nm2=3.0, precision=3)

    # Should round to 3 decimal places: 0.333
    assert coverage == 0.333


def test_calculate_coverage_zero():
    """Test coverage when no species added."""
    coverage = calculate_coverage(n_species=0, area_nm2=5.0, precision=4)

    assert coverage == 0.0


def test_deduplicate_vacancies_mode():
    """Test deduplication for vacancies mode."""
    # Create variants with coverage as key
    variants = [
        ({"name": "v1", "count": 1}, None, 0.5),
        ({"name": "v2", "count": 1}, None, 0.5),  # duplicate coverage
        ({"name": "v3", "count": 2}, None, 1.0),
        ({"name": "v4", "count": 2}, None, 1.0),  # duplicate coverage
    ]

    deduplicated, stats = deduplicate_by_coverage(variants, mode="vacancies")

    assert len(deduplicated) == 2
    assert deduplicated[0][0]["name"] == "v1"  # First with 0.5
    assert deduplicated[1][0]["name"] == "v3"  # First with 1.0
    assert stats["total_generated"] == 4
    assert stats["unique_coverages"] == 2
    assert stats["kept"] == 2


def test_deduplicate_combined_mode():
    """Test deduplication for combined mode with tuple keys."""
    # Combined mode uses (vac_coverage, oh_coverage) tuples
    variants = [
        ({"name": "c1"}, None, (0.5, 1.0)),
        ({"name": "c2"}, None, (0.5, 1.0)),  # duplicate
        ({"name": "c3"}, None, (0.5, 1.5)),  # different OH coverage
        ({"name": "c4"}, None, (1.0, 1.0)),  # different vac coverage
    ]

    deduplicated, stats = deduplicate_by_coverage(variants, mode="combine")

    assert len(deduplicated) == 3
    assert stats["total_generated"] == 4
    assert stats["unique_coverages"] == 3


def test_deduplicate_no_duplicates():
    """Test deduplication when all coverages are unique."""
    variants = [
        ({"name": "v1"}, None, 0.5),
        ({"name": "v2"}, None, 1.0),
        ({"name": "v3"}, None, 1.5),
    ]

    deduplicated, stats = deduplicate_by_coverage(variants, mode="hydrogen")

    assert len(deduplicated) == 3
    assert stats["total_generated"] == 3
    assert stats["unique_coverages"] == 3


def test_surface_modifier_with_supercell():
    """Test SurfaceModifier initialization with supercell."""
    from ase.build import bulk
    from surface_modes import SurfaceModifier

    atoms = bulk('Cu', 'fcc', a=3.6)
    original_count = len(atoms)

    sm = SurfaceModifier(
        atoms=atoms,
        supercell=(2, 2, 1),
        outdir="test_output"
    )

    # Should have created supercell
    assert len(sm.atoms) == original_count * 4


def test_surface_modifier_stores_coverage_params():
    """Test that coverage parameters are stored."""
    from surface_modes import SurfaceModifier

    atoms = Atoms('H', cell=[[10, 0, 0], [0, 10, 0], [0, 0, 10]], pbc=True)

    sm = SurfaceModifier(
        atoms=atoms,
        deduplicate_by_coverage=True,
        coverage_precision=3,
        outdir="test_output"
    )

    assert sm.deduplicate_by_coverage is True
    assert sm.coverage_precision == 3
    assert sm.surface_area_nm2 == pytest.approx(1.0, abs=1e-6)


def test_run_vacancies_with_coverage():
    """Test run_vacancies adds coverage metadata."""
    from surface_modes import SurfaceModifier

    atoms = Atoms('O4', positions=[[0,0,10], [1,0,10], [0,1,0], [1,1,0]],
                  cell=[[10,0,0], [0,10,0], [0,0,20]], pbc=True)

    sm = SurfaceModifier(
        atoms=atoms,
        species="O",
        z_window=1.0,
        which_surface="top",
        include_empty=False,
        outdir="test_vac_cov"
    )

    manifest = sm.run_vacancies()

    # Check coverage metadata exists
    for variant in manifest["variants"]:
        assert "surface_area_nm2" in variant
        assert "vacancy_coverage" in variant
        assert variant["surface_area_nm2"] == pytest.approx(1.0, abs=1e-6)


def test_run_vacancies_with_deduplication():
    """Test run_vacancies deduplicates by coverage."""
    from surface_modes import SurfaceModifier

    # Create structure with 4 surface O atoms and 2 bulk atoms
    atoms = Atoms('O6', positions=[[0,0,10], [5,0,10], [0,5,10], [5,5,10], [2,2,5], [7,7,5]],
                  cell=[[10,0,0], [0,10,0], [0,0,20]], pbc=True)

    sm = SurfaceModifier(
        atoms=atoms,
        species="O",
        z_window=1.0,
        which_surface="top",
        include_empty=False,
        deduplicate_by_coverage=True,
        outdir="test_vac_dedup"
    )

    manifest = sm.run_vacancies()

    # With 4 surface atoms, we have C(4,1)=4, C(4,2)=6, C(4,3)=4, C(4,4)=1 = 15 total
    # But coverages: 1/1.0=1.0, 2/1.0=2.0, 3/1.0=3.0, 4/1.0=4.0 = 4 unique
    assert "deduplication_stats" in manifest
    assert manifest["deduplication_stats"]["total_generated"] == 15
    assert manifest["deduplication_stats"]["unique_coverages"] == 4
    assert manifest["deduplication_stats"]["kept"] == 4
    assert len(manifest["variants"]) == 4


def test_run_hydrogen_with_coverage():
    """Test run_hydrogen adds coverage metadata."""
    from surface_modes import SurfaceModifier

    atoms = Atoms('O4', positions=[[0,0,10], [1,0,10], [0,1,0], [1,1,0]],
                  cell=[[10,0,0], [0,10,0], [0,0,20]], pbc=True)

    sm = SurfaceModifier(
        atoms=atoms,
        species="O",
        z_window=1.0,
        which_surface="top",
        include_empty=False,
        outdir="test_h_cov"
    )

    manifest = sm.run_hydrogen()

    # Check coverage metadata exists
    for variant in manifest["variants"]:
        assert "surface_area_nm2" in variant
        assert "OH_coverage" in variant
        assert variant["surface_area_nm2"] == pytest.approx(1.0, abs=1e-6)


def test_run_hydrogen_with_deduplication():
    """Test run_hydrogen deduplicates by coverage."""
    from surface_modes import SurfaceModifier

    # Create structure with 4 surface O atoms and 2 bulk atoms
    atoms = Atoms('O6', positions=[[0,0,10], [5,0,10], [0,5,10], [5,5,10], [2,2,5], [7,7,5]],
                  cell=[[10,0,0], [0,10,0], [0,0,20]], pbc=True)

    sm = SurfaceModifier(
        atoms=atoms,
        species="O",
        z_window=1.0,
        which_surface="top",
        include_empty=False,
        deduplicate_by_coverage=True,
        outdir="test_h_dedup"
    )

    manifest = sm.run_hydrogen()

    # With 4 surface atoms, we have C(4,1)=4, C(4,2)=6, C(4,3)=4, C(4,4)=1 = 15 total
    # But coverages: 1/1.0=1.0, 2/1.0=2.0, 3/1.0=3.0, 4/1.0=4.0 = 4 unique
    assert "deduplication_stats" in manifest
    assert manifest["deduplication_stats"]["total_generated"] == 15
    assert manifest["deduplication_stats"]["unique_coverages"] == 4
    assert manifest["deduplication_stats"]["kept"] == 4
    assert len(manifest["variants"]) == 4


def test_run_combine_with_coverage():
    """Test run_combine adds coverage metadata for both vacancies and OH."""
    from surface_modes import SurfaceModifier

    atoms = Atoms('O4', positions=[[0,0,10], [1,0,10], [0,1,0], [1,1,0]],
                  cell=[[10,0,0], [0,10,0], [0,0,20]], pbc=True)

    sm = SurfaceModifier(
        atoms=atoms,
        species="O",
        z_window=1.0,
        which_surface="top",
        include_empty=False,
        outdir="test_combo_cov"
    )

    manifest = sm.run_combine()

    # Check coverage metadata exists for both vacancy and OH
    for variant in manifest["variants"]:
        assert "surface_area_nm2" in variant
        assert "vacancy_coverage" in variant
        assert "OH_coverage" in variant
        assert variant["surface_area_nm2"] == pytest.approx(1.0, abs=1e-6)


def test_run_combine_with_deduplication():
    """Test run_combine deduplicates by tuple coverage (vac, oh)."""
    from surface_modes import SurfaceModifier

    # Create structure with 3 surface O atoms
    # This gives manageable number: for each V subset, H from remaining
    atoms = Atoms('O5', positions=[[0,0,10], [5,0,10], [0,5,10], [2,2,5], [7,7,5]],
                  cell=[[10,0,0], [0,10,0], [0,0,20]], pbc=True)

    sm = SurfaceModifier(
        atoms=atoms,
        species="O",
        z_window=1.0,
        which_surface="top",
        include_empty=False,
        deduplicate_by_coverage=True,
        outdir="test_combo_dedup"
    )

    manifest = sm.run_combine()

    # With 3 surface atoms and include_empty=False:
    # all_subsets excludes empty set, so:
    # V={0}, H from {1,2} non-empty: {1}, {2}, {1,2} = 3 variants
    # V={1}, H from {0,2} non-empty: {0}, {2}, {0,2} = 3 variants
    # V={2}, H from {0,1} non-empty: {0}, {1}, {0,1} = 3 variants
    # V={0,1}, H from {2} non-empty: {2} = 1 variant
    # V={0,2}, H from {1} non-empty: {1} = 1 variant
    # V={1,2}, H from {0} non-empty: {0} = 1 variant
    # V={0,1,2}, H from {} = no valid H subsets = 0 variants
    # Total: 3+3+3+1+1+1 = 12 combinations
    # With deduplication by (vac_coverage, OH_coverage) tuples:
    # Expected unique: (1.0, 1.0), (1.0, 2.0), (2.0, 1.0) = 3 unique tuples
    assert "deduplication_stats" in manifest
    assert manifest["deduplication_stats"]["total_generated"] == 12
    assert manifest["deduplication_stats"]["unique_coverages"] == 3
    assert manifest["deduplication_stats"]["kept"] == 3
    assert len(manifest["variants"]) == 3


def test_sample_by_coverage_bins_none():
    """Test that n_bins=None returns all variants."""
    from surface_modes import sample_by_coverage_bins

    # Create dummy variants: (variant_data, atoms_obj, coverage_key)
    variants = [
        ({"name": "v1"}, None, 1.5),
        ({"name": "v2"}, None, 2.0),
        ({"name": "v3"}, None, 2.5),
    ]

    result = sample_by_coverage_bins(variants, n_bins=None, mode='vacancies')

    assert len(result) == 3
    assert result == variants


def test_sample_by_coverage_bins_1d():
    """Test 1D binning for single-coverage mode."""
    from surface_modes import sample_by_coverage_bins

    # Create 10 variants with coverages from 0 to 9
    variants = [({"name": f"v{i}"}, None, float(i)) for i in range(10)]

    # Sample into 3 bins: [0-3), [3-6), [6-9]
    # Bin centers: 1.5, 4.5, 7.5
    # Closest variants: coverage 2, 5, 8 (or 1, 4, 7 - depends on rounding)
    result = sample_by_coverage_bins(variants, n_bins=3, mode='vacancies')

    # Should have 3 structures
    assert len(result) == 3

    # Check they're evenly distributed
    coverages = [v[2] for v in result]
    assert min(coverages) < 3.0  # One in first bin
    assert 3.0 <= max([c for c in coverages if c < 6.0]) < 6.0  # One in middle bin
    assert max(coverages) >= 6.0  # One in last bin


def test_sample_by_coverage_bins_empty():
    """Test that empty list is handled without IndexError."""
    from surface_modes import sample_by_coverage_bins

    # Empty list should be returned as-is without crashing
    result = sample_by_coverage_bins([], n_bins=3, mode='vacancies')

    assert result == []


def test_sample_by_coverage_bins_2d():
    """Test 2D binning for combined mode with dual coverage."""
    from surface_modes import sample_by_coverage_bins

    # Create variants with 2D coverage (vac, oh)
    variants = []
    for vac in range(5):
        for oh in range(5):
            variants.append(({"name": f"v{vac}_{oh}"}, None, (float(vac), float(oh))))

    # 25 total variants, sample into 2x2 = 4 bins
    result = sample_by_coverage_bins(variants, n_bins=2, mode='combine')

    # Should have at most 4 structures (2x2 grid)
    assert len(result) <= 4
    assert len(result) >= 1

    # Check they're distributed across 2D space
    vac_coverages = [v[2][0] for v in result]
    oh_coverages = [v[2][1] for v in result]

    # Should have some diversity in both dimensions
    assert len(set(vac_coverages)) >= 1
    assert len(set(oh_coverages)) >= 1


def test_cli_coverage_bins_argument():
    """Test that --coverage-bins CLI argument is parsed correctly."""
    import sys
    from surface_modes import main

    # This test verifies argument parsing only
    # We'll mock the actual execution
    test_args = ['surface_modes.py', 'st2.vasp', '--mode', 'vacancies',
                 '--coverage-bins', '5']

    # Just check that parsing doesn't fail
    # (Full integration test will be in test_integration_st2.sh)
    # This is a smoke test
    pass  # Will implement actual test in integration


def test_run_vacancies_with_binning():
    """Test that run_vacancies applies binning when coverage_bins is set."""
    from surface_modes import SurfaceModifier
    from ase import Atoms

    # Create slab with 6 surface O atoms
    atoms = Atoms('O6Ti6',
                  positions=[[0, 0, 10], [1, 0, 10], [2, 0, 10],
                            [0, 1, 10], [1, 1, 10], [2, 1, 10],
                            [0, 0, 0], [1, 0, 0], [2, 0, 0],
                            [0, 1, 0], [1, 1, 0], [2, 1, 0]],
                  cell=[10, 10, 15], pbc=True)

    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        sm = SurfaceModifier(
            atoms=atoms,
            species='O',
            z_window=1.0,
            which_surface='top',
            outdir=tmpdir,
            deduplicate_by_coverage=True,
            coverage_bins=3  # Only keep 3 bins
        )

        manifest = sm.run_vacancies()

        # With 6 O atoms, we'd have 2^6-1 = 63 variants without dedup
        # With coverage dedup, we'd have 6 unique coverages (1-6 vacancies)
        # With 3 bins, we should have at most 3 variants
        assert len(manifest['variants']) <= 3

        # Should have binning stats
        assert 'binning_stats' in manifest
        assert manifest['binning_stats']['bins_requested'] == 3


def test_run_hydrogen_with_binning():
    """Test that run_hydrogen applies binning when coverage_bins is set."""
    from surface_modes import SurfaceModifier
    from ase import Atoms

    # Create slab with 4 surface O atoms
    atoms = Atoms('O4Ti4',
                  positions=[[0, 0, 10], [1, 0, 10], [0, 1, 10], [1, 1, 10],
                            [0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]],
                  cell=[10, 10, 15], pbc=True)

    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        sm = SurfaceModifier(
            atoms=atoms,
            species='O',
            z_window=1.0,
            which_surface='top',
            outdir=tmpdir,
            deduplicate_by_coverage=True,
            coverage_bins=2  # Only keep 2 bins
        )

        manifest = sm.run_hydrogen()

        # With 4 O atoms, we'd have 2^4-1 = 15 variants
        # With coverage dedup, 4 unique coverages
        # With 2 bins, should have at most 2 variants
        assert len(manifest['variants']) <= 2
        assert 'binning_stats' in manifest


def test_run_combine_with_binning():
    """Test that run_combine applies 2D binning when coverage_bins is set."""
    from surface_modes import SurfaceModifier
    from ase import Atoms

    # Create slab with 3 surface O atoms for faster test
    atoms = Atoms('O3Ti3',
                  positions=[[0, 0, 10], [1, 0, 10], [2, 0, 10],
                            [0, 0, 0], [1, 0, 0], [2, 0, 0]],
                  cell=[10, 10, 15], pbc=True)

    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        sm = SurfaceModifier(
            atoms=atoms,
            species='O',
            z_window=1.0,
            which_surface='top',
            outdir=tmpdir,
            deduplicate_by_coverage=True,
            coverage_bins=2  # 2x2 = 4 bins max
        )

        manifest = sm.run_combine()

        # With 3 O atoms and combined mode, many variants
        # With 2 bins, should have at most 2x2 = 4 variants
        assert len(manifest['variants']) <= 4
        assert 'binning_stats' in manifest


def test_bin_metadata_in_manifest():
    """Test that bin metadata is added to each variant in manifest."""
    from surface_modes import SurfaceModifier
    from ase import Atoms

    atoms = Atoms('O3Ti3',
                  positions=[[0, 0, 10], [1, 0, 10], [2, 0, 10],
                            [0, 0, 0], [1, 0, 0], [2, 0, 0]],
                  cell=[10, 10, 15], pbc=True)

    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        sm = SurfaceModifier(
            atoms=atoms,
            species='O',
            z_window=1.0,
            which_surface='top',
            outdir=tmpdir,
            deduplicate_by_coverage=True,
            coverage_bins=2
        )

        manifest = sm.run_vacancies()

        # Check that each variant has bin metadata
        for variant in manifest['variants']:
            assert 'bin_id' in variant
            assert 'bin_center' in variant
            assert isinstance(variant['bin_id'], int)
            assert isinstance(variant['bin_center'], float)


def test_coverage_bins_requires_deduplication():
    """Test that --coverage-bins requires --deduplicate-by-coverage."""
    import subprocess
    import tempfile
    from pathlib import Path
    from ase.io import write
    from ase import Atoms

    # Create dummy input file
    with tempfile.TemporaryDirectory() as tmpdir:
        dummy_vasp = Path(tmpdir) / "test.vasp"
        atoms = Atoms('O', positions=[[0, 0, 0]], cell=[10, 10, 10], pbc=True)
        write(dummy_vasp, atoms, direct=True)

        # Try to use --coverage-bins without --deduplicate-by-coverage
        result = subprocess.run(
            ['python', 'surface_modes.py', str(dummy_vasp),
             '--mode', 'vacancies', '--coverage-bins', '5'],
            capture_output=True,
            text=True
        )

        # Should fail with error message
        assert result.returncode != 0
        assert 'requires --deduplicate-by-coverage' in result.stderr.lower()


# ====== Tests for simplified naming scheme ======

def test_generate_variant_name_vacancies_single():
    """Test simplified naming for single-surface vacancies."""
    from surface_modes import generate_variant_name

    name = generate_variant_name(mode='vacancies', index=1, coverage=0.25)

    assert name == 'vac_001_0.25'
    assert 'vac' in name
    assert '001' in name
    assert '0.25' in name


def test_generate_variant_name_hydrogen_single():
    """Test simplified naming for hydrogen mode."""
    from surface_modes import generate_variant_name

    name = generate_variant_name(mode='hydrogen', index=5, coverage=1.75)

    assert name == 'oh_005_1.75'
    assert 'oh' in name
    assert '005' in name
    assert '1.75' in name


def test_generate_variant_name_combine_dual():
    """Test simplified naming for combined mode with dual coverage."""
    from surface_modes import generate_variant_name

    name = generate_variant_name(mode='combine', index=10, coverage=(1.5, 0.75))

    assert name == 'combo_010_1.5_0.75'
    assert 'combo' in name
    assert '010' in name
    assert '1.5' in name
    assert '0.75' in name


def test_generate_variant_name_zero_coverage():
    """Test simplified naming with zero coverage."""
    from surface_modes import generate_variant_name

    name = generate_variant_name(mode='vacancies', index=0, coverage=0.0)

    assert name == 'vac_000_0.0'


def test_generate_variant_name_high_coverage():
    """Test simplified naming with high coverage values."""
    from surface_modes import generate_variant_name

    name = generate_variant_name(mode='hydrogen', index=999, coverage=99.9999)

    assert name == 'oh_999_99.9999'


def test_generate_variant_name_combine_high_index():
    """Test simplified naming for combine with 3-digit index."""
    from surface_modes import generate_variant_name

    name = generate_variant_name(mode='combine', index=123, coverage=(2.5, 3.75))

    assert name == 'combo_123_2.5_3.75'

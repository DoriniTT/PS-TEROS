"""
Charge-Neutral Termination Tests

Tier 1 tests use pymatgen only and check the termination finder on small,
well-understood cases: Ag3PO4(110) (the electron-count problem that motivated
the module), rutile SnO2, rock-salt MgO and wurtzite ZnO. One tier 2 test runs
``generate_slab_structures`` in ``charge_neutral`` mode inside a WorkGraph.
"""

import os

import numpy as np
import pytest

pymatgen = pytest.importorskip('pymatgen')

from pymatgen.core import Lattice, Structure  # noqa: E402
from pymatgen.core.surface import SlabGenerator  # noqa: E402

from psteros.core.terminations import (  # noqa: E402
    NoChargeNeutralTerminationError,
    classify_slab,
    find_charge_neutral_terminations,
    find_face_operation,
    formal_charge,
    has_face_reversing_operation,
    normal_repeat,
    resolve_oxidation_states,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AG3PO4_CIF = os.path.join(ROOT, 'examples', 'vasp', 'structures', 'ag3po4.cif')
AG3PO4_STATES = {'Ag': 1, 'P': 5, 'O': -2}
PO4_BONDS = {('P', 'O'): 1.9}


@pytest.fixture(scope='module')
def ag3po4():
    return Structure.from_file(AG3PO4_CIF)


@pytest.fixture(scope='module')
def ag3po4_110(ag3po4):
    return find_charge_neutral_terminations(
        ag3po4, (1, 1, 0), 10.0, 15.0,
        oxidation_states=AG3PO4_STATES, unit_bonds=PO4_BONDS,
    )


@pytest.fixture(scope='module')
def sno2():
    return Structure.from_file(os.path.join(ROOT, 'tests', 'fixtures', 'structures', 'sno2_rutile.vasp'))


@pytest.fixture(scope='module')
def mgo():
    return Structure.from_spacegroup('Fm-3m', Lattice.cubic(4.21), ['Mg', 'O'],
                                     [[0, 0, 0], [0.5, 0.5, 0.5]])


@pytest.fixture(scope='module')
def zno():
    return Structure.from_spacegroup('P6_3mc', Lattice.hexagonal(3.25, 5.21), ['Zn', 'O'],
                                     [[1 / 3, 2 / 3, 0], [1 / 3, 2 / 3, 0.382]])


def _intact_po4(structure):
    """Every P has exactly four O within 1.9 A and every O belongs to one P."""
    phosphorus = [index for index, site in enumerate(structure) if site.specie.symbol == 'P']
    owners = []
    for index in phosphorus:
        oxygens = [n.index for n in structure.get_neighbors(structure[index], 1.9)
                   if n.specie.symbol == 'O']
        if len(oxygens) != 4:
            return False
        owners.extend(oxygens)
    n_oxygen = sum(site.specie.symbol == 'O' for site in structure)
    return len(owners) == len(set(owners)) == n_oxygen


@pytest.mark.tier1
class TestOxidationStates:

    def test_explicit_states_are_checked_for_neutrality(self, ag3po4):
        assert resolve_oxidation_states(ag3po4, AG3PO4_STATES) == {'Ag': 1.0, 'O': -2.0, 'P': 5.0}
        with pytest.raises(ValueError, match='not neutral'):
            resolve_oxidation_states(ag3po4, {'Ag': 2, 'P': 5, 'O': -2})

    def test_missing_element_is_rejected(self, ag3po4):
        with pytest.raises(ValueError, match='missing'):
            resolve_oxidation_states(ag3po4, {'Ag': 1, 'O': -2})

    def test_states_are_guessed_when_omitted(self, mgo):
        assert resolve_oxidation_states(mgo) == {'Mg': 2.0, 'O': -2.0}


@pytest.mark.tier1
class TestGeometry:

    def test_normal_repeat_includes_centring(self, ag3po4, mgo):
        assert normal_repeat(mgo, (1, 0, 0)) == pytest.approx(4.21 / 2)
        assert normal_repeat(mgo, (1, 1, 1)) == pytest.approx(4.21 / np.sqrt(3))
        assert normal_repeat(ag3po4, (1, 1, 0)) == pytest.approx(ag3po4.lattice.a / np.sqrt(2))

    def test_polar_directions(self, ag3po4, zno, mgo):
        assert not has_face_reversing_operation(ag3po4, (1, 1, 1))
        assert not has_face_reversing_operation(zno, (0, 0, 1))
        assert has_face_reversing_operation(ag3po4, (1, 1, 0))
        assert has_face_reversing_operation(mgo, (1, 1, 1))

    def test_polar_direction_raises_with_explanation(self, zno):
        with pytest.raises(NoChargeNeutralTerminationError, match='polar direction'):
            find_charge_neutral_terminations(zno, (0, 0, 1), 10.0)

    def test_extended_unit_bonds_are_rejected(self, ag3po4):
        with pytest.raises(ValueError, match='extended network'):
            find_charge_neutral_terminations(
                ag3po4, (1, 1, 0), 10.0, oxidation_states=AG3PO4_STATES,
                unit_bonds={('Ag', 'O'): 2.6, ('P', 'O'): 1.9},
            )


@pytest.mark.tier1
class TestAg3PO4:

    def test_pymatgen_symmetric_slabs_are_all_charged(self, ag3po4):
        """The motivating problem: no symmetrized SlabGenerator slab is neutral."""
        generator = SlabGenerator(ag3po4, (1, 1, 0), 10.0, 15.0, center_slab=True,
                                  lll_reduce=True, in_unit_planes=False)
        charges = [formal_charge(slab, AG3PO4_STATES) for slab in generator.get_slabs(symmetrize=True)]
        assert charges and all(abs(charge) > 0.5 for charge in charges)

    def test_110_reproduces_ag18p6o24(self, ag3po4_110):
        """Ag-rich Ag20P6O24 (+2) minus one Ag per face is the closed-shell slab."""
        by_formula = {termination.formula: termination for termination in ag3po4_110}
        assert set(by_formula) == {'Ag18P6O24', 'Ag12P4O16'}
        silver = by_formula['Ag18P6O24']
        assert silver.parent_formula == 'Ag20P6O24'
        assert silver.parent_formal_charge == pytest.approx(2.0)
        assert silver.removed_per_face == ('Ag1',)
        phosphate = by_formula['Ag12P4O16']
        assert phosphate.removed_per_face == ('P1O4',)

    def test_110_slabs_are_neutral_symmetric_and_keep_po4(self, ag3po4, ag3po4_110):
        for termination in ag3po4_110:
            report = classify_slab(termination.structure, ag3po4, AG3PO4_STATES)
            assert report['is_charge_neutral'] and report['is_symmetric']
            assert find_face_operation(termination.structure) is not None
            assert _intact_po4(termination.structure)
            assert 10.0 <= termination.thickness < 10.0 + normal_repeat(ag3po4, (1, 1, 0))
            normal = termination.structure.lattice.matrix[2]
            assert np.allclose(termination.structure.lattice.matrix[:2] @ normal, 0, atol=1e-6)

    def test_100_needs_a_supercell(self, ag3po4):
        """Window charges are +-1, so a 1x1 cell cannot be neutralised."""
        with pytest.raises(NoChargeNeutralTerminationError, match=r'-1, \+1'):
            find_charge_neutral_terminations(ag3po4, (1, 0, 0), 10.0, supercell=(1, 1),
                                             oxidation_states=AG3PO4_STATES, unit_bonds=PO4_BONDS)
        terminations = find_charge_neutral_terminations(
            ag3po4, (1, 0, 0), 10.0, oxidation_states=AG3PO4_STATES, unit_bonds=PO4_BONDS,
        )
        assert terminations.supercell == (2, 1)
        _check_all(terminations, ag3po4, AG3PO4_STATES, (1, 0, 0), 10.0)
        for termination in terminations:
            assert _intact_po4(termination.structure)


@pytest.mark.tier1
class TestSimpleOxides:

    def test_sno2_110_is_the_stoichiometric_o_terminated_slab(self, sno2):
        terminations = find_charge_neutral_terminations(sno2, (1, 1, 0), 12.0)
        neutral = [t for t in terminations if not t.removed_per_face]
        assert len(neutral) == 1
        assert neutral[0].is_stoichiometric
        heights = neutral[0].structure.cart_coords[:, 2]
        top = neutral[0].structure[int(np.argmax(heights))]
        assert top.specie.symbol == 'O'

    def test_mgo_100_appears_once(self, mgo):
        terminations = find_charge_neutral_terminations(mgo, (1, 0, 0), 8.0)
        assert [t.formula for t in terminations] == ['Mg5O5']

    def test_mgo_111_half_occupied_outer_layer(self, mgo):
        """The polar rock-salt (111) surface is neutralised by half a layer."""
        terminations = find_charge_neutral_terminations(mgo, (1, 1, 1), 8.0)
        assert terminations.supercell == (2, 1)
        assert sorted(t.removed_per_face for t in terminations) == [('Mg1',), ('O1',)]
        for termination in terminations:
            assert termination.parent_formal_charge != 0
            assert abs(termination.formal_charge) < 1e-6

    def test_cubic_111_uses_the_primitive_surface_cell(self, mgo):
        """SlabGenerator's (111) cell for a cubic bulk is 4x too large."""
        terminations = find_charge_neutral_terminations(mgo, (1, 1, 1), 8.0, supercell=(2, 1))
        primitive_area = np.sqrt(3) / 2 * (4.21 / np.sqrt(2)) ** 2
        assert terminations[0].area == pytest.approx(2 * primitive_area, rel=1e-6)


# =============================================================================
# SURVEY ACROSS MATERIAL CLASSES
# =============================================================================

def _sg(group, lattice, species, coords):
    return Structure.from_spacegroup(group, lattice, species, coords)


MATERIALS = {
    'Cu': lambda: _sg('Fm-3m', Lattice.cubic(3.615), ['Cu'], [[0, 0, 0]]),
    'Mg': lambda: _sg('P6_3/mmc', Lattice.hexagonal(3.209, 5.211), ['Mg'], [[1 / 3, 2 / 3, 1 / 4]]),
    'Cu3Au': lambda: _sg('Pm-3m', Lattice.cubic(3.75), ['Au', 'Cu'], [[0, 0, 0], [0, .5, .5]]),
    'NiAl': lambda: _sg('Pm-3m', Lattice.cubic(2.887), ['Ni', 'Al'], [[0, 0, 0], [.5, .5, .5]]),
    'Si': lambda: _sg('Fd-3m', Lattice.cubic(5.431), ['Si'], [[0, 0, 0]]),
    'GaAs': lambda: _sg('F-43m', Lattice.cubic(5.653), ['Ga', 'As'], [[0, 0, 0], [.25, .25, .25]]),
    'ZnO': lambda: _sg('P6_3mc', Lattice.hexagonal(3.25, 5.21), ['Zn', 'O'],
                       [[1 / 3, 2 / 3, 0], [1 / 3, 2 / 3, .382]]),
    'CeO2': lambda: _sg('Fm-3m', Lattice.cubic(5.411), ['Ce', 'O'], [[0, 0, 0], [.25, .25, .25]]),
    'SrTiO3': lambda: _sg('Pm-3m', Lattice.cubic(3.905), ['Sr', 'Ti', 'O'],
                          [[0, 0, 0], [.5, .5, .5], [.5, .5, 0]]),
    'Al2O3': lambda: _sg('R-3c', Lattice.hexagonal(4.759, 12.991), ['Al', 'O'],
                         [[0, 0, .3523], [.3064, 0, .25]]),
    'CaCO3': lambda: _sg('R-3c', Lattice.hexagonal(4.99, 17.06), ['Ca', 'C', 'O'],
                         [[0, 0, 0], [0, 0, .25], [.257, 0, .25]]),
}
CALCITE = {'oxidation_states': {'Ca': 2, 'C': 4, 'O': -2}, 'unit_bonds': {('C', 'O'): 1.5}}

# material, (hkl), options, expected number of terminations, expected surface cell,
# whether every slab is stoichiometric (None: not checked), physical meaning
SURVEY = [
    ('Cu', (1, 1, 1), {}, 1, (1, 1), True, 'fcc metal, close-packed'),
    ('Mg', (1, 0, 0), {}, 2, (1, 1), True, 'hcp (10-10): short and long interlayer cut'),
    ('Cu3Au', (1, 0, 0), {}, 2, (1, 1), False, 'L1_2 alloy: pure Cu or mixed CuAu planes'),
    ('NiAl', (1, 1, 0), {}, 1, (1, 1), True, 'B2 alloy, mixed planes'),
    ('Si', (1, 1, 1), {}, 2, (1, 1), True, 'diamond (111): shuffle and glide cuts'),
    ('GaAs', (1, 1, 0), {}, 1, (1, 1), True, 'zinc blende cleavage plane'),
    ('GaAs', (1, 0, 0), {}, 2, (2, 1), True, 'Ga or As half layer'),
    ('ZnO', (1, 0, 0), {}, 2, (1, 1), True, 'wurtzite (10-10): one or two bonds cut'),
    ('CeO2', (1, 1, 1), {}, 1, (1, 1), True, 'O-terminated O-Ce-O trilayer'),
    ('CeO2', (1, 0, 0), {}, 1, (1, 1), True, 'half the outer O removed'),
    ('SrTiO3', (1, 0, 0), {}, 2, (1, 1), False, 'SrO- and TiO2-terminated'),
    ('Al2O3', (0, 0, 1), {}, 1, (1, 1), True, 'single-Al termination'),
    ('CaCO3', (1, 0, 4), CALCITE, 1, (1, 1), True, 'calcite cleavage plane, CO3 intact'),
]


def _check_all(terminations, bulk, states, miller, thickness):
    """Invariants every returned slab must satisfy."""
    assert terminations, 'no terminations returned'
    labels = [termination.label for termination in terminations]
    assert labels == [f"term_{index}" for index in range(len(terminations))]
    for termination in terminations:
        slab = termination.structure
        assert abs(formal_charge(slab, states)) < 1e-6
        assert find_face_operation(slab) is not None
        assert thickness <= termination.thickness < thickness + normal_repeat(bulk, miller) + 1e-6
        matrix = slab.lattice.matrix
        assert np.allclose(matrix[:2] @ matrix[2], 0, atol=1e-6)
        distances = slab.distance_matrix + np.eye(len(slab)) * 10
        assert distances.min() > 0.9, 'overlapping atoms'


@pytest.mark.tier1
@pytest.mark.parametrize('name, miller, options, count, cell, stoichiometric, meaning', SURVEY,
                         ids=[f"{row[0]}{''.join(map(str, row[1]))}" for row in SURVEY])
def test_survey(name, miller, options, count, cell, stoichiometric, meaning):
    bulk = MATERIALS[name]()
    terminations = find_charge_neutral_terminations(bulk, miller, 10.0, **options)
    _check_all(terminations, bulk, terminations.oxidation_states, miller, 10.0)
    assert len(terminations) == count, terminations.summary()
    assert terminations.supercell == cell
    if stoichiometric is True:
        assert all(t.is_stoichiometric for t in terminations)
    elif stoichiometric is False:
        assert not all(t.is_stoichiometric for t in terminations)


@pytest.mark.tier1
@pytest.mark.parametrize('name, miller', [('GaAs', (1, 1, 1)), ('ZnO', (0, 0, 1))])
def test_survey_polar(name, miller):
    with pytest.raises(NoChargeNeutralTerminationError, match='polar direction'):
        find_charge_neutral_terminations(MATERIALS[name](), miller, 10.0)


@pytest.mark.tier1
def test_intermetallics_are_treated_as_metallic():
    terminations = find_charge_neutral_terminations(MATERIALS['Cu3Au'](), (1, 0, 0), 10.0)
    assert terminations.is_metallic
    assert 'only face symmetry' in terminations.summary()


@pytest.mark.tier1
def test_calcite_keeps_carbonate_whole():
    terminations = find_charge_neutral_terminations(MATERIALS['CaCO3'](), (1, 0, 4), 10.0, **CALCITE)
    slab = terminations[0].structure
    for site in slab:
        if site.specie.symbol == 'C':
            assert sum(n.specie.symbol == 'O' for n in slab.get_neighbors(site, 1.5)) == 3


# =============================================================================
# USER INTERFACE
# =============================================================================

@pytest.mark.tier1
class TestInterface:

    def test_summary_table(self, ag3po4_110):
        text = ag3po4_110.summary()
        assert text == str(ag3po4_110)
        assert text.splitlines()[0] == 'Ag3PO4(110) | 1x1 surface cell | Ag+1 P+5 O-2 | kept whole: P-O'
        assert 'term_1  Ag18P6O24' in text
        assert 'Ag20P6O24 (+2) minus Ag per face' in text
        assert '<table>' in ag3po4_110._repr_html_()

    def test_origin_text(self, ag3po4_110):
        assert ag3po4_110[0].origin == 'Ag12P6O24 (-6) minus PO4 per face'

    def test_write(self, ag3po4_110, tmp_path):
        import json

        paths = ag3po4_110.write(str(tmp_path / 'slabs'))
        assert sorted(os.path.basename(path) for path in paths) == [
            'term_0_Ag12P4O16.vasp', 'term_1_Ag18P6O24.vasp', 'terminations.json',
        ]
        back = Structure.from_file(str(tmp_path / 'slabs' / 'term_1_Ag18P6O24.vasp'))
        assert back.composition == ag3po4_110[1].structure.composition
        summary = json.loads((tmp_path / 'slabs' / 'terminations.json').read_text())
        assert summary['terminations'][1]['parent_formula'] == 'Ag20P6O24'

    def test_plot(self, ag3po4_110, tmp_path):
        figure = ag3po4_110.plot(str(tmp_path / 'slabs.png'))
        assert (tmp_path / 'slabs.png').stat().st_size > 10000
        assert len([axis for axis in figure.axes if axis.get_visible()]) == 2

    def test_removed_sites_are_recorded(self, ag3po4_110):
        silver = ag3po4_110[1]
        assert [symbol for symbol, _ in silver.removed_sites] == ['Ag', 'Ag']

    def test_command_line(self, tmp_path, capsys):
        from psteros.core.terminations import main

        code = main([AG3PO4_CIF, '110', '--oxidation', 'Ag=1,P=5,O=-2', '--keep', 'P-O:1.9',
                     '--write', str(tmp_path / 'out'), '--plot', str(tmp_path / 'out.png')])
        output = capsys.readouterr().out
        assert code == 0
        assert 'Ag18P6O24' in output and 'Wrote 2 slabs' in output
        assert (tmp_path / 'out.png').exists()

    def test_command_line_reports_polar_surface(self, capsys):
        from psteros.core.terminations import main

        assert main([AG3PO4_CIF, '1', '1', '1', '--oxidation', 'Ag=1,P=5,O=-2']) == 1
        assert 'polar direction' in capsys.readouterr().err


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_generate_slab_structures_charge_neutral_mode():
    """The calcfunction returns the same slabs as the pure-Python finder."""
    from aiida import orm
    from aiida_workgraph import WorkGraph
    from psteros.core.slabs import generate_slab_structures

    bulk = orm.StructureData(pymatgen=Structure.from_file(AG3PO4_CIF))
    wg = WorkGraph('charge_neutral_slabs')
    wg.add_task(
        generate_slab_structures,
        name='generate_slab_structures',
        bulk_structure=bulk,
        miller_indices=orm.List(list=[1, 1, 0]),
        min_slab_thickness=orm.Float(10.0),
        min_vacuum_thickness=orm.Float(15.0),
        lll_reduce=orm.Bool(True),
        center_slab=orm.Bool(True),
        symmetrize=orm.Bool(True),
        primitive=orm.Bool(False),
        termination_mode=orm.Str('charge_neutral'),
        oxidation_states=orm.Dict(dict=AG3PO4_STATES),
        unit_bonds=orm.List(list=[['P', 'O', 1.9]]),
    )
    wg.run()
    task_ = wg.tasks['generate_slab_structures']
    slabs = {label: getattr(task_.outputs.slabs, label).value for label in ('term_0', 'term_1')}
    assert sorted(slab.get_formula() for slab in slabs.values()) == ['Ag12O16P4', 'Ag18O24P6']
    report = task_.outputs.termination_report.value.get_dict()
    assert report['terminations']['term_1']['parent_formula'] == 'Ag20P6O24'
    assert 'Ag18P6O24' in report['summary']


@pytest.mark.tier2
@pytest.mark.requires_aiida
def test_generate_slab_structures_default_mode_is_unchanged():
    """Without termination_mode, the SlabGenerator path returns no report."""
    from aiida import orm
    from aiida_workgraph import WorkGraph
    from psteros.core.slabs import generate_slab_structures

    wg = WorkGraph('pymatgen_slabs')
    wg.add_task(
        generate_slab_structures,
        name='generate_slab_structures',
        bulk_structure=orm.StructureData(pymatgen=Structure.from_file(AG3PO4_CIF)),
        miller_indices=orm.List(list=[1, 1, 0]),
        min_slab_thickness=orm.Float(10.0),
        min_vacuum_thickness=orm.Float(15.0),
        lll_reduce=orm.Bool(True),
        center_slab=orm.Bool(True),
        symmetrize=orm.Bool(True),
        primitive=orm.Bool(False),
    )
    wg.run()
    node = wg.tasks['generate_slab_structures'].process
    assert node.is_finished_ok
    assert 'termination_report' not in node.outputs
    assert len(node.outputs.slabs) == 6

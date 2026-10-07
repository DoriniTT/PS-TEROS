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
        with pytest.raises(NoChargeNeutralTerminationError, match=r'\[-1\.0, 1\.0\]'):
            find_charge_neutral_terminations(ag3po4, (1, 0, 0), 10.0,
                                             oxidation_states=AG3PO4_STATES, unit_bonds=PO4_BONDS)
        terminations = find_charge_neutral_terminations(
            ag3po4, (1, 0, 0), 10.0, oxidation_states=AG3PO4_STATES,
            unit_bonds=PO4_BONDS, supercell=(1, 2),
        )
        assert terminations
        for termination in terminations:
            assert abs(formal_charge(termination.structure, AG3PO4_STATES)) < 1e-6
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
        removed = sorted(set(t.removed_per_face[0] for t in terminations))
        assert removed == ['Mg1', 'O1']
        for termination in terminations:
            assert termination.parent_formal_charge != 0
            assert abs(termination.formal_charge) < 1e-6


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
    outputs = wg.tasks['generate_slab_structures'].outputs.slabs
    slabs = {label: getattr(outputs, label).value for label in ('term_0', 'term_1')}
    assert sorted(slab.get_formula() for slab in slabs.values()) == ['Ag12O16P4', 'Ag18O24P6']

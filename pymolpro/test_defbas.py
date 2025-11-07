from .defbas import mixed_core_correlation_assert, element_ranges, elements_from_xyz
from .defbas import periodic_table
import pytest


def test_mixed_core_correlation_only_valence():
    assert periodic_table.index('Zn') + 1 == 30

    for range in [
        'ga',
        'GA',
        'Ga',
        31,
        'Ga-Kr',
        999,
        -1,
        'Zn-Ga',
        'Zn-Kr',
        'Ga-Xe',
    ]:
        assert mixed_core_correlation_assert(range, False)

    for range in [
        'Zn',
        30,
        'Sc-Zn',
        'Na-Mg',
        'Zn-Ga',
        'Zn-Kr',
        'Ga-Xe',
    ]:
        assert mixed_core_correlation_assert(range)

    for range in [
        'Kr-Ga',
    ]:
        assert not mixed_core_correlation_assert(range)
        assert not mixed_core_correlation_assert(range, False)

    for range in [
        'bad',
        'bad-worse',
        4.0,
        {'a': 'b'},
        True,
    ]:
        with pytest.raises(ValueError):
            mixed_core_correlation_assert(range)


def test_element_ranges():
    for case in [
        (['Cu'], ['Cu']),
        (['Cu', 'Zn'], ['Cu-Zn']),
        (['Cu', 'Ga', 'Zn'], ['Cu-Zn', 'Ga']),
        (['Ga', 'Ge', 'As'], ['Ga-As']),
        (['Ge', 'Ga', 'As'], ['Ga-As']),
        (['Ga', 'Ge', 'Ga', 'As'], ['Ga-As']),
        (['Ga', 'Zn', 'As'], ['Zn', 'Ga', 'As']),
    ]:
        assert element_ranges(case[0]) == case[1]


def test_element_ranges_from_xyz():
    for xyz, ranges in {
        '2\n\nHe 0 0 0\nHe 2 0 0\n': ['He'],
        'He 0 0 0\n Ne 0 0 3': ['He', 'Ne'],
        'He 0 0 0\n H 0 0 3': ['H-He'],
        'C 0 0 0\n O 0 0 -1\nFe 0 0 2': ['C','O','Fe'],
    }.items():
        assert (element_ranges(elements_from_xyz(xyz)) == ranges)

def test_element_ranges_from_xyz_heavy():
    for xyz, ranges in {
        '2\n\nHe 0 0 0\nHe 2 0 0\n': ['He'],
        '3\n\nFe 0 0 0\nCo 2 0 0\nO 4 0 0\n': ['Fe-Co'],
        '3\n\nW 0 0 0\nGa 2 0 0\nO 4 0 0\n': ['W'],
    }.items():
        assert (element_ranges(elements_from_xyz(xyz),light=False) == ranges)
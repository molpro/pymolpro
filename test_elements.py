from pymolpro.elements import ElementRange, element, elements_include_heavy
from pymolpro.elements import mixed_core_correlation_assert, element_ranges, elements_from_xyz
from pymolpro.elements import periodic_table
import pytest


@pytest.mark.parametrize("given,expected", [
    (1, 'H'),
    ("He", 'He'),
    ("HE", 'He'),
    ("he", 'He'),
])
def test_element_symbol(given, expected):
    assert element(given).symbol == expected


def test_element_valid_init():
    # print(len(periodic_table))
    pytest.raises(ValueError, element, -1)
    pytest.raises(ValueError, element, 127)
    pytest.raises(ValueError, element, 'garbage')


def test_mixed_core_correlation_only_valence():
    assert periodic_table.index('Zn') + 1 == 30

    for range in [
        'al',
        'AL',
        'Al',
        13,
        'Al-Ar',
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
        (['Cu', 'Ga', 'Zn'], ['Cu-Ga']),
        (['Ga', 'Ge', 'As'], ['Ga-As']),
        (['Ge', 'Ga', 'As'], ['Ga-As']),
        (['Ga', 'Ge', 'Ga', 'As'], ['Ga-As']),
    ]:
        # print([str(r) for r in element_ranges(case[0])], case[1])
        assert element_ranges(case[0]) == [ElementRange(expected) for expected in case[1]]


def test_elements_from_xyz():
    for xyz, elements in {
        '2\n\nHe 0 0 0\nHe 2 0 0\n': ['He'],
        'He 0 0 0\n Ne 0 0 3': ['He', 'Ne'],
        'He 0 0 0\n H 0 0 3': ['H', 'He'],
        'C 0 0 0\n O 0 0 -1\nFe1 0 0 2': ['C', 'O', 'Fe'],
    }.items():
        # print (elements_from_xyz(xyz) , [element(instance) for instance in elements])
        assert (elements_from_xyz(xyz) == [element(instance) for instance in elements])


def test_element_range_str():
    assert str(ElementRange('H-He')) == 'H-He'
    assert str(ElementRange(2, 4)) == 'He-Be'


def test_element_ranges_from_xyz():
    for xyz, ranges in {
        '2\n\nHe 0 0 0\nHe 2 0 0\n': ['He'],
        'He 0 0 0\n Ne 0 0 3': ['He', 'Ne'],
        'He 0 0 0\n H 0 0 3': ['H-He'],
        'C 0 0 0\n O 0 0 -1\nFe 0 0 2': ['C', 'O', 'Fe'],
    }.items():
        # print ([str(range) for range in element_ranges(elements_from_xyz(xyz))] , [str(ElementRange(range)) for range in ranges])
        assert (element_ranges(elements_from_xyz(xyz)) == [ElementRange(range) for range in ranges])


def test_element_ranges_from_xyz_heavy():
    for xyz, ranges in {
        '2\n\nHe 0 0 0\nHe 2 0 0\n': ['He'],
        '3\n\nFe 0 0 0\nCo 2 0 0\nO 4 0 0\n': ['Fe-Co'],
        '3\n\nW 0 0 0\nGa 2 0 0\nO 4 0 0\n': ['Ga', 'W'],
    }.items():
        # print (element_ranges(elements_from_xyz(xyz),light=False) , [str(ElementRange(range)) for range in ranges])
        assert (element_ranges(elements_from_xyz(xyz), light=False) == [ElementRange(range) for range in ranges])


@pytest.mark.parametrize("given,expected", [
    (['H', 'He', 'B', 'C', 'N', 'O', 'Ne'], False),
    (['H', 'He', 'B', 'C', 'N', 'O', 'Ne', 'Na'], True),
    (['Ar'], False),
    (['Mg'], True),
    ([126], True),
    (['Og'], True),
])
def test_elements_include_heavy(given, expected):
    assert elements_include_heavy(given) == expected

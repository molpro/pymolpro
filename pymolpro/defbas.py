import pathlib
import re


class Defbas:
    def __init__(self, molpro_root=None):
        self.contents = open(pathlib.Path(molpro_root) / 'lib' / 'defbas', 'r').readlines()

    def search(self, element=None, key=None, type=None, context=None):
        r"""
        Search for matching entries

        :param element: Restrict to a particular element
        :type element: str
        :param key: Restrict to a particular key
        :type key: str
        :param type: Type - can be 'aux' or 'poisson' or undefined, meaning orbital
        :type type: str
        :param context: Restrict to a particular context
        :type context: str
        :return:
        :rtype list[dict]:
        """
        result = []
        contexts = context.split(' ') if context else None
        fully_wild = not (element or key or type or context)
        n_extra = 0
        for line in self.contents:
            line = line.strip(' \n')
            if line and line[0] == '!':
                if fully_wild:
                    result.append({})
                    result[-1]['name'] = '!'
                    result[-1]['comment'] = line
                    result[-1]['minz'] = 0
                    result[-1]['maxz'] = 0
                    result[-1]['minang'] = 0
                    result[-1]['maxang'] = 0
                continue
            split_line = line.split(':')
            if n_extra > 0:
                result[-1]['extra'].append(line.strip())
                n_extra -= 1
                continue
            if len(split_line) <= 1: continue
            colon1 = re.sub('  *', ' ', split_line[1].strip()).split(' ')
            assert len(colon1) >= 1
            r = {}
            r['name'] = colon1[0]
            if len(colon1) >= 6:
                r['minz'] = int(colon1[1])
                r['maxz'] = int(colon1[2])
                r['minang'] = int(colon1[3])
                r['maxang'] = int(colon1[4])
                n_extra_ = int(colon1[5])
            else:
                n_extra_ = 0
            if n_extra_ > 0:
                r['extra'] = []
            if len(colon1) >= 7: r['type'] = colon1[6]
            r['keys'] = split_line[0].strip().split(' ')
            if len(split_line) > 2:
                r['contexts'] = re.sub('  *', ' ', split_line[2].strip(' ')).split(' ')
            if len(split_line) > 3:
                r['attributes'] = re.sub('  *', ' ', split_line[3].strip(' ')).split(' ')
            if len(split_line) > 3:
                r['comment'] = split_line[4].strip(' ')
            if element and (len(colon1) < 6 or periodic_table.index(element) + 1 < r['minz'] or periodic_table.index(
                    element) + 1 > r['maxz']): continue
            if key and (len(colon1) < 6 or not any([key.lower() == k.lower() for k in r['keys']])): continue
            if context and (
                    'contexts' not in r or not any([context.lower() == k.lower() for k in r['contexts']])): continue
            if not context and 'contexts' in r and not any(['orbital' == k.lower() for k in r['contexts']]): continue
            if type and 'type' in r and type != r['type']: continue
            result.append(r)
            n_extra = n_extra_
        return result


periodic_table = [
    "H", "He",
    "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni",
    "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
    "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt",
    "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu",
    "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
    "Ua", "Ub", "Uc", "Ud", "Ue", "Uf", "Ug", "Uh",
]

small_core_ranges = [
    (1, 4),
    (11, 12),
    (19, 30),
    (37, 48),
    (55, 80),
    (87, 112),
]


def element(given) -> str:
    """
    Return the properly-cased chemical element symbol corresponding to a given atomic number or case-insensitive chemical element symbol
    """
    if isinstance(given, str):
        return given[0].upper() + given[1:].lower()
    elif type(given) is int:
        return periodic_table[given - 1]
    else:
        raise ValueError


def atomic_number(given) -> int:
    """
    Return the atomic number corresponding to a given chemical element symbol
    """
    if type(given) is int:
        return given
    elif isinstance(given, str):
        return periodic_table.index(element(given)) + 1
    else:
        raise ValueError


def core_correlation_segment(element):
    z = atomic_number(element)
    core_correlation_boundaries = ([range[0] - 1 for range in small_core_ranges] + [range[1] for range in
                                                                                    small_core_ranges]) + [9999]
    core_correlation_boundaries.sort()
    for segment in range(len(core_correlation_boundaries) - 1):
        if z > core_correlation_boundaries[segment] and z <= core_correlation_boundaries[segment + 1]:
            return segment


def elements_from_xyz(xyz: str) -> list[str]:
    """
    Extract the unique chemical elements in an xyz file
    :param xyz: The xyz
    :return: A list of unique chemical elements
    """
    lines = xyz.strip().split('\n')
    try:
        if lines and int(lines[0].strip()) == len(lines) - 2:
            lines = lines[2:]
    except:
        pass
    elements = []
    for line in lines:
        elements.append(element(line.strip().split(' ')[0]))
    elements = list(set(elements))
    elements.sort()
    return elements


def element_ranges(elements: list, heavy: bool = True, light: bool = True) -> list[str]:
    """
    Format a list of chemical elements into ranges of contiguous elements not spanning a core-correlation segment boundary
    :param elements: Atomic numbers or chemical element symbols
    :param heavy: Whether to include 'heavy' elements, ie those for which core correlation is active in the mixed core-correlation scheme
    :param light: Whether to include 'light' elements, ie those for which core correlation is not active in the mixed core-correlation scheme
    """
    result = []
    this = []
    atomic_numbers = list(set([atomic_number(element) for element in elements]))
    atomic_numbers.sort()
    for k, atno in enumerate(atomic_numbers):
        if (heavy and core_correlation_segment(atno) % 2 == 0) or (light and core_correlation_segment(atno) % 2 == 1):
            if this and (core_correlation_segment(atno) != core_correlation_segment(this[-1]) or atno != this[-1] + 1):
                result.append(element(this[0]) + ('-' + element(this[-1]) if len(this) > 1 else ''))
                this = []
            this.append(atno)
    if this:
        result.append(element(this[0]) + ('-' + element(this[-1]) if len(this) > 1 else ''))
    return result


def mixed_core_correlation_assert(element_range: str, core_correlation: bool = True) -> bool:
    """
    Determine whether the given range of chemical elements overlaps with the set of elements for which in the 'mixed' core correlation model, core correlation is active or inactive, depending on option.
    :param element_range: A single element symbol, or a range given as, e.g., Li-Ne
    :param core_correlation: If True, then the function answers the question whether, for any element in the range of elements, core correlation is active in mixed core; if false, then the function answers the question whether, for any element in the range of elements, core correlation is inactive in mixed core
    """

    if isinstance(element_range, str) and '-' in element_range:
        start, end = element_range.split('-')
        start = atomic_number(start)
        end = atomic_number(end)
        if start > end:
            return False
        return any([mixed_core_correlation_assert(element, core_correlation) for element in range(start, end + 1)])
    else:
        element = atomic_number(element_range)
        return not core_correlation ^ any([element >= range[0] and element <= range[1] for range in small_core_ranges])

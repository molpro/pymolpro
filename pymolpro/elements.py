import re
from dataclasses import dataclass


@dataclass(frozen=True)
class Element:
    atomic_number: int

    def __str__(self):
        return self.symbol

    def __eq__(self, other):
        return self.atomic_number == other.atomic_number

    @property
    def symbol(self):
        return periodic_table[self.atomic_number - 1].capitalize()

    def __lt__(self, other):
        return self.atomic_number < other.atomic_number


def element(given: int | str | Element) -> Element:
    """
    Factory for Element class

    :param given: An atomic number, or case-insensitive element symbol, or an instance of Element
    """
    if type(given) is int:
        if given < 1 or given > len(periodic_table):
            raise ValueError
        return Element(given)
    elif isinstance(given, str):
        try:
            return Element(periodic_table.index(given.capitalize()) + 1)
        except ValueError:
            return element(int(given))
    elif isinstance(given, Element):
        return Element(given.atomic_number)
    else:
        raise ValueError


class ElementRange:
    def __init__(self, first, second=None):
        if isinstance(first, ElementRange):
            self.first = first.first
            self.second = first.second
        elif second is not None:
            self.first = element(first)
            self.second = element(second)
        elif isinstance(first, str) and '-' in first:
            self.first, self.second = first.split('-')
            self.first = element(self.first)
            self.second = element(self.second)
        else:
            self.first = element(first)
            self.second = self.first

    def __str__(self):
        if self.first == self.second:
            return str(self.first)
        else:
            return str(self.first) + '-' + str(self.second)

    def __eq__(self, other):
        return self.first == other.first and self.second == other.second


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
    # the ranges of elements for which, in the mixed core-correlation scheme, core correlation is active, ie large != mixed in https://www.molpro.net/manual/doku.php?id=general_program_structure#defining_orbital_subspaces_occ_closed_core_frozen
    (3, 4),
    (11, 12),
    (19, 126),
]


def core_correlation_segment(el: Element | str | int):
    element_ = element(el)
    core_correlation_boundaries = ([range[0] - 1 for range in small_core_ranges] + [range[1] for range in
                                                                                    small_core_ranges]) + [
                                      len(periodic_table)]
    core_correlation_boundaries.sort()
    for segment in range(len(core_correlation_boundaries) - 1):
        if element_.atomic_number > core_correlation_boundaries[segment] and element_.atomic_number <= \
                core_correlation_boundaries[segment + 1]:
            return segment
    return 0


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
        split_ = re.sub('[0-9]*$','',line.strip().split(' ')[0])
        instance = element(split_)
        elements.append(instance)
    elements = list(set(elements))
    elements.sort()
    return elements


def element_ranges(elements: list[Element], heavy: bool = True, light: bool = True) -> list[ElementRange]:
    """
    Format a list of chemical elements into ranges of contiguous elements not spanning a core-correlation segment boundary
    :param elements:
    :param heavy: Whether to include 'heavy' elements, ie those for which core correlation is active in the mixed core-correlation scheme
    :param light: Whether to include 'light' elements, ie those for which core correlation is not active in the mixed core-correlation scheme
    """
    result = []
    this = []
    atomic_numbers = list(set([element(instance).atomic_number for instance in elements]))
    atomic_numbers.sort()
    for k, atno in enumerate(atomic_numbers):
        if (heavy and core_correlation_segment(atno) % 2 == 0) or (light and core_correlation_segment(atno) % 2 == 1):
            if this and (core_correlation_segment(atno) != core_correlation_segment(this[-1]) or atno != this[-1] + 1):
                result.append(ElementRange(this[0], this[-1]))
                this = []
            this.append(atno)
    if this:
        result.append(ElementRange(this[0], this[-1]))
    return result


def mixed_core_correlation_assert(element_range: ElementRange | str, core_correlation: bool = True) -> bool:
    """
    Determine whether the given range of chemical elements overlaps with the set of elements for which in the 'mixed' core correlation model, core correlation is active or inactive, depending on option.
    :param element_range: A single element symbol, or a range given as, e.g., Li-Ne
    :param core_correlation: If True, then the function answers the question whether, for any element in the range of elements, core correlation is active in mixed core; if false, then the function answers the question whether, for any element in the range of elements, core correlation is inactive in mixed core
    """

    start = ElementRange(element_range).first.atomic_number
    end = ElementRange(element_range).second.atomic_number
    if start == end:
        return not core_correlation ^ any([start >= range[0] and start <= range[1] for range in small_core_ranges])
    elif start < end:
        return any([mixed_core_correlation_assert(element, core_correlation) for element in range(start, end + 1)])
    else:
        return False


def elements_include_heavy(elements: list[str | int | Element]) -> bool:
    """
    Whether any of the provided elements lies within the set that have core correlation active in the mixed core correlation scheme

    :param elements:
    """
    result = False
    for el in elements:
        result |= mixed_core_correlation_assert(element(el), True)
    return result

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

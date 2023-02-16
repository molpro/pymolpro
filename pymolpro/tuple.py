import numpy as np
from pysjef import xpath
from pymolpro.orbital import Orbital
from lxml import etree


class Tuple:
    """
    Container for a tuple of orbitals
    """

    def __init__(self, node):
        """
        Initialise from a node on a Molpro output xml tree

        :param node: lxml.etree.Element holding a correlation single or pair descriptor
        """
        self.node = node
        orbital_labels = [(node.get("orbital1") if node.get("orbital") is None else node.get("orbital"))]
        orbital2 = node.get("orbital2")
        if orbital2 is not None:
            orbital_labels.append(orbital2)
        orbital3 = node.get("orbital3")
        if orbital3 is not None:
            orbital_labels.append(orbital3)
        self.len_ = len(orbital_labels)
        self.spins = [1 for orbital in orbital_labels]
        for particle in range(self.len_):
            if orbital_labels[particle][0] == '-':
                self.spins[particle] = -1
                orbital_labels[particle] = orbital_labels[particle][1:]
        # TODO UHF case

        # look for orbitals preceding
        self.orbitals = []
        for orbital_label in orbital_labels:
            xpath_results = xpath(self.node,
                                  '../preceding-sibling::molecule/orbitals/orbital[@ID="' + orbital_label + '"]')
            if len(xpath_results) == 0:
                break
            self.orbitals.append(Orbital(xpath_results[-1]))

        if len(self.orbitals) < self.len_:
            # look for orbitals following
            self.orbitals = []
            for orbital_label in orbital_labels:
                xpath_results = xpath(self.node,
                                      '../following-sibling::molecule/orbitals/orbital[@ID="' + orbital_label + '"]')
                if len(xpath_results) == 0:
                    raise Exception("xml output does not contain orbitals for pair")
                self.orbitals.append(Orbital(xpath_results[0]))

        self.energy = None if self.node.get('energy') is None else float(self.node.get('energy'))

    def __len__(self):
        return self.len_

class Pair(Tuple):
    """
    Container for a pair of orbitals
    """


class Single(Tuple):
    """
    Container for a single orbital
    """

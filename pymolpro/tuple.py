import numpy as np
from pysjef import xpath
from pymolpro.orbital import Orbital


class Tuple:
    """
    Container for a tuple of orbitals

    :param node: Node holding a correlation single or pair descriptor
    :type node: lxml.etree.Element
    """

    def __init__(self, node):
        """
        Initialise from a node on a Molpro output xml tree

        """
        self.node = node
        self.spins = [1] #: 1(alpha) or -1(beta) for each orbital
        orbital1 = node.get("orbital1") if node.get("orbital") is None else node.get("orbital")
        if orbital1[0] == '-':
            self.spins[0] = -1
            orbital1 = orbital1[1:]
        orbital2 = node.get("orbital2")
        if orbital2 is None:
            orbital2 = "nomatch"
        else:
            self.spins.append(1)
            if orbital2[0] == '-':
                self.spins[-1] = -1
                orbital2 = orbital2[1:]
        orbital3 = node.get("orbital3")
        if orbital3 is None:
            orbital3 = "nomatch"
        else:
            self.spins.append(1)
            if orbital3[0] == '-':
                self.spins[-1] = -1
                orbital3 = orbital3[1:]
        # TODO UHF case
        self.orbitals = [Orbital(node) for node in xpath(self.node,
                                                         '/molpro//orbitals/orbital[@ID="' + orbital1 + '" or @ID="' + orbital2 + '" or @ID="' + orbital3 + '"]')] #: Orbital objects forming the tuple
        if orbital2 == orbital1: self.orbitals.append(self.orbitals[0])
        self.energy = None if self.node.get('energy') is None else float(self.node.get('energy')) #: Correlation energy contribution from a Single or Pair

class Pair(Tuple):
    """
    Container for a pair of orbitals
    """


class Single(Tuple):
    """
    Container for a single orbital
    """

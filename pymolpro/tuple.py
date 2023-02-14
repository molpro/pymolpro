import numpy as np
import scipy as sp
import math
import scipy.special
from lxml import etree

import pymolpro
from pymolpro.orbital import Orbital


class Tuple:
    """
    Container for a pair of orbitals
    """

    def __init__(self, node):
        """
        Initialise from a node on a Molpro output xml tree

        :param node: lxml.etree.Element holding a single orbital
        """
        self.node = node
        self.spins = [1]
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
        self.orbitals = [Orbital(node) for node in pymolpro.xpath(self.node,
                                                                  '/molpro//orbitals/orbital[@ID="' + orbital1 + '" or @ID="' + orbital2 + '" or @ID="' + orbital3 + '"]')]
        self.energy = float(self.node.get('energy'))

    def data(self, grid=2):
        """
        Construct a set of floating point data describing the orbital tuple together with the tuple energy

        :param grid: Number of orbital grid points in each direction
        :return: flat numpy array
        """
        result = np.array([], dtype=float)
        for particle in range(len(self.spins)):
            result = np.append(result, self.spins[particle])
            orbital = self.orbitals[particle]
            orbital_grid = np.reshape(orbital.grid(grid), (-1, 4))[:, :3]
            print("orbital_grid: ", orbital_grid)
            result = np.append(result, orbital_grid)
            for other_particle in range(particle + 1, len(self.spins)):
                other_orbital = self.orbitals[other_particle]
                other_orbital_grid = np.reshape(other_orbital.grid(grid), (-1, 4))[:, :3]
                for point in orbital_grid:
                    for other_point in other_orbital_grid:
                        dist = np.linalg.norm(point - other_point)
                        result = np.append(result, 1 / dist)
        return [result, self.energy]


class Pair(Tuple):
    """
    Container for a pair of orbitals
    """


class Single(Tuple):
    """
    Container for a single orbital
    """

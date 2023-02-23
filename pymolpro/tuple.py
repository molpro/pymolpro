import numpy as np
from pysjef import xpath
from pymolpro.orbital import Orbital
from lxml import etree


class Tuple:
    r"""
    Container for a tuple of orbitals

    :param node: Node holding a correlation single or pair descriptor
    :type node: lxml.etree.Element
    """

    def __init__(self, node):
        """
        Initialise from a node on a Molpro output xml tree

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
        self.spins = [1 for orbital in orbital_labels]  #: 1(alpha) or -1(beta) for each orbital
        for particle in range(self.len_):
            if orbital_labels[particle][0] == '-':
                self.spins[particle] = -1
                orbital_labels[particle] = orbital_labels[particle][1:]
        # TODO UHF case

        # look for orbitals preceding
        self.orbitals = []  #: Orbital objects forming the tuple
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

        self.energy = None if self.node.get('energy') is None else float(
            self.node.get('energy'))  #: Correlation energy contribution from a Single or Pair

    def __len__(self):
        return self.len_


class Pair(Tuple):
    """
    Container for a pair of orbitals
    """

    @property
    def axes(self):
        r"""
        The coordinate system for the pair of orbitals.  The :math:`z` axis is the vector
        :math:`\vec e_{z,10}=(\vec r_1-\vec r_0)/|\vec r_1-\vec r_0|` from
        orbital 0 to orbital 1.
        The :math:`x` axis is the along the vector
        :math:`\vec e_{x,10} = \vec e_{z,10} \times \vec e_z` (or :math:`\vec e_y` if necessary), normalised.
        Then :math:`\vec e_{y,10}=\vec e_{z,10}\times\vec e_{x,10}`.
        If :math:`|\vec r_1-\vec r_0|` is small, the global coordinate axes are adopted.

        :return: The axes, with the first index labelling the component in the base coordinate system, and the second index specifying which pair axis.
        :rtype: np.array(3,3)
        """
        axes = np.empty([3, 3], dtype=float)
        axes[:, 2] = self.orbitals[1].centroid - self.orbitals[0].centroid
        if np.linalg.norm(axes[:, 2]) < 1e-8: return np.identity(3)
        axes[:, 2] = axes[:, 2] / np.linalg.norm(axes[:, 2])
        axes[0, 0] = axes[1, 2]
        axes[1, 0] = -axes[0, 2]
        axes[2, 0] = 0
        if np.linalg.norm(axes[:, 0]) < 1e-6:
            axes[0, 0] = axes[2, 2]
            axes[1, 0] = 0
            axes[2, 0] = -axes[0, 2]
        axes[:, 0] = axes[:, 0] / np.linalg.norm(axes[:, 0])
        axes[:, 1] = np.cross(axes[:, 2], axes[:, 0])
        return axes


class Single(Tuple):
    """
    Container for a single orbital
    """

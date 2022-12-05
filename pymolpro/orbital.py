import numpy as np
import scipy as sp
import math


class Orbital:
    """
    Container for an orbital (usually molecular).
    """

    @property
    def local_second_moments(self):
        global_second_moments = [float(self.node['moments'].split()[3+k]) for k in range(6)]
        second_moments_matrix = np.zeros((3, 3))
        second_moments_matrix[0][0] = global_second_moments[0]
        second_moments_matrix[0][1] = global_second_moments[3]
        second_moments_matrix[0][2] = global_second_moments[4]
        second_moments_matrix[1][0] = global_second_moments[3]
        second_moments_matrix[1][1] = global_second_moments[1]
        second_moments_matrix[1][2] = global_second_moments[5]
        second_moments_matrix[2][0] = global_second_moments[4]
        second_moments_matrix[2][1] = global_second_moments[5]
        second_moments_matrix[2][2] = global_second_moments[2]
        second_moments_matrix -= np.outer(self.centroid,self.centroid)
        return second_moments_matrix

    @property
    def kinetic_energy(self):
        return float(self.node['moments'].split()[9])

    def attribute(self, key):
        return self.node[key]

    @property
    def ID(self):
        return self.attribute('ID')

    def global_to_local(self, global_coordinates):
        return self.second_moment_eigenvectors.transpose() * (global_coordinates - self.centroid)

    def local_to_global(self, local_coordinates):
        return self.centroid + self.second_moment_eigenvectors * local_coordinates

    def grid(self, npt, method='erfinv'):
        if type(npt) != list:
            npt = [npt for k in range(3)]
        if method == 'erfinv':
            points = [
                [sp.special.erfinv(k / float(npt[i] + 1) - 0.5) / math.sqrt(self.second_moment_eigenvalues[i]) for k in
             range(npt[i])] for i in range(3)]
        else:
            assert False
        return self.local_to_global(np.array(points))

    def __init__(self, node):
        """
        Initialise from a node on a Molpro output xml tree

        :param node: lxml.etree.Element holding a single orbital
        """
        self.node = node
        self.occupation = float(node['occupation'])
        self.centroid = [float(self.node['moments'].split()[k]) for k in range(3)]
        self.second_moment_eigenvalues, self.second_moment_eigenvectors = np.linalg.eigh(self.local_second_moments)
        self.coefficients = node['#text']


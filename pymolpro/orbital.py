import numpy as np
import scipy as sp
import math
import scipy.special

import pymolpro.grid


class Orbital:
    """
    Container for an orbital (usually molecular).
    """

    @property
    def kinetic_energy(self):
        return float(self.attribute('moments').split()[9])

    def attribute(self, key):
        return self.node.get(key)

    @property
    def ID(self):
        return self.attribute('ID')

    def grid(self, npt, method='erfinv', integration_weights=False, scale=1.0):
        """
        Generate a grid centred on the orbital.

        :param npt: Number of desired points in each coordinate. If a single number is specified, it is used for all three directions.
        :param method: Algorithm for grid generation.
        :param integration_weights: Whether integration weights are required.
        :param scale: Scale the grid by this factor.
        :return: if integration_weights, a tuple containing the points (numpy array [det(npt),3]) and weights (numpy array). Otherwise, just the points.
        """
        if type(npt) != list:
            npt = [npt for k in range(3)]
        points = []
        weights = []
        if method == 'erfinv':
            points = [
                [sp.special.erfinv(2 * (k + 1) / float(npt[i] + 1) - 1) * scale * math.sqrt(
                    self.second_moment_eigenvalues[i])
                 for k in range(npt[i])] for i in range(3)]
            weights = [[1.0 / npt[i] for k in range(npt[i])] for i in range(3)]
        elif method == 'Gauss-Hermite':
            for i in range(3):
                x, w = scipy.special.roots_hermite(npt[i])
                betamh = scale * math.sqrt(2 * self.second_moment_eigenvalues[i])
                points.append([x[k] * betamh for k in range(npt[i])])
                weights.append([math.exp(x[k] * x[k]) * w[k] * betamh for k in range(npt[i])])
        else:
            assert False
        points3d = np.empty([npt[0] * npt[1] * npt[2], 3])
        weights3d = np.ones(npt[0] * npt[1] * npt[2])
        for j in range(npt[0]):
            for k in range(npt[1]):
                for l in range(npt[2]):
                    points3d[j + npt[0] * (k + npt[1] * l), 0] = points[0][j]
                    points3d[j + npt[0] * (k + npt[1] * l), 1] = points[1][k]
                    points3d[j + npt[0] * (k + npt[1] * l), 2] = points[2][l]
                    weights3d[j + npt[0] * (k + npt[1] * l)] = weights[0][j] * weights[1][k] * weights[2][l]
        global_points = self.__local_to_global(points3d)
        return [global_points, weights3d] if integration_weights else global_points

    def evaluate(self, points, values=False):
        """
        Evaluate orbital on a grid of points

        :param points: List of geometries specified as 3-list of values in bohr
        :param values:
        :return: array of dictionaries giving the occupation and values on the grid, or if ID is specified, a single dictionary, or if values==True, a numpy array
        """
        return pymolpro.grid.evaluateOrbitals(self.node.xpath('./parent::*/parent::*')[-1], points, ID=self.ID,
                                              values=values)
    def __init__(self, node):
        """
        Initialise from a node on a Molpro output xml tree

        :param node: lxml.etree.Element holding a single orbital
        """
        self.node = node
        self.occupation = float(self.attribute('occupation'))
        self.centroid = [float(self.attribute('moments').split()[k]) for k in range(3)]
        self.second_moment_eigenvalues, self.second_moment_eigenvectors = np.linalg.eigh(self.local_second_moments)
        self.coefficients = self.node.text.split()
        self.coefficients = [float(c) for c in self.coefficients]

    @property
    def local_second_moments(self):
        global_second_moments = [float(self.attribute('moments').split()[3 + k]) for k in range(6)]
        second_moments_matrix = np.zeros((3, 3))
        second_moments_matrix[0][0] = global_second_moments[0]
        second_moments_matrix[1][0] = second_moments_matrix[0][1] = global_second_moments[3]
        second_moments_matrix[2][0] = second_moments_matrix[0][2] = global_second_moments[4]
        second_moments_matrix[1][1] = global_second_moments[1]
        second_moments_matrix[1][2] = second_moments_matrix[2][1] = global_second_moments[5]
        second_moments_matrix[2][2] = global_second_moments[2]
        # print("raw second moments matrix", second_moments_matrix)
        second_moments_matrix -= np.outer(self.centroid, self.centroid)
        return second_moments_matrix

    def __local_to_global(self, local_coordinates):
        return np.array([self.centroid for k in local_coordinates[:, 0]]) + np.matmul(
            self.second_moment_eigenvectors, local_coordinates.transpose()).transpose()
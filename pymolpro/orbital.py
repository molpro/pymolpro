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
            _npt = [npt for k in range(3)]
        else:
            _npt = npt
        points = []
        weights = []
        # grids for second moment eigenvalues unity
        if method == 'erfinv':
            points = [
                [sp.special.erfinv(2 * (k + 1) / float(_npt[i] + 1) - 1)
                 for k in range(_npt[i])] for i in range(3)]
            weights = [[1.0 / _npt[i] for k in range(_npt[i])] for i in range(3)]
        elif method == 'Gauss-Hermite':
            for i in range(3):
                x, w = scipy.special.roots_hermite(_npt[i])
                points.append([x[k] * math.sqrt(2) for k in range(_npt[i])])
                weights.append([math.exp(x[k] * x[k]) * w[k] * math.sqrt(2) for k in range(_npt[i])])
        elif method == 'Gauss-Laguerre-Lebedev':
            import scipy.special
            nptang = self.__lebedev_select(_npt[1], True)
            import numgrid
            gridang = numgrid.get_angular_grid(nptang)
            npt3 = _npt[0] * nptang
            # print("radial size:", _npt[0], " angular size:", nptang, " total:", npt3)
            points3d = np.empty([npt3, 3])
            weights3d = np.ones(npt3)
            radial_points, radial_weights = scipy.special.roots_laguerre(_npt[0])
            for i in range(len(radial_weights)):
                radial_weights[i] *= pow(radial_points[i], 2) * math.exp(radial_points[i]) * 0.5 * np.pi
            radial_points *= 0.5
            for j in range(_npt[0]):
                for k in range(nptang):
                    for i in range(3):
                        points3d[j + _npt[0] * k, i] = radial_points[j] * gridang[i][k]
                    weights3d[j + _npt[0] * k] = radial_weights[j] * gridang[3][k]
        else:
            assert False

        if len(points) > 0:
            # expand outer product grid to three dimensions
            npt3 = _npt[0] * _npt[1] * _npt[2]
            points3d = np.empty([npt3, 3])
            weights3d = np.ones(npt3)
            for j in range(_npt[0]):
                for k in range(_npt[1]):
                    for l in range(_npt[2]):
                        points3d[j + _npt[0] * (k + _npt[1] * l), 0] = points[0][j]
                        points3d[j + _npt[0] * (k + _npt[1] * l), 1] = points[1][k]
                        points3d[j + _npt[0] * (k + _npt[1] * l), 2] = points[2][l]
                        weights3d[j + _npt[0] * (k + _npt[1] * l)] = weights[0][j] * weights[1][k] * weights[2][l]

        coordinate_scaling = np.array([scale * math.sqrt(e) for e in self.second_moment_eigenvalues])
        # print("coordinate scaling:", coordinate_scaling)
        jacobian = coordinate_scaling[0] * coordinate_scaling[1] * coordinate_scaling[2]
        weights3d = jacobian * weights3d
        global_points = self.__local_to_global(points3d, coordinate_scaling)
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
        second_moments_matrix -= np.outer(self.centroid, self.centroid)
        return second_moments_matrix

    def __local_to_global(self, local_coordinates, scaling=[1., 1., 1.]):
        return np.array([self.centroid for k in local_coordinates[:, 0]]) + np.matmul(
            np.matmul(self.second_moment_eigenvectors, np.diagflat(scaling)), local_coordinates.transpose()).transpose()

    def __lebedev_select(self, l, size=False):
        lebedev = {3: 6, 5: 14, 7: 26, 9: 38, 11: 50, 13: 74, 15: 86, 17: 110, 19: 146,
                   21: 170, 23: 194, 25: 230, 27: 266, 29: 302, 31: 350, 35: 434, 41: 590, 47: 770,
                   53: 974, 59: 1202, 65: 1454, 71: 1730, 77: 2030, 83: 2354, 89: 2702, 95: 3074, 101: 3470,
                   107: 3890, 113: 4334, 119: 4802, 125: 5294, 131: 5810}
        for k, v in lebedev.items():
            if l <= k: return v if size else k
        return lebedev[131] if size else 131

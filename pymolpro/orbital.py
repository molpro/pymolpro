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

    def grid(self, npt, method='erfinv', integration_weights=False, scale=1.0, grid_parameters=[],
             spherical_average=False):
        """
        Generate a grid centred on the orbital.

        :param npt: Number of desired points in each coordinate.
        :param method: Algorithm for grid generation.
        :param integration_weights: Whether integration weights are required.
        :param scale: Scale the grid by this factor.
        :return: if integration_weights, a tuple containing the points (numpy array [det(npt),3]) and weights (numpy array). Otherwise, just the points.
        """
        # grids for second moment eigenvalues unity
        if method == 'erfinv':
            assert type(npt) != list
            points = [sp.special.erfinv(2 * (k + 1) / float(npt + 1) - 1) for k in range(npt)]
            weights = [1.0 / npt for k in range(npt)]
            points3d, weights3d = pymolpro.grid.cubical_grid(points, weights)
        elif method == 'Gauss-Hermite':
            assert type(npt) != list
            import scipy.special
            x, w = scipy.special.roots_hermite(npt)
            points = [x[k] * math.sqrt(2) for k in range(npt)]
            weights = [math.exp(x[k] * x[k]) * w[k] * math.sqrt(2) for k in range(npt)]
            points3d, weights3d = pymolpro.grid.cubical_grid(points, weights)
        elif 'Lebedev' in method:
            import scipy.special
            if 'Laguerre' in method:
                radial_points, radial_weights = scipy.special.roots_laguerre(npt[0] if type(npt) == list else npt)
                for i in range(len(radial_weights)):
                    radial_weights[i] *= math.exp(radial_points[i]) / 2
                radial_points *= 0.5  # why?
            elif 'Mura' in method:
                n1 = npt[0] if type(npt) == list else npt
                m = grid_parameters[0] if len(grid_parameters) > 1 else 3
                xpoints = [(float(i) + 0.5) / n1 for i in range(n1)]
                radial_weights = [m * scale * pow(x, m - 1) / (1 - pow(x, m)) / n1 for x in xpoints]
                radial_points = np.array([- scale * math.log(1 - pow(x, m)) for x in xpoints])
            else:
                assert False
            points3d, weights3d = pymolpro.grid.spherical_grid(radial_points, radial_weights,
                                                               npt[1] if type(npt) == list and len(npt) > 1 else len(
                                                                   radial_weights))
        else:
            assert False

        coordinate_scaling = np.array([scale * math.sqrt(e) for e in self.second_moment_eigenvalues])
        if spherical_average and 'Lebedev' in method:
            coordinate_scaling[:] = scale * pow(
                self.second_moment_eigenvalues[0] * self.second_moment_eigenvalues[1] * self.second_moment_eigenvalues[
                    2],
                1.0 / 6.0)
        # print("before __local_to_global points3d", points3d)
        weights3d = coordinate_scaling[0] * coordinate_scaling[1] * coordinate_scaling[2] * weights3d
        global_points = self.__local_to_global(points3d, coordinate_scaling)
        return tuple([global_points, weights3d]) if integration_weights else global_points

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

import pathlib

import numpy as np
import scipy as sp
import math
import scipy.special

import pymolpro.grid
from . import sparse_dump

class Orbital:
    """
    Container for an orbital (usually molecular).
    """

    @property
    def kinetic_energy(self):
        """
        Kinetic energy expectation value for the orbital.


        """
        return float(self.attribute('moments').split()[9])

    def attribute(self, key):
        return self.node.get(key)

    @property
    def ID(self):
        return self.attribute('ID')

    def grid(self, npt, method='erfinv', scale=1.0, grid_parameters=[],
             spherical_average=False):
        """
        Generate a grid centred on the orbital.

        :param npt: Number of desired points in each coordinate.
        :param method: Algorithm for grid generation.
        :param scale: Scale the grid by this factor.
        :return: points and weights (numpy array [npt,4])
        """
        # grids for second moment eigenvalues unity
        if method == 'erfinv':
            assert type(npt) is not list
            points = [sp.special.erfinv(2 * (k + 1) / float(npt + 1) - 1) for k in range(npt)]
            weights = [1.0 / npt for k in range(npt)]
            points3d = pymolpro.grid.cubical_grid(points, weights)
        elif method == 'Gauss-Hermite':
            assert type(npt) is not list
            import scipy.special
            x, w = scipy.special.roots_hermite(npt)
            points = [x[k] * math.sqrt(2) for k in range(npt)]
            weights = [math.exp(x[k] * x[k]) * w[k] * math.sqrt(2) for k in range(npt)]
            points3d = pymolpro.grid.cubical_grid(points, weights)
        elif 'Lebedev' in method:
            import scipy.special
            if 'Laguerre' in method:
                radial_points, radial_weights = scipy.special.roots_laguerre(npt[0] if type(npt) is list else npt)
                for i in range(len(radial_weights)):
                    radial_weights[i] *= math.exp(radial_points[i]) / 2
                radial_points *= 0.5  # why?
            elif 'Mura' in method:
                scalem = grid_parameters[1] if len(grid_parameters) > 1 else 10
                n1 = npt[0] if type(npt) is list else npt
                m = grid_parameters[0] if len(grid_parameters) > 0 else 3
                xpoints = [(float(i) + 0.5) / n1 for i in range(n1)]
                radial_weights = [m * scalem * pow(x, m - 1) / (1 - pow(x, m)) / n1 for x in xpoints]
                radial_points = np.array([- scalem * math.log(1 - pow(x, m)) for x in xpoints])
            else:
                assert False
            points3d = pymolpro.grid.spherical_grid(radial_points, radial_weights,
                                                    npt[1] if type(npt) is list and len(npt) > 1 else len(
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
        points3d[:, 3] = coordinate_scaling[0] * coordinate_scaling[1] * coordinate_scaling[2] * points3d[:, 3]
        points3d[:, :3] = self.__local_to_global(points3d, coordinate_scaling)
        return points3d

    def evaluate(self, points, values=False):
        """
        Evaluate orbital on a grid of points

        :param points: List of geometries specified as 3-list of values in bohr
        :param values:
        :return: array of dictionaries giving the occupation and values on the grid, or if ID is specified, a single dictionary, or if values==True, a numpy array
        """
        return pymolpro.grid.evaluateOrbitals(self.node.xpath('./parent::*/parent::*')[-1], points, ID=self.ID,
                                              values=values)

    def __init__(self, node, directory=None):
        """
        Initialise from a node on a Molpro output xml tree

        :param node: lxml.etree.Element holding a single orbital
        :param directory: the directory in which the xml file, and its sidecar, live
        """
        self.node = node
        if self.attribute('energy') not in ['',None]:
            self.energy = float(self.attribute('energy'))  #: energy of the orbital
        self.occupation = float(self.attribute('occupation'))  #: Occupation of the orbital
        self.centroid = np.array(
            [float(self.attribute('moments').split()[k]) for k in range(3)])  #: Centroid of the orbital
        #: Eigenvalues of the orbital second-moment tensor (origin at centre of charge) in ascending order
        self.second_moment_eigenvalues, self._second_moment_eigenvectors = np.linalg.eigh(self.local_second_moments)
        for i in range(3):
            if max([abs(c) for c in self._second_moment_eigenvectors[:, i]]) < 0:
                self._second_moment_eigenvectors[:, i] *= -1
        coefficients_text = self.node.text
        if coefficients_text is not None and coefficients_text.strip() != '':
            self.coefficients = coefficients_text.split()
            self.coefficients = np.array([float(c) for c in self.coefficients])
        elif self.attribute('sidecar_offset') not in (None, '') and directory is not None:
            sidecar_file = pathlib.Path(directory) / pathlib.Path(node.getparent().get('sidecar'))
            self.coefficients,nothing = sparse_dump.sparse_dump_get(sidecar_file,int(self.attribute('sidecar_offset')))

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
            np.matmul(self.axes, np.diagflat(scaling)),
            local_coordinates[:, :3].transpose()).transpose()

    @property
    def axes(self):
        r"""
        The coordinate system for the orbital.  The :math:`z` axis is the normalised vector
        :math:`\vec e_{z}` that is the eigenvector of the second-moment tensor of largest eigenvalue,
        similarly for :math:`\vec e_y, \vec e_x`.
        The phase of each axis is chosen such that the largest of the coefficients specifying the axis in the global coordinate system is positive.

        :return: The axes, with the first index labelling the component in the base coordinate system, and the second index specifying which orbital axis.
        :rtype: np.array(3,3)
        """
        return self._second_moment_eigenvectors

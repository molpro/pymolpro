import numpy as np
import scipy as sp
import math

import scipy.special


class Orbital:
    """
    Container for an orbital (usually molecular).
    """

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

    @property
    def kinetic_energy(self):
        return float(self.attribute('moments').split()[9])

    def attribute(self, key):
        return self.node.get(key)

    @property
    def ID(self):
        return self.attribute('ID')

    def global_to_local(self, global_coordinates):
        return np.matmul(self.second_moment_eigenvectors.transpose(),
                         global_coordinates - np.array([self.centroid for k in global_coordinates[0, :]]).transpose())

    def local_to_global(self, local_coordinates):
        # print("local_to_global centroid shift", np.array([self.centroid for k in local_coordinates[0, :]]).transpose())
        # print("local_coordinates", local_coordinates)
        # print("local_coordinates rotated", np.matmul(
        #     self.second_moment_eigenvectors, local_coordinates))
        # print("local_coordinates rotated and shifted",
        #       np.array([self.centroid for k in local_coordinates[0, :]]).transpose() + np.matmul(
        #           self.second_moment_eigenvectors.transpose(), local_coordinates))
        # print("local_coordinates rotated and shifted",
        #       np.array([self.centroid for k in local_coordinates[0, :]]).transpose() + np.matmul(
        #           self.second_moment_eigenvectors, local_coordinates))
        return np.array([self.centroid for k in local_coordinates[0, :]]).transpose() + np.matmul(
            self.second_moment_eigenvectors, local_coordinates)

    def grid(self, npt, method='erfinv', integration_weights=False, scale=1.0):
        if type(npt) != list:
            npt = [npt for k in range(3)]
        points = []
        weights = []
        if method == 'erfinv':
            points = [
                [sp.special.erfinv(k / float(npt[i] + 1) - 0.5) * scale * math.sqrt(self.second_moment_eigenvalues[i])
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
        # print("local frame 1D points", points)
        # print("weights", weights)
        points3d = np.empty([3, npt[0] * npt[1] * npt[2]])
        weights3d = np.ones(npt[0] * npt[1] * npt[2])
        for j in range(npt[0]):
            for k in range(npt[1]):
                for l in range(npt[2]):
                    points3d[0, j + npt[0] * (k + npt[1] * l)] = points[0][j]
                    points3d[1, j + npt[0] * (k + npt[1] * l)] = points[1][k]
                    points3d[2, j + npt[0] * (k + npt[1] * l)] = points[2][l]
                    weights3d[j + npt[0] * (k + npt[1] * l)] = weights[0][j] * weights[1][k] * weights[2][l]
        # print("local frame 3D points", points3d)
        global_points = self.local_to_global(points3d)
        # print("global frame 3D points", global_points)
        result = [[global_points[i, j] for i in range(3)] for j in range(len(weights3d))]
        return [result, [w for w in weights3d]] if integration_weights else result

    def __init__(self, node):
        """
        Initialise from a node on a Molpro output xml tree

        :param node: lxml.etree.Element holding a single orbital
        """
        self.node = node
        self.occupation = float(self.attribute('occupation'))
        self.centroid = [float(self.attribute('moments').split()[k]) for k in range(3)]
        print("centroid in orbital.__init__", self.centroid)
        self.second_moment_eigenvalues, self.second_moment_eigenvectors = np.linalg.eigh(self.local_second_moments)
        self.coefficients = self.node.text.split()
        self.coefficients = [float(c) for c in self.coefficients]

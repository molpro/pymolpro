import pathlib

import numpy as np
import scipy as sp
import math
import scipy.special

import pymolpro.grid
from . import sparse_dump
from .cube_data import CubeData
from . elements import element
from . grid import cuboidal_grid


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
             spherical_average=False,resolution=None, orbital_cutoff=.001):
        """
        Generate a grid centred on the orbital.

        :param npt: Number of desired points in each coordinate.
        :param method: Algorithm for grid generation.
        :param scale: Scale the grid by this factor.
        :param resolution: Resolution of the grid in Bohr, for the uniform grid method.
        :return: points and weights (numpy array [npt,4])
        """
        if method == 'uniform':
            raise NotImplementedError("uniform grid not implemented yet")
            if resolution is None: resolution = 0.2
            # first construct a coarse grid to eliminate space where the orbital is everywhere less than orbital_cutoff
            coarse_resolution = 0.6
            cutoff = 0.001
            #initially expand to make a box that contains some orbital amplitude at all
            for nscale in range(1,100):
                bounds = np.array([[-nscale*coarse_resolution, nscale*coarse_resolution] for k in range(3)])
                points = []
                for i in range(3):
                    points.append(np.linspace(bounds[i][0],bounds[i][1],round(bounds[i][1]-bounds[i][0])/coarse_resolution)+1)
                grid = cuboidal_grid([points[0][:], points[1][:], points[2][:]])
                values = self.evaluate(grid,values=True)
                if np.max(values) > 0.01:
                    break
            # now extend the box in all directions until the orbital value is less than the cutoff on all faces
            while True:
                points = []
                for i in range(3):
                    points.append(np.linspace(bounds[i][0],bounds[i][1],round(bounds[i][1]-bounds[i][0])/coarse_resolution)+1)
                # +x
                grid = cuboidal_grid([points[0][-1], points[1][:], points[2][:]])
                values = self.evaluate(grid,values=True)
                if np.max(values) > cutoff:
                    bounds[0][1] += coarse_resolution
                    continue
                points3d = regular_grid(bounds, coarse_resolution)
                for axis in range(3):
                    for direction in range(2):
                        if np.linalg.norm(points3d[:, axis] - self.centroid[axis]) > orbital_cutoff:
                            bounds[axis][direction] = points3d[np.argmin(np.linalg.norm(points3d[:, axis] - self.centroid[axis], axis=1)), axis]
        # grids for second moment eigenvalues unity
        elif method == 'erfinv':
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

    @property
    def atoms(self):
        angstrom = 1.8897161646321
        atom_elements = self.node.xpath('ancestor::*[local-name() = "molecule"]/*[local-name() = "molecule"]//*[local-name() = "atom"]')
        geometry=[]
        for atom in atom_elements:
            type_ = atom.xpath('@elementType')[0]
            element_ = element(type_)
            geometry.append({'xyz': tuple([angstrom*float(*atom.xpath('@x3')), angstrom*float(*atom.xpath('@y3')), angstrom*float(*atom.xpath('@z3'))]),'atomic_number': element_.atomic_number, 'charge': element_.atomic_number})
        return geometry

    def evaluate(self, points, values=False):
        """
        Evaluate orbital on a grid of points

        :param points: List of geometries specified as 3-list of values in bohr
        :param values:
        :return: array of dictionaries giving the occupation and values on the grid, or if ID is specified, a single dictionary, or if values==True, a numpy array
        """
        return pymolpro.grid.evaluateOrbitals(self.node.xpath('./parent::*/parent::*')[-1], points, ID=self.ID,
                                              values=values, directory=self.directory)

    def __init__(self, node, directory=None):
        """
        Initialise from a node on a Molpro output xml tree

        :param node: lxml.etree.Element holding a single orbital
        :param directory: the directory in which the xml file, and its sidecar, live
        """
        self.node = node
        self.directory = directory
        try:
            self.energy = float(self.attribute('energy'))  #: energy of the orbital
        except:
            self.energy = 0.0
        try:
          self.occupation = float(self.attribute('occupation'))  #: Occupation of the orbital
        except:
          self.occupation = 2.0
        try:
          self.centroid = np.array(
              [float(self.attribute('moments').split()[k]) for k in range(3)])  #: Centroid of the orbital
        #: Eigenvalues of the orbital second-moment tensor (origin at centre of charge) in ascending order
          self.second_moment_eigenvalues, self._second_moment_eigenvectors = np.linalg.eigh(self.local_second_moments)
          for i in range(3):
              if max([abs(c) for c in self._second_moment_eigenvectors[:, i]]) < 0:
                  self._second_moment_eigenvectors[:, i] *= -1
        except:
          pass
        coefficients_text = self.node.text
        if coefficients_text is not None and coefficients_text.strip().replace('\n', '').strip() != '':
            coefficients_text = coefficients_text.strip().replace('\n', '').strip()
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

    def cube_data(self, resolution:float=.2, border=4.25)->CubeData:
        """
        Generates a 3D data cube representation of the molecular system.

        The method computes a 3D data grid based on the atomic coordinates in the
        molecular system. The grid is defined by a specified resolution and includes
        an additional border around the atomic region. The computed grid points are
        evaluated for property values to create the data cube.

        Parameters:
            resolution: float
                The spacing between grid points in Bohr.
            border: float
                The distance to extend the grid beyond the minimum and maximum atomic
                coordinates in all dimensions, in Bohr.

        Returns:
            CubeData
                An object containing details about the generated data cube, including
                atomic information, grid properties, and evaluated data.
        """
        atoms = self.atoms
        xyz = np.array([atom['xyz'] for atom in atoms])
        origin = np.min(xyz,axis=0)-border
        far_corner = np.max(xyz,axis=0)+border
        cells = np.diag([resolution,resolution,resolution])
        dimensions = [int((far_corner[i]-origin[i])/resolution+1) for i in range(3)]
        points3d = np.empty([dimensions[0]*dimensions[1]*dimensions[2],4],np.double)
        for i in range(dimensions[0]):
            for j in range(dimensions[1]):
                for k in range(dimensions[2]):
                    points3d[i*dimensions[1]*dimensions[2]+j*dimensions[2]+k,0] = origin[0]+i*resolution
                    points3d[i*dimensions[1]*dimensions[2]+j*dimensions[2]+k,1] = origin[1]+j*resolution
                    points3d[i*dimensions[1]*dimensions[2]+j*dimensions[2]+k,2] = origin[2]+k*resolution
                    points3d[i*dimensions[1]*dimensions[2]+j*dimensions[2]+k,3] = 1.0
        data = self.evaluate(points3d,values=True)

        return CubeData({
            'atoms':atoms,
            'natoms':len(atoms),
            'origin':origin,
            'cells':cells,
            'dimensions':dimensions,
            'orbitals':True,
            'data':np.reshape(data,dimensions),
            'title':[self.attribute('ID'),''],
            'orbital_identifiers':[self.attribute('ID')],
        })

def regular_grid(bounds,resolution):
    points = []
    for i in range(3):
        points.append(np.linspace(bounds[i][0],bounds[i][1],round(bounds[i][1]-bounds[i][0])/resolution)+1)
    points3d = np.ndarray((points[0].size,points[1].size,points[2].size,3),order='F')
    offset = 0
    for ix in range(points[0].size):
        for iy in range(points[1].size):
            for iz in range(points[2].size):
                points3d[offset,:] = [points[0][ix],points[1][iy],points[2][iz]]
                offset += 1
    return points, points3d

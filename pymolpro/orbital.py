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

    def grid(self,npt, method='erfinv'):
        return []

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


        self.linear_gaussians = []
        for i in range(3):
            self.linear_gaussians.append(self.construct_linear_gaussian(i))

        # self.linear_gaussian_areas = []
        # for i in range(3):
        #     self.linear_gaussian_areas.append(self.linear_gaussian_area(i))

        # self.linear_gaussian_region_deltas = []
        # for i in range(3):
        #     delta = self.auto_delta(i)
        #     self.linear_gaussian_region_deltas.append(delta)

    # def auto_delta(self, axis):
    #     area = self.linear_gaussian_area(axis)
    #     delta = math.sqrt(area)
    #     delta = int(delta)

    #     return delta

    def linear_gaussian_width(self, axis):
        second_moment = self.second_moments_matrix[axis][axis]
        width = 1 / (2 * (second_moment ** 2))
        # width = 1/(2*second_moment)

        return width

    # construct a mixture of Gaussians from a combination of linear Gaussian distributions
    def construct_linear_gaussian(self, axis):
        center = self.dipole[axis]
        width = self.linear_gaussian_width(axis)

        gaussian = lambda x: (2 * width / math.pi) * math.exp(-width * (x - center) ** 2)

        return gaussian

    def linear_gaussian_area(self, axis):
        width = self.linear_gaussian_width(axis)
        # area = math.sqrt(2*math.pi/width)
        # allow for negative values
        print("width: " + str(width))
        area = math.sqrt(2 * math.pi / width)
        return area

    def linear_gaussian_region_area(self, axis, delta):
        area = self.linear_gaussian_area(axis)
        area = area / delta
        return area

    def linear_gaussian_region_widths(self, axis, delta):
        width = self.linear_gaussian_width(axis)
        region_widths = []
        for i in range(delta):
            region_width = math.sqrt(2 * math.pi / (width * (delta - i)))
            region_widths.append(region_width)

        return region_widths

    def linear_gaussian_region_centers(self, axis, delta):
        # center = self.dipole[axis]
        center = 0

        centers = []
        # center_count = center
        center_count = 0
        for i in range(delta):
            width = self.linear_gaussian_region_widths(axis, delta)[i]
            center_count = center_count + width
            centers.append(center_count)

        # reflect the centers around the center
        reflected_centers = [i for i in centers[::-1]]
        reflected_points = []
        for i in range(len(reflected_centers)):
            # reflected_points.append(center - reflected_centers[i])
            reflected_points.append(-reflected_centers[i])

        print("reflected_points len: " + str(len(reflected_points)))
        print("centers len: " + str(len(centers)))
        print("center len: " + str(len([center])))

        centers = reflected_points + [center] + centers

        print("centers len: " + str(len(centers)))

        return centers

    # def linear_gaussian_region_centers(self, axis, delta):
    #     width = self.linear_gaussian_width(axis)
    #     half_width = width/2
    #     low_bound = -half_width
    #     high_bound = half_width

    #     # use the inverse error function to calculate points
    #     print("inputs: ", low_bound, high_bound, delta)
    #     linspace = np.linspace(low_bound, high_bound, num=delta)
    #     points = sp.special.erfinv(linspace)

    def gaussian_region_centers(self, delta):
        half_delta = int(delta / 2)
        # half_delta = delta
        region_centers = []
        # for i in range(delta):
        #     for j in range(delta):
        #         for k in range(delta):
        #             region_centers.append([self.linear_gaussian_region_centers(0, half_delta)[i], self.linear_gaussian_region_centers(1, half_delta)[j], self.linear_gaussian_region_centers(2, half_delta)[k]])

        x_centers = self.linear_gaussian_region_centers(0, half_delta)
        y_centers = self.linear_gaussian_region_centers(1, half_delta)
        z_centers = self.linear_gaussian_region_centers(2, half_delta)

        delta += 1
        for i in range(delta):
            for j in range(delta):
                for k in range(delta):
                    region_centers.append([x_centers[i], y_centers[j], z_centers[k]])

        # for i in range(3):
        #     region_centers.append(self.linear_gaussian_region_centers(i, half_delta))

        # for i in range(3):
        #     region_centers[i] = region_centers[i] + [0] + region_centers[i][::-1]

        return region_centers

    def get_gaussian_region_centers(self, delta):
        try:
            error_message = "delta must be an integer or a list of integers"
            if type(delta) == int:
                return self.gaussian_region_centers(delta)
            elif type(delta) == list:
                region_centers = []
                for i in range(len(delta)):
                    if type(i) == int:
                        region_centers.append(self.gaussian_region_centers(delta[i]))
                    else:
                        raise Exception(error_message)
                return region_centers
        except:
            raise Exception(error_message)

    # we must rotate the translate the orbital back to its original position
    def get_points(self, delta):
        points = self.get_gaussian_region_centers(delta)

        # we need to rotate the points back to the original coordinate system
        rotated_points = []
        rotation_matrix_inverse = self.inverse_rotation_matrix
        for i in range(len(points)):
            point = points[i]
            point = np.array(point)
            point = np.matmul(rotation_matrix_inverse, point)
            point = point.tolist()
            rotated_points.append(point)

        # return rotated_points

        # rotated_points = points

        # rotated_points = points

        transformed_points = []
        vector = self.dipole

        # # print("vector: ", vector)
        for i in range(len(rotated_points)):
            point = rotated_points[i]
            point = np.array(point)
            point = point + vector
            point = point.tolist()
            transformed_points.append(point)

        return transformed_points

    def print_linear_gaussian_data(self, axis):
        print("center: " + str(self.dipole[axis]))
        print("width: " + str(self.linear_gaussian_width(axis)))
        print("area: " + str(self.linear_gaussian_area(axis)))

    def print_gaussian_data_xyz(self):
        for i in range(3):
            self.print_linear_gaussian_data(i)

    # for each axis print the linear_gaussian_width
    def print_linear_gaussian_widths(self):
        for i in range(3):
            print("linear_gaussian_width for axis " + str(i) + " is " + str(self.linear_gaussian_width(i)))

    # for each axis print the linear_gaussian_area
    def print_linear_gaussian_areas(self):
        for i in range(3):
            print("linear_gaussian_areas " + str(i) + " is " + str(self.linear_gaussian_area(i)))

    def print_linear_gaussian_region_widths(self, delta):
        for i in range(3):
            print(
                "linear_gaussian_region_widths " + str(i) + " is " + str(self.linear_gaussian_region_widths(i, delta)))

    def print_linear_gaussian_region_areas(self, delta):
        for i in range(3):
            print("linear_gaussian_region_areas " + str(i) + " is " + str(self.linear_gaussian_region_area(i, delta)))

    def print_linear_gaussian_region_centers(self, delta):
        for i in range(3):
            print("linear_gaussian_region_centers " + str(i) + " is " + str(
                self.linear_gaussian_region_centers(i, delta)))

    def print_gaussian_region_centers(self, delta):
        print("gaussian_region_centers is " + str(self.get_gaussian_region_centers(delta)))

    # for jupyter notebook
    import math
    from matplotlib import pyplot as plt

    def plot_coordinates_3d(coordinates):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for coordinate in coordinates:
            x = coordinate[0]
            y = coordinate[1]
            z = coordinate[2]
            ax.scatter(x, y, z)

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        plt.show()

    # plot charge density as [[x1, y1, z1, c1], [x2, y2, z2, c2], ...]
    def plot_charge_density_3d(charge_density):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for point in charge_density:
            x = point[0]
            y = point[1]
            z = point[2]
            charge = point[3]
            ax.scatter(x, y, z, s=charge * 100)

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        plt.show()

    def plot_3d_matrix_shape(matrix):
        # matrix is a np.array like [[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]]
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for i in range(3):
            for j in range(3):
                x = i
                y = j
                z = matrix[i][j]
                ax.scatter(x, y, z)

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        plt.show()

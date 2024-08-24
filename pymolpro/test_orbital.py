import unittest
import pymolpro
import math
import numpy as np
from lxml import etree


class TestOrbital(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = etree.Element("molpro")
        orbital_node = etree.SubElement(cls.root, "orbital", occupation='2.0', ID='1.1',
                                        moments='0.0 0.0 0.0 1.0 2.0 3.0 1.0 0.7 -0.4 0.0')
        orbital_node.text = '1.0 2.0'
        cls.orbital = pymolpro.Orbital(orbital_node)
        orbital_node = etree.SubElement(cls.root, "orbital", occupation='2.0', ID='1.1',
                                        moments='0.0 0.0 0.0 1.0 1.0 1.0 0.0 0.0 0.0 0.0')
        orbital_node.text = '1.0 2.0'
        cls.aligned_orbital = pymolpro.Orbital(orbital_node)

    def test_construction(self):
        self.assertEqual(self.orbital.ID, '1.1')  # add assertion here
        self.assertEqual(self.orbital.occupation, 2.0)  # add assertion here
        self.assertEqual(self.orbital.coefficients, [1.0, 2.0])

    def test_axes(self):
        self.assertEqual(self.orbital.ID, '1.1')
        reconstructed_second_moments = np.matmul(self.orbital.axes,
                                                 np.matmul(np.diagflat(self.orbital.second_moment_eigenvalues),
                                                           self.orbital.axes.transpose()))
        np.testing.assert_almost_equal(reconstructed_second_moments, self.orbital.local_second_moments)

    def test_large_grid(self):
        for grid in ['Gauss-Hermite', 'Lebedev-Mura', 'Lebedev-Gauss-Laguerre']:
            points = self.aligned_orbital.grid(47, grid,
                                               scale=2.0 if 'Mura' in grid else 0.5 if 'Hermite' in grid else 1.0)
            integral = 0.0
            for i in range(len(points)):
                # integral += weights[i] * math.exp(-np.linalg.norm(points[i, :])) / 4 / math.pi
                integral += points[i, 3] * math.exp(-pow(np.linalg.norm(points[i, :3]), 2)) / pow(math.pi, 1.5)
            self.assertAlmostEqual(integral, 1.0, 11, msg=grid)

    def test_gauss_hermite(self):
        from scipy.special import factorial2
        size = 15
        points = self.aligned_orbital.grid(size, 'Gauss-Hermite', scale=1.0)
        integrals = []
        for power in range(0, size // 2):
            f2 = factorial2(2 * power - 1) if power > 0 else 1  # circumvent bug in scipy 1.11.1
            integral002 = 0.0
            integral022 = 0.0
            for p in points:
                x2 = pow(p[0], 2)
                y2 = pow(p[1], 2)
                z2 = pow(p[2], 2)
                r2 = x2 + y2 + z2
                integral002 += p[3] * pow(z2, power) * math.exp(-r2 / 2)
                integral022 += p[3] * pow(y2 * z2, power) * math.exp(-r2 / 2)
            integrals.append(integral002 / f2 / pow(2 * math.pi, 1.5))
            integrals.append(integral022 / pow(f2, 2) / pow(2 * math.pi, 1.5))
        for integral in integrals:
            self.assertAlmostEqual(integral, 1, 13, msg=integrals)

    def test_mura(self):
        integrals = []
        for size in range(1, 10):
            for m in range(1, 5):
                points = self.aligned_orbital.grid(size, 'Lebedev-Mura', scale=1.0, grid_parameters=[m, 1])
                alpha = 0.7475238520470884
                alpha = 1
                integral = 0.0
                integralx = 0.0
                for p in points:
                    x2 = pow(p[0], 2)
                    y2 = pow(p[1], 2)
                    z2 = pow(p[2], 2)
                    r2 = x2 + y2 + z2
                    r = math.sqrt(r2)
                    x = pow(1 - math.exp(-r / alpha), 1.0 / float(m))
                    integrand = math.exp(-r) / r2 / 4 / math.pi / (
                            pow(x, m - 1) * alpha * m)  # should integrate to 1 exactly on this grid
                    integral += p[3] * integrand
                    integralx += 2 * x * p[3] * integrand
                integrals.append(integral)
                integrals.append(integralx)
        for integral in integrals:
            self.assertAlmostEqual(integral, 1, 12, msg=integrals)

    def test_mura_exponential(self):
        beta = 2 * math.sqrt(3.0)
        size = 20
        alpha = 1 / beta
        alpha = 1
        scales = []
        errors = []
        smerrors = []
        ldaerrors = []
        coulomberrors = []
        for scalelog in range(-10, 10):
            scale = math.exp(scalelog / 100.0)
            alpha = scale / beta
            for m in range(2, 5):
                points = self.aligned_orbital.grid(size, 'Lebedev-Mura', scale=alpha, grid_parameters=[m])
                # print("alpha=",alpha)
                # print(points)
                norm = 0
                sm = 0
                lda = 0
                coulomb = 0
                for p in points:
                    r = np.linalg.norm(p[:3])
                    norm += p[3] * (pow(beta, 3) / 8 / math.pi) * math.exp(-beta * r)
                    lda += p[3] * pow((pow(beta, 3) / 8 / math.pi) * math.exp(-beta * r), 4.0 / 3.0)
                    coulomb += p[3] / r * (pow(beta, 3) / 8 / math.pi) * math.exp(-beta * r)
                    sm += p[3] * r * r * (pow(beta, 3) / 8 / math.pi) * math.exp(-beta * r)
                # print(norm)
                scales.append(scale)
                errors.append(abs(norm - 1))
                ldaerrors.append(abs(lda - 27 * beta / 128.0 / pow(math.pi, 1.0 / 3.0)))
                coulomberrors.append(abs(coulomb - beta / 2))
                smerrors.append(abs(sm - 1))
                self.assertLessEqual(smerrors[-1], 1e-4)
                self.assertLessEqual(errors[-1], 2e-6)
                self.assertLessEqual(ldaerrors[-1], 1e-7)
                self.assertLessEqual(coulomberrors[-1], 1e-4)
        if False:
            import matplotlib.pyplot as plt
            plt.loglog(scales, errors)
            plt.loglog(scales, smerrors)
            plt.loglog(scales, ldaerrors)
            plt.loglog(scales, coulomberrors)
            plt.legend(['norm', 'second moment', 'lda', 'coulomb'])
            plt.show()


if __name__ == '__main__':
    unittest.main()

import unittest
import pymolpro
import math
import numpy as np
from lxml import etree


class TestOrbital(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = etree.Element("molpro")
        orbital_node = etree.SubElement(root, "orbital", occupation='2.0', ID='1.1',
                                        moments='0.0 0.0 0.0 1.0 1.0 1.0 0.0 0.0 0.0 0.0')
        orbital_node.text = '1.0 2.0'
        cls.orbital = pymolpro.Orbital(orbital_node)

    def test_construction(self):
        root = etree.Element("molpro")
        orbital_node = etree.SubElement(root, "orbital", occupation='2.0', ID='1.1',
                                        moments='0.0 0.0 0.0 1.0 1.0 1.0 0.0 0.0 0.0 0.0')
        orbital_node.text = '1.0 2.0'
        orbital = pymolpro.Orbital(orbital_node)
        self.assertEqual(orbital.ID, '1.1')  # add assertion here
        self.assertEqual(orbital.occupation, 2.0)  # add assertion here
        self.assertEqual(orbital.coefficients, [1.0, 2.0])

    def test_large_grid(self):
        self.assertEqual(self.orbital.ID, '1.1')
        for grid in ['Gauss-Hermite', 'Lebedev-Mura', 'Lebedev-Gauss-Laguerre']:
            points = self.orbital.grid(55, grid, scale=2.0 if 'Mura' in grid else 1.0)
            integral = 0.0
            for i in range(len(points)):
                # integral += weights[i] * math.exp(-np.linalg.norm(points[i, :])) / 4 / math.pi
                integral += points[i, 3] * math.exp(-pow(np.linalg.norm(points[i, :3]), 2)) / pow(math.pi, 1.5)
            self.assertAlmostEqual(integral, 1.0, 12, msg=grid)


if __name__ == '__main__':
    unittest.main()

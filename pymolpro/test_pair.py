import unittest

import lxml.etree

import pymolpro
import numpy as np


class TestPair(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = lxml.etree.fromstring('''<?xml version="1.0"?>
<molpro xmlns="http://www.molpro.net/schema/molpro-output"
  xmlns:xsd="http://www.w3.org/1999/XMLSchema"
  xmlns:cml="http://www.xml-cml.org/schema"
  xmlns:stm="http://www.xml-cml.org/schema"
  xmlns:xhtml="http://www.w3.org/1999/xhtml">
  <jobstep>
  <pair orbital1="1.1" orbital2="-2.1" energy="3.2345"/>
  <pair orbital1="1.1" orbital2="2.1" energy="4.2345"/>
  </jobstep>
  <molecule>
  <orbitals>
  <orbital occupation="2.0" ID="1.1" moments="0.0 0.0 -1.0 2.0 2.0 2.0 0.0 0.0 0.0 0.0">-1.0 -2.0</orbital>
  <orbital occupation="2.0" ID="2.1" moments="0.0 0.0 1.0 2.0 2.0 2.0 0.0 0.0 0.0 0.0">-3.0 -4.0</orbital>
  </orbitals>
  </molecule>
  <molecule>
  <orbitals>
  <orbital occupation="2.0" ID="1.1" moments="0.0 0.0 -1.0 2.0 2.0 2.0 0.0 0.0 0.0 0.0">1.0 2.0</orbital>
  <orbital occupation="2.0" ID="2.1" moments="0.0 0.0 1.0 2.0 2.0 2.0 0.0 0.0 0.0 0.0">3.0 4.0</orbital>
  </orbitals>
  </molecule>
  <jobstep>
  <pair orbital1="1.1" orbital2="-2.1" energy="1.2345"/>
  <pair orbital1="1.1" orbital2="2.1" energy="0.2345"/>
  </jobstep>
  <molecule>
  <orbitals>
  <orbital occupation="2.0" ID="1.1" moments="0.0 0.0 -1.0 2.0 2.0 2.0 0.0 0.0 0.0 0.0">-1.0 -2.0</orbital>
  <orbital occupation="2.0" ID="2.1" moments="0.0 0.0 1.0 2.0 2.0 2.0 0.0 0.0 0.0 0.0">-3.0 -4.0</orbital>
  </orbitals>
  </molecule>
  <molecule>
  <orbitals>
  <orbital occupation="2.0" ID="1.1" moments="0.0 0.0 -1.0 2.0 2.0 2.0 0.0 0.0 0.0 0.0">-11.0 -12.0</orbital>
  <orbital occupation="2.0" ID="2.1" moments="0.0 0.0 1.0 2.0 2.0 2.0 0.0 0.0 0.0 0.0">-13.0 -14.0</orbital>
  </orbitals>
  </molecule>
  </molpro>''')
        # print(etree.tostring(root, pretty_print=True))
        cls.pair = pymolpro.Pair(pymolpro.xpath(root, "//pair")[-2])
        cls.root = root

    def test_construction(self):
        pair = pymolpro.Pair(pymolpro.xpath(self.root, "//pair")[-2])
        self.assertEqual(len(pair.orbitals), 2)  # add assertion here
        self.assertEqual(pair.orbitals[0].ID, '1.1')  # add assertion here
        self.assertEqual(pair.orbitals[1].ID, '2.1')  # add assertion here
        self.assertEqual(pair.spins, [1, -1])
        self.assertAlmostEqual(float(pair.node.get('energy')), 1.2345)
        self.assertAlmostEqual(pair.energy, 1.2345)
        self.assertEqual(pair.orbitals[0].node.text, "1.0 2.0")
        self.assertEqual(len(pair), 2)

    def test_construction_forward_orbitals(self):
        pair = pymolpro.Pair(pymolpro.xpath(self.root, "//pair")[1])
        self.assertEqual(len(pair.orbitals), 2)  # add assertion here
        self.assertEqual(pair.orbitals[0].ID, '1.1')  # add assertion here
        self.assertEqual(pair.orbitals[1].ID, '2.1')  # add assertion here
        self.assertEqual(pair.spins, [1, 1])
        self.assertAlmostEqual(float(pair.node.get('energy')), 4.2345)
        self.assertAlmostEqual(pair.energy, 4.2345)
        self.assertEqual(pair.orbitals[0].node.text, "-1.0 -2.0")
        self.assertEqual(len(pair), 2)

    def test_axes(self):
        pair = pymolpro.Pair(pymolpro.xpath(self.root, "//pair")[1])
        self.assertTrue(np.allclose(np.matmul(pair.axes, pair.axes.transpose()), np.identity(3)))
        self.assertTrue(np.allclose(np.matmul(pair.axes.transpose(), pair.axes), np.identity(3)))
        self.assertTrue(np.allclose(pair.axes, np.identity(3)))


if __name__ == '__main__':
    unittest.main()

import unittest

import lxml.etree

import pymolpro
import math
import numpy as np
from lxml import etree


class TestPair(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = lxml.etree.fromstring('''<?xml version="1.0"?>
<molpro xmlns="http://www.molpro.net/schema/molpro-output"
  xmlns:xsd="http://www.w3.org/1999/XMLSchema"
  xmlns:cml="http://www.xml-cml.org/schema"
  xmlns:stm="http://www.xml-cml.org/schema"
  xmlns:xhtml="http://www.w3.org/1999/xhtml">
  <orbitals>
  <orbital occupation="2.0" ID="1.1" moments="0.0 0.0 -1.0 2.0 2.0 2.0 0.0 0.0 0.0 0.0">1.0 2.0 </orbital>
  <orbital occupation="2.0" ID="2.1" moments="0.0 0.0 1.0 2.0 2.0 2.0 0.0 0.0 0.0 0.0">3.0 4.0 </orbital>
  </orbitals>
  <pair orbital1="1.1" orbital2="-2.1" energy="1.2345"/>
  <pair orbital1="1.1" orbital2="2.1" energy="0.2345"/>
  </molpro>''')
        print(etree.tostring(root, pretty_print=True))
        cls.pair = pymolpro.Pair(pymolpro.xpath(root, "//pair")[-2])
        cls.root = root

    def test_construction(self):
        pair = self.pair
        self.assertEqual(len(pair.orbitals), 2)  # add assertion here
        self.assertEqual(pair.orbitals[0].ID, '1.1')  # add assertion here
        self.assertEqual(pair.orbitals[1].ID, '2.1')  # add assertion here
        self.assertEqual(pair.spins,[1,-1])
        self.assertAlmostEqual(float(pair.node.get('energy')),1.2345)
        self.assertAlmostEqual(pair.energy,1.2345)
        print(self.pair.data())


if __name__ == '__main__':
    unittest.main()

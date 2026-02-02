import subprocess
import unittest

import pytest

import pymolpro
import os
import shutil


class TestRegistry(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_procedures_registry(self):
        proc_reg = pymolpro.procedures_registry()
        if proc_reg:
            assert proc_reg['PNO-UCCSD']['gradient'] == -1
            assert 'pairopt' in proc_reg['PNO-UCCSD']['options']

    def test_basis_registry(self):
        basis_reg = pymolpro.basis_registry()
        assert None not in basis_reg.keys()
        if basis_reg:
            assert basis_reg['cc-pV(D+d)Z']['quality'] == 'DZ'

    def test_molpro_root(self):
        # print('registry',pymolpro.run_local_molpro(['--registry']).stdout)
        if pymolpro.local_molpro_root():
            assert os.path.exists(pymolpro.local_molpro_root() / 'lib' / 'defbas')

    def test_registry(self):
        # print('registry', pymolpro.registry(), '\n\n')
        # print('GMB registry', pymolpro.registry('GMB'), '\n\n')
        # print('EOM registry', pymolpro.registry('EOM'), '\n\n')
        if pymolpro.registry():
            assert 'GMB' in pymolpro.registry()
            assert 'STATES' in pymolpro.registry('GMB')
            assert pymolpro.registry('GMB')['STATES']['set'] == 'GMB'


if __name__ == '__main__':
    unittest.main()

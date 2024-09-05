import unittest

import pytest

import pymolpro
import os


class TestProject(unittest.TestCase):
    def setUp(self):
        self.project = pymolpro.Project("TestProject", location=os.path.dirname(os.path.abspath(__file__)))

    def test_local_molpro(self):
        print('local molpro root',self.project.local_molpro_root)

    def test_pairs_discovery(self):
        pair_energies = 0.0
        for tuple in self.project.singles() + self.project.pairs():
            pair_energies += tuple.energy
            # print(tuple.energy, tuple.spins, [orbital.ID for orbital in tuple.orbitals])
        # print(self.project.energies())
        correlation_energy = self.project.energies()[-1] - self.project.energies()[0]
        # print(correlation_energy, pair_energies)
        self.assertAlmostEqual(correlation_energy, pair_energies)
        # print(self.project.properties('correlation energy')) #TODO fix xpath search where attribute has embedded spaces

    def test_copy(self):
        newproject = self.project.copy('copied', location=os.path.dirname(os.path.abspath(__file__)))
        newproject.erase()

    def test_corrupt_project(self):
        import os, shutil
        import pysjef
        pname = 'check_empty.molpro'
        shutil.rmtree(pname,ignore_errors=True)
        os.mkdir(pname)
        with open(pname + '/Info.plist', 'w') as f:
            f.write('')
        with pytest.raises(FileNotFoundError, match=r'Cannot open project') as e:
            p = pysjef.Project(pname)
        assert e.type is FileNotFoundError

    def test_procedures_registry(self):
        proc_reg = self.project.procedures_registry()
        if proc_reg:
            assert proc_reg['PNO-UCCSD']['gradient'] == -1
            assert 'pairopt' in proc_reg['PNO-UCCSD']['options']
    def test_basis_registry(self):
        basis_reg = self.project.basis_registry()
        assert None not in basis_reg.keys()
        if basis_reg:
            assert basis_reg['cc-pV(D+d)Z']['quality'] == 'DZ'

    def test_molpro_root(self):
        # print('registry',self.project.run_local_molpro(['--registry']).stdout)
        if self.project.local_molpro_root:
            assert os.path.exists(self.project.local_molpro_root / 'lib' / 'defbas')

    def test_registry(self):
        if self.project.registry():
            assert 'GMB' in self.project.registry()
            assert 'STATES' in self.project.registry('GMB')
            assert self.project.registry('GMB')['STATES']['set'] == 'GMB'

    def test_orbitals(self):
        assert len(self.project.orbitals()) == 2

if __name__ == '__main__':
    unittest.main()

import unittest
import pymolpro
import os


class TestProject(unittest.TestCase):
    def setUp(self):
        self.project = pymolpro.Project("TestProject", location=os.path.dirname(os.path.abspath(__file__)))

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

    # def test_corrupt_project(self):
    #     import os, shutil
    #     import pysjef
    #     pname = 'check_empty.molpro'
    #     shutil.rmtree(pname,ignore_errors=True)
    #     os.mkdir(pname)
    #     open(pname + '/Info.plist', 'w').write('')
    #     try:
    #         p = pysjef.Project(pname)
    #     except:
    #         print("exception caught")

    def test_procedures_registry(self):
        proc_reg = self.project.procedures_registry()
        if proc_reg:
            assert proc_reg['PNO-UCCSD']['gradient'] == -1
            assert 'pno' in proc_reg['PNO-UCCSD']['options']
    def test_basis_registry(self):
        basis_reg = self.project.basis_registry()
        if basis_reg:
            assert basis_reg['cc-pV(D+d)Z']['quality'] == 'DZ'

    def test_molpro_root(self):
        molpro_root = self.project.local_molpro_root()
        if molpro_root:
            assert os.path.exists(molpro_root/'lib'/'defbas')

    def test_registry(self):
        assert 'GMB' in self.project.registry()
        assert 'STATES' in self.project.registry('GMB')
        assert self.project.registry('GMB')['STATES']['set'] == 'GMB'

if __name__ == '__main__':
    unittest.main()

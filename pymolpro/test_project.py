import subprocess
import unittest

import pytest

import pymolpro
import os
import shutil


class TestProject(unittest.TestCase):
    def setUp(self):
        self.project = pymolpro.Project("TestProject", location=os.path.dirname(os.path.abspath(__file__)))
        self.projects = []

    def tearDown(self):
        for project in self.projects:
            project.erase()
        pass

    def new_project(self, *args, **kwargs):
        self.projects.append(pymolpro.Project(*args, **kwargs))
        return self.projects[-1]

    # def test_project_from_files(self):
    #     shutil.rmtree('test_project_from_files.molpro', ignore_errors=True)
    #     self.projects.append(pymolpro.Project("test_project_from_files", files=['/tmp/molpro.xml','/tmp/whatsit.xyz','/tmp/whatsit.molden']))
    #     subprocess.run(['ls','-lR','test_project_from_files.molpro'])

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
        shutil.rmtree(os.path.dirname(os.path.abspath(__file__))+'/copied.molpro',ignore_errors=True)
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

    def test_run_needed(self):
        print(self.project.run_needed())
        print(pymolpro.__version__)
        # %%
        dir = 'run_needed'
        shutil.rmtree(dir, ignore_errors=True)
        os.mkdir(dir)
        # %%
        p = pymolpro.Project('one', location=dir)
        p.write_input('geometry={He};rhf')
        print('Should be True:', p.run_needed())
        # %%
        p.run(wait=True)
        print('Should be False:', p.run_needed())
        # %%
        p = pymolpro.Project('one', location=dir)
        p.write_input('geometry={He};rhf;mp3')
        print('Should be True:', p.run_needed())
        # %%
        p2 = pymolpro.Project('two', location=dir)
        p2.write_input('geometry={He};rhf;mp3')
        print('Should be True:', p2.run_needed())
        # %%
        shutil.rmtree(dir, ignore_errors=True)

    def test_ansatz(self):
        for ansatz in [
            'B3LYP/cc-pVTZ',
            'CCSD(T)-F12A/cc-pVTZ//B3LYP/cc-pVDZ',
            'DF-MP2/aug-cc-pvdz',
            'UmP2/aug-cc-pvdz',
            'HF/cc-pvdz',
        ]:
            p1 = self.new_project('test_ansatz' + ansatz.replace('/', '_'), ansatz=ansatz, geometry='He')
            self.assertEqual(ansatz, p1.ansatz, open(p1.filename('inp'), 'r').read())

    def test_force_constants(self):
        import numpy as np
        for r in [2.0, 4.0]:
            p = self.new_project('test_hessian')
            dr=.001
            p.write_input(f'basis,svp;bohr;geometry={{F,0,0,0;H,0,0,{r}}};rhf;frequencies,analytic')
            if p.local_molpro_root is not None:
                p.run(wait=True)
                # print('gradient', p.gradient())
                hessian_0 = p.vibrations['force_constants']
                try:
                    # print('hessian', hessian_0 )
                    pass
                except:
                    print(p.out)
                pp = self.new_project('test_hessianp')
                pp.write_input(f'basis,svp;bohr;geometry={{F,0,0,0;H,0,0,{r+dr}}};rhf;forces')
                pm = self.new_project('test_hessianm')
                pm.write_input(f'basis,svp;bohr;geometry={{F,0,0,0;H,0,0,{r-dr}}};rhf;forces')
                pp.run()
                pm.run()
                pp.wait()
                pm.wait()
                numhess = (pp.gradient()-pm.gradient())/(2*dr)
                np.testing.assert_almost_equal(numhess, hessian_0[5,:],decimal=6)
                np.testing.assert_almost_equal(-numhess, hessian_0[2,:],decimal=6)
                # print('numhess', numhess)
                # print(p.gradient())
                # print(pp.gradient())
                # print(pm.gradient())
                # print((pp.gradient()+pm.gradient())/2)

if __name__ == '__main__':
    unittest.main()

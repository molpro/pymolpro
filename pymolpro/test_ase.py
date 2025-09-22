import unittest
import pymolpro
import ase.build
import ase.optimize


class TestASE(unittest.TestCase):

    def test_opt(self):
        if pymolpro.Project('test').local_molpro_root:
            atoms = ase.build.molecule('H2O')

            atoms.calc = pymolpro.ASEMolpro(method='df-hf', basis='cc-pVDZ')

            with ase.optimize.BFGS(atoms) as opt:
                assert opt.run(fmax=0.0001)

            atoms.calc.clean()

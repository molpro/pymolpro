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

            E = atoms.get_potential_energy()
            print(E)
            atoms.set_distance(0,1,1.1)
            E = atoms.get_potential_energy()
            print(E)

            atoms.calc.clean()

# from pymolpro import ASEMolpro
# from ase.optimize import BFGS
#
# h2 = ase.build.molecule('H2')
# h2.calc = ASEMolpro(ansatz='df-lmp2/aug-cc-pVTZ')
#
# with BFGS(h2) as opt:
#     opt.run(fmax=0.0001)
#
# E = h2.get_potential_energy()      # compute potential energy
#
# h2.set_distance(0, 1, 0.7)           # set H–H distance
# E = h2.get_potential_energy()      # compute potential energy
# print(f"H–H distance: {0.7} Å, Energy: {E:.6f} eV")

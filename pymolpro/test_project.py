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


if __name__ == '__main__':
    unittest.main()

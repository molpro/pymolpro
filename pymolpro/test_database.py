import unittest

import pymolpro
from pymolpro.database import Database
import os
import shutil


class TestDatabase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    def test_empty_construction(self):
        db = Database()
        self.assertEqual(len(db.molecules), 0)  # add assertion here
        self.assertEqual(len(db.reactions), 0)  # add assertion here
        db.add_molecule('H2', 'H 0 0 0\nH 0 0 .7')
        self.assertEqual(len(db.molecules), 1)  # add assertion here

    def test_filled_construction(self):
        db = Database(molecules={
            'H2': 'H 0 0 0\nH 0 0 .8',
            'F2': 'F 0 0 0\nF 0 0 1.43',
            'HF': 'H 0 0 0\nF 0 0 1',
        })
        self.assertEqual(len(db.molecules), 3)  # add assertion here
        self.assertEqual(len(db.reactions), 0)  # add assertion here

    def test_json(self):
        db = Database(
            molecules={
                'H2': 'H 0 0 0\nH 0 0 .8',
                'F2': 'F 0 0 0\nF 0 0 1.43',
                'HF': 'H 0 0 0\nF 0 0 1',
                'HFHF': """
H          0.0000000000        0.0000000000       -0.1672052678
F          0.0000000000        0.0000000000        0.7584550750
H          0.0000000000        0.0000000000        2.7403780100
F          0.0000000000        0.0000000000        3.6683721829"""
            },
        )
        db.add_reaction('H2+F2',
                        {'HF': 2,
                         'H2': -1,
                         'F2': -1,
                         }, description='Hydrogenation of difluorine')
        db.add_reaction('HFHF',
                        {'HF': 2,
                         'HFHF': -1,
                         }, description='Binding of hydrogen fluoride dimer')
        db.add_subset('non-covalent', ['HFHF'])
        db.dump("sample.json")
        db2 = Database()
        db2.load("sample.json")
        os.remove("sample.json")
        self.assertEqual(db.molecules, db2.molecules)
        self.assertEqual(db.reactions, db2.reactions)

    def test_library(self):
        db = pymolpro.database.load('sample')
        self.assertEqual(len(db.molecules), 4)
        self.assertEqual(len(db.reactions), 2)
        raw_path = os.path.realpath(os.path.join(__file__, '..', '..', 'pymolpro', 'share', 'database', 'sample'))
        db_raw_path = pymolpro.database.load(raw_path)
        self.assertEqual(db, db_raw_path)
        raw_path = os.path.realpath(os.path.join(__file__, '..', '..', 'pymolpro', 'share', 'database', 'sample.json'))
        db_raw_path = pymolpro.database.load(raw_path)
        self.assertEqual(db, db_raw_path)
        db_string = pymolpro.database.load(' '.join(open(raw_path, 'r').readlines()))
        self.assertEqual(db, db_string)

    def test_run_database(self):
        if shutil.which('molpro'):
            db = pymolpro.database.load('sample')
            results = pymolpro.database.run(db, method='df-lmp2', basis='aug-cc-pVDZ', func="energy")
            results2 = pymolpro.database.run(db, method='df-lmp2', basis='aug-cc-pVDZ', func="energy")
            # del results['projects']
            # del results2['projects']
            self.assertEqual(results.molecule_energies, results2.molecule_energies)
            self.assertEqual(results, results2)
            results2.molecule_energies['HF'] = 12345.0
            self.assertNotEqual(results.molecule_energies, results2.molecule_energies)
            self.assertNotEqual(results, results2)
            pymolpro.database.run(db, method='df-lmp2', basis='aug-cc-pVDZ', func="energy", clean=True)
            # print(results)

    def test_compare_database_runs(self):
        if shutil.which('molpro'):
            db = pymolpro.database.load('sample')
            results = [
                pymolpro.database.run(db, method='df-lmp2', basis='aug-cc-pVDZ', func="energy"),
                pymolpro.database.run(db, method='df-lmp2', basis='aug-cc-pVTZ', func="energy"),
                pymolpro.database.run(db, method='df-lmp2', basis='aug-cc-pVQZ', func="energy"),
            ]
            outputs = pymolpro.database.analyse(results, results[-1])
            self.assertEqual(len(outputs['reaction statistics'].index), 4)
            self.assertEqual(len(outputs['molecule statistics'].index), 4)
            self.assertEqual(len(outputs['reaction energies'].index), len(db))
            for title, output in outputs.items():
                self.assertEqual(len(output.columns), 3)
                # print(title)
                # print(output)
            db.reaction_energies = results[-1].reaction_energies
            db.molecule_energies = {}
            self.assertEqual(len(pymolpro.database.analyse(results, db)), 4)

    def test_fail(self):
        if shutil.which('molpro'):
            db = pymolpro.database.load('sample')

            result = pymolpro.database.run(db, method='hf', basis='minao')
            self.assertFalse(result.failed)
            shutil.rmtree(result.project_directory)

            result = pymolpro.database.run(db, method='hf', preamble='memory,1,k', basis='minao')
            self.assertTrue(result.failed)
            self.assertEqual(len(result.failed), len(db.molecules))
            self.assertEqual(result.failed['HF'].status, 'failed')
            shutil.rmtree(result.project_directory)

            result = pymolpro.database.run(db, method='ccsd', preamble='charge=-10', basis='cc-pvdz')
            self.assertTrue(result.failed)
            self.assertEqual(len(result.failed), 1)
            self.assertEqual(result.failed['HF'].status, 'failed')
            shutil.rmtree(result.project_directory)

            result = pymolpro.database.run(db, method='bad-method', basis='minao')
            self.assertEqual(len(result.failed), len(db.molecules))
            self.assertEqual(result.failed['HF'].status, 'failed')
            shutil.rmtree(result.project_directory)

    def test_subset(self):
        db = pymolpro.database.load('sample')
        self.assertEqual(len(db.subset('non-covalent')), 1)
        self.assertEqual(len(db.subset('non-covalent').molecules), 2)
        # print(db.subset('non-covalent'))

        subset = db.subset(['HFHF'])
        # print(subset)
        self.assertEqual(len(subset), 1)
        self.assertEqual(len(subset.molecules), 2)

        db.add_subset('vdW', 'HFHF')
        self.assertEqual(len(db.subset('vdW')), 1)
        self.assertEqual(len(db.subset('vdW').molecules), 2)
        # print(db.subset('vdW'))

        if shutil.which('molpro'):
            result = pymolpro.database.run(subset, clean=True)
            self.assertEqual(len(result.molecule_energies), 2)
            self.assertEqual(len(result.reaction_energies), 1)

    def test_extrapolate(self):
        if shutil.which('molpro'):
            db = pymolpro.database.load('sample').subset('non-covalent')
            hfresults = []
            results = []
            for basis in ['cc-pvdz', 'cc-pvtz']:
                hfresults.append(pymolpro.database.run(db, 'hf', basis))
                results.append(pymolpro.database.run(db, 'mp2', basis))
            results += pymolpro.database.basis_extrapolate(results, hfresults, [2, 3])
            results += pymolpro.database.basis_extrapolate(results, hfresults)
            # print(pymolpro.database.analyse(results))
            self.assertEqual(len(results), 4)

            hfresults = []
            results = []
            for basis in ['cc-pvdz', 'cc-pvtz', 'cc-pvqz']:
                hfresults.append(pymolpro.database.run(db, 'hf', basis))
                results.append(pymolpro.database.run(db, 'mp2', basis))
            results += pymolpro.database.basis_extrapolate(results, hfresults)
            # print(pymolpro.database.analyse(results))
            self.assertEqual(len(results), 5)

    def test_units(self):
        from pymolpro.database import units
        self.assertEqual(units[None], 1.0)
        self.assertEqual(units[1.2345], 1.2345)
        self.assertAlmostEqual(units['kJ/mol'], 1/2625.499)
        self.assertAlmostEqual(units['kj/mol'], 1/2625.499)
        self.assertAlmostEqual(units['KJ/MOL'], 1/2625.499)
        self.assertAlmostEqual(units[1/2625.499], 1/2625.499)
        with self.assertRaises(KeyError):
            units['badunit']


if __name__ == '__main__':
    unittest.main()

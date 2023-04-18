import unittest

import pymolpro
from pymolpro import database
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
        db = database.load('sample')
        self.assertEqual(len(db.molecules), 4)
        self.assertEqual(len(db.reactions), 2)
        raw_path = os.path.realpath(os.path.join(__file__, '..', '..', 'pymolpro', 'share', 'database', 'sample'))
        db_raw_path = database.load(raw_path)
        self.assertEqual(db, db_raw_path)
        raw_path = os.path.realpath(os.path.join(__file__, '..', '..', 'pymolpro', 'share', 'database', 'sample.json'))
        db_raw_path = database.load(raw_path)
        self.assertEqual(db, db_raw_path)
        db_string = database.load(' '.join(open(raw_path, 'r').readlines()))
        self.assertEqual(db, db_string)

    def test_run_database(self):
        if shutil.which('molpro'):
            db = database.load('sample')
            results = database.run(db, method='df-lmp2', basis='aug-cc-pVDZ', func="energy")
            results2 = database.run(db, method='df-lmp2', basis='aug-cc-pVDZ', func="energy")
            # del results['projects']
            # del results2['projects']
            self.assertEqual(results.molecule_energies, results2.molecule_energies)
            self.assertEqual(results, results2)
            results2.molecule_energies['HF'] = 12345.0
            self.assertNotEqual(results.molecule_energies, results2.molecule_energies)
            self.assertNotEqual(results, results2)
            database.run(db, method='df-lmp2', basis='aug-cc-pVDZ', func="energy", clean=True)
            # print(results)

    def test_preamble(self):
        if shutil.which('molpro'):
            db = Database()
            db.preamble = 'angstrom'
            db.add_molecule('H2', 'H;H,H,0.7'
                            # , preamble='angstrom'
                            )
            # print(db)
            db2 = database.load(db.dump())
            # print(db2)
            results = database.run(db2)
            # print(results)
            self.assertAlmostEqual(results.molecule_energies['H2'], -1.13207566548333)

    def test_compare_database_runs(self):
        if shutil.which('molpro'):
            db = database.load('sample')
            results = [
                database.run(db, method='df-lmp2', basis='aug-cc-pVDZ', func="energy"),
                database.run(db, method='df-lmp2', basis='aug-cc-pVTZ', func="energy"),
                database.run(db, method='df-lmp2', basis='aug-cc-pVQZ', func="energy"),
            ]
            outputs = database.analyse(results, results[-1])
            self.assertEqual(len(outputs['reaction statistics'].index), 4)
            self.assertEqual(len(outputs['molecule statistics'].index), 4)
            self.assertEqual(len(outputs['reaction energies'].index), len(db))
            for title, output in outputs.items():
                self.assertEqual(len(output.columns), 3)
                # print(title)
                # print(output)
            db.reaction_energies = results[-1].reaction_energies
            db.molecule_energies = {}
            self.assertEqual(len(database.analyse(results, db)), 4)

    def test_fail(self):
        if shutil.which('molpro'):
            db = database.load('sample')

            result = database.run(db, method='hf', basis='minao')
            self.assertFalse(result.failed)
            shutil.rmtree(result.project_directory)

            result = database.run(db, method='hf', preamble='memory,1,k', basis='minao')
            self.assertTrue(result.failed)
            self.assertEqual(len(result.failed), len(db.molecules))
            self.assertEqual(result.failed['HF'].status, 'failed')
            shutil.rmtree(result.project_directory)

            result = database.run(db, method='ccsd', preamble='charge=-10', basis='cc-pvdz')
            self.assertTrue(result.failed)
            self.assertEqual(len(result.failed), 1)
            self.assertEqual(result.failed['HF'].status, 'failed')
            shutil.rmtree(result.project_directory)

            result = database.run(db, method='bad-method', basis='minao')
            self.assertEqual(len(result.failed), len(db.molecules))
            self.assertEqual(result.failed['HF'].status, 'failed')
            shutil.rmtree(result.project_directory)

    def test_subset(self):
        db = database.load('sample')
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
            result = database.run(subset, clean=True)
            self.assertEqual(len(result.molecule_energies), 2)
            self.assertEqual(len(result.reaction_energies), 1)

    def test_extrapolate(self):
        if shutil.which('molpro'):
            db = database.load('sample').subset('non-covalent')
            hfresults = []
            results = []
            for basis in ['cc-pvdz', 'cc-pvtz']:
                hfresults.append(database.run(db, 'hf', basis))
                results.append(database.run(db, 'mp2', basis))
            results += database.basis_extrapolate(results, hfresults, [2, 3])
            results += database.basis_extrapolate(results, hfresults)
            # print(database.analyse(results))
            self.assertEqual(len(results), 4)

            hfresults = []
            results = []
            for basis in ['cc-pvdz', 'cc-pvtz', 'cc-pvqz']:
                hfresults.append(database.run(db, 'hf', basis))
                results.append(database.run(db, 'mp2', basis))
            results += database.basis_extrapolate(results, hfresults)
            # print(database.analyse(results))
            self.assertEqual(len(results), 5)

    def test_units(self):
        from pymolpro.database import units
        self.assertEqual(units[None], 1.0)
        self.assertEqual(units[1.2345], 1.2345)
        self.assertAlmostEqual(units['kJ/mol'], 1 / 2625.499)
        self.assertAlmostEqual(units['kj/mol'], 1 / 2625.499)
        self.assertAlmostEqual(units['KJ/MOL'], 1 / 2625.499)
        self.assertAlmostEqual(units[1 / 2625.499], 1 / 2625.499)
        with self.assertRaises(KeyError):
            units['badunit']
        # print(units)

    def test_add(self):
        db = database.load('sample')
        db += db
        # print(db)
        self.assertEqual(db, db)
        db1 = database.load('sample').subset(['H2+F2'])
        db2 = database.load('sample').subset(['HFHF'])
        db12 = db1.copy()
        db12 += db2
        # print(db12)
        self.assertEqual(len(db12), 2)
        db12 = db1 + db2
        # print(db12)
        self.assertEqual(len(db12), 2)

    def test_library(self):
        from pymolpro.database import library
        self.assertEqual(len(library('GMTKN55')), 55)
        self.assertEqual(len(library('^.MTKN55')), 55)
        self.assertEqual(len(library('GMTKN55$')), 0)
        self.assertEqual(len(library('sample')), 1)
        self.assertGreater(len(library()), 55)
        self.assertEqual(len(library('^((?!GMTKN55).)*$')) + len(library('GMTKN55')), len(library()))

    def test_no_reactions(self):
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
        db.molecule_energies['H2'] = 1.0
        db.molecule_energies['F2'] = 2.0
        db.molecule_energies['HF'] = 3.0
        db.molecule_energies['HFHF'] = 4.0
        db2 = db.copy()
        for molecule in db.molecules:
            db2.molecule_energies[molecule] = db.molecule_energies[molecule] + 0.01
        # print(db)
        # print(db2)
        statistics_ = pymolpro.database.analyse(db2, db)['molecule statistics']
        # print(statistics_[0])
        self.assertAlmostEqual(statistics_[0]['stdev'],0.0)
        self.assertAlmostEqual(statistics_[0]['maxabs'],1e-2)
        statistics_ = pymolpro.database.analyse(db2, db)['molecule energies']
        # print(statistics_)
        self.assertAlmostEqual(statistics_[0]['H2'],1.01)


if __name__ == '__main__':
    unittest.main()

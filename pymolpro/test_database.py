import unittest

import pymolpro
import pandas
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

    def test_bad_library(self):
        import pytest
        with pytest.raises(Exception) as e_info:
            database.load('junk')
        self.assertTrue(isinstance(e_info.value, ValueError))
        self.assertEqual(str(e_info.value),
                         'Cannot resolve "junk" as a library key, library file name, or library-defining json string')

    def test_run_database(self):
        if pymolpro.Project('test').local_molpro_root:
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

    def test_run_opt(self):
        if pymolpro.Project('test').local_molpro_root:
            db = database.run(database.load('sample'), method='hf', basis='minao')
            db_opt = database.run(db, method='hf', basis='minao', func='opt')
            db_opt_test = database.run(db_opt, method='hf', basis='minao', func='energy')
            # print('H2 db', db.project_directory, db.molecules['H2']['geometry'])
            # print('H2 db_opt', db_opt.project_directory, db_opt.molecules['H2']['geometry'])
            # print('H2 db_opt_test', db_opt_test.project_directory, db_opt_test.molecules['H2']['geometry'])
            self.assertNotEqual(db.molecules['H2']['geometry'], db_opt_test.molecules['H2']['geometry'])
            self.assertEqual(db_opt.molecules['H2']['geometry'], db_opt_test.molecules['H2']['geometry'])

    def test_preamble(self):
        if pymolpro.Project('test').local_molpro_root:
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

    def test_orbitals(self):
        if pymolpro.Project('test').local_molpro_root:
            db = Database()
            db.preamble = 'angstrom'
            db.add_molecule('H2', 'H;H,H,0.7'
                            # , preamble='angstrom'
                            )
            results = database.run(db, postamble='put,xml')
            p_ = results.projects['H2']
            m = p_.xpath('//molecule')[-1]
            print(m)
            orbitalSets = p_.xpath('orbitals', m)
            print(orbitalSets)
            print(p_.xml)
            print(p_.orbitals())


    def test_compare_database_runs(self):
        if pymolpro.Project('test').local_molpro_root:
            db = database.load('sample')
            results = [
                database.run(db, method='df-lmp2', basis='aug-cc-pVDZ', func="energy"),
                database.run(db, method='df-lmp2', basis='aug-cc-pVTZ', func="energy"),
                database.run(db, method='df-lmp2', basis='aug-cc-pVQZ', func="energy"),
            ]
            outputs = database.analyse(results, results[-1])
            # print(outputs)
            self.assertEqual(len(outputs['reaction statistics'].index), 5)
            self.assertEqual(len(outputs['molecule statistics'].index), 5)
            self.assertEqual(len(outputs['reaction energies'].index), len(db))
            rmsd = outputs['molecule statistics'].iloc[:, 0]['RMSD']
            msd = outputs['molecule statistics'].iloc[:, 0]['MSD']
            stdevd = outputs['molecule statistics'].iloc[:, 0]['STDEVD']
            self.assertAlmostEqual((len(db.molecules) - 1) * pow(stdevd, 2),
                                   len(db.molecules) * (pow(rmsd, 2) - pow(msd, 2)))
            for title, output in outputs.items():
                if type(output) == pandas.DataFrame: self.assertEqual(len(output.columns), 3)
                # print(title)
                # print(output)
            db.reaction_energies = results[-1].reaction_energies
            db.molecule_energies = {}
            self.assertEqual(len(database.analyse(results, db)), 7)

    def test_fail(self):
        if pymolpro.Project('test').local_molpro_root:
            db = database.load('sample')

            result = database.run(db, method='hf', basis='minao')
            self.assertFalse(result.failed)
            shutil.rmtree(result.project_directory)

            result = database.run(db, method='hf', preamble='memory,1,m', basis='minao')
            self.assertTrue(result.failed)
            self.assertEqual(len(result.failed), len(db.molecules))
            self.assertEqual(result.failed['HF'].status, 'failed')
            shutil.rmtree(result.project_directory)

            try:
                result = database.run(db, method='ccsd', preamble='charge=-10', basis='cc-pvdz')
            except:
                pass
            self.assertTrue(result.failed)
            self.assertEqual(len(result.failed), 1)
            self.assertEqual(result.failed['HF'].status, 'failed')
            shutil.rmtree(result.project_directory)

            try:
                result = database.run(db, method='bad-method', basis='minao')
            except:
                pass
            self.assertNotEqual(len(result.failed), 0)
            # self.assertEqual(result.failed['HF'].status, 'failed')
            try:
                shutil.rmtree(result.project_directory)
            except:
                pass

    def test_subset(self):
        db = database.load('sample')
        self.assertEqual(len(db.subset('non-covalent')), 1)
        self.assertEqual(len(db.subset('non-covalent').molecules), 2)
        # print(db.subset('non-covalent'))

        subset = db.subset()
        self.assertEqual(subset, db)

        subset = db.subset(max_atoms=2)
        # print(db)
        # print(subset)
        self.assertEqual(len(subset), 1)
        # print("description:", db.description)
        # print("description:", subset.description)

        subset = db.subset(['HFHF'])
        # print(subset)
        self.assertEqual(len(subset), 1)
        self.assertEqual(len(subset.molecules), 2)

        db.add_subset('vdW', 'HFHF')
        self.assertEqual(len(db.subset('vdW')), 1)
        self.assertEqual(len(db.subset('vdW').molecules), 2)
        # print(db.subset('vdW'))

        if pymolpro.Project('test').local_molpro_root:
            result = database.run(subset, clean=True)
            self.assertEqual(len(result.molecule_energies), 2)
            self.assertEqual(len(result.reaction_energies), 1)

        subset = pymolpro.database.load('GMTKN55_BH76').subset(open_shell=False)
        self.assertEqual(len(subset), 20)

        # print(pymolpro.database.elements(list(subset.molecules.values())[0]['geometry']))
        # print(pymolpro.database.electrons(list(subset.molecules.values())[0]['geometry']))
        self.assertEqual(len(pymolpro.database.elements(list(subset.molecules.values())[0])),5)
        self.assertEqual(pymolpro.database.electrons(list(subset.molecules.values())[0]['geometry']),18)
        self.assertEqual(pymolpro.database.electrons(list(subset.molecules.values())[0]),18)

        subset = pymolpro.database.load('GMTKN55_BH76').subset(max_electrons=10)
        # print(subset)
        # for k,m in subset.molecules.items():
        #     print(k,pymolpro.database.electrons(m))
        self.assertEqual(len(subset), 4)

    def test_extrapolate(self):
        if pymolpro.Project('test').local_molpro_root:
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

    def test_cache(self):
        if pymolpro.Project('test').local_molpro_root:
            db=pymolpro.database.Database()
            db.add_molecule('H2','H 0 0 .35\nH 0 0 -.35')
            newdb=pymolpro.database.run(db, check=False)
            first_filename = newdb.projects['H2'].filename()
            newdb=pymolpro.database.run(db, check=False)
            second_filename = newdb.projects['H2'].filename()
            shutil.rmtree(newdb.project_directory)
            self.assertEqual(first_filename,second_filename) # because the second run should not have been done because input identical

    def test_vector_energy(self):
        if pymolpro.Project('test').local_molpro_root:
            db=pymolpro.database.Database()
            db.add_molecule('CH2-a','C;H,C,2;H,C,2,H,110')
            db.add_molecule('CH2-b','C;H,C,2;H,C,2,H,120')
            db.add_reaction('bend',{'CH2-a':-1,'CH2-b':1})
            pymolpro.database.run(db, method='ccsd(t)-F12',basis='cc-pVDZ') # to cover the case that the Molpro ENERGY variable is a vector

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
        self.assertAlmostEqual((len(db.molecules) - 1) * pow(statistics_[0]['STDEVD'], 2),
                               len(db.molecules) * (pow(statistics_[0]['RMSD'], 2) - pow(statistics_[0]['MSD'], 2)))
        self.assertAlmostEqual(statistics_[0]['STDEVD'], 0.0)
        self.assertAlmostEqual(statistics_[0]['MAXD'], 1e-2)
        energies_ = pymolpro.database.analyse(db2, db)['molecule energies']
        # print(energies_)
        self.assertAlmostEqual(energies_[0]['H2'], 1.01)


if __name__ == '__main__':
    unittest.main()

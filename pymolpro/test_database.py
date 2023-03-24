import unittest

import pymolpro
from pymolpro import Database
import os


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
            'H2': 'H 0 0 0\nH 0 0 .7',
            'F2': 'F 0 0 0\nF 0 0 1',
            'HF': 'H 0 0 0\nF 0 0 1',
        })
        self.assertEqual(len(db.molecules), 3)  # add assertion here
        self.assertEqual(len(db.reactions), 0)  # add assertion here

    def test_json(self):
        db = Database(
            molecules={
                'H2': 'H 0 0 0\nH 0 0 .7',
                'F2': 'F 0 0 0\nF 0 0 1',
                'HF': 'H 0 0 0\nF 0 0 1',
            },
        )
        db.add_reaction('H2+F2',
                        {'HF': 2,
                         'H2': -1,
                         'F2': -1,
                         }, description='Hydrogenation of difluorine')
        db.dump("sample.json")
        db2 = Database()
        db2.load("sample.json")
        os.remove("sample.json")
        self.assertEqual(db.molecules, db2.molecules)
        self.assertEqual(db.reactions, db2.reactions)

    def test_library(self):
        db = pymolpro.library_database('sample')
        self.assertEqual(len(db.molecules), 3)
        self.assertEqual(len(db.reactions), 1)


if __name__ == '__main__':
    unittest.main()

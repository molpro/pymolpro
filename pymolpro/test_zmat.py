import unittest
import pymolpro


class TestZmat(unittest.TestCase):
    def test_zmat_string(self):
        z = pymolpro.xyz_to_zmat('2\ntest\nF 0.0 0.0 0.0\nH 0.0 0.0 1.0')
        assert z == 'F1\nH2 F1 1.0\n'

    def test_zmat_file(self):
        import tempfile
        with tempfile.NamedTemporaryFile(suffix='.xyz') as f:
            f.write('2\ntest\nF 0.0 0.0 0.0\nH 0.0 0.0 1.0'.encode('utf-8'))
            f.flush()
            z = pymolpro.xyz_to_zmat(f.name)
            assert z == 'F1\nH2 F1 1.0\n'

    def test_zmat_1atom(self):
        z = pymolpro.xyz_to_zmat('1\ntest with 1 atom\nO 0.000000 0.000000  0.000000\n')
        assert z == 'O1\n'

    def test_zmat_2atoms(self):
        z = pymolpro.xyz_to_zmat('2\ntest with 2 atoms\nO 0.000000 0.000000  0.000000\nH 0.758602 0.000000  0.504284\n')
        assert z == 'O1\nH2 O1 0.910922\n'

    def test_zmat_3atoms(self):
        z = pymolpro.xyz_to_zmat(
            '3\ntest with 3 atoms\nO 0.000000 0.000000  0.000000\nH 0.758602 0.000000  0.504284\nH 0.260455 0.000000 -0.872893\n')
        assert z == 'O1\nH2 O1 0.910922\nH3 O1 0.910922 H2 107.000024\n'

    def test_zmat_4atoms(self):
        z = pymolpro.xyz_to_zmat(
            '4\ntest with 4 atoms\nO 0.000000 0.000000  0.000000\nH 0.758602 0.000000  0.504284\nH 0.260455 0.000000 -0.872893\nO 3.000000 0.500000  0.000000\n')
        assert z == 'O1\nH2 O1 0.910922\nH3 O1 0.910922 H2 107.000024\nO4 H2 2.351206 O1 132.466298 H3 -16.755013\n'

    def test_zmat_5atoms(self):
        z = pymolpro.xyz_to_zmat(
            '5\ntest with 5 atoms\nO 0.000000 0.000000  0.000000\nH 0.758602 0.000000  0.504284\nH 0.260455 0.000000 -0.872893\nO 3.000000 0.500000  0.000000\nH 3.758602 0.500000  0.504284\n')
        assert z == 'O1\nH2 O1 0.910922\nH3 O1 0.910922 H2 107.000024\nO4 H2 2.351206 O1 132.466298 H3 -16.755013\nH5 O4 0.910922 H2 132.466298 O1 180.000000\n'

    def test_zmat_6atoms(self):
        z = pymolpro.xyz_to_zmat(
            '6\ntest with 6 atoms\nO 0.000000 0.000000  0.000000\nH 0.758602 0.000000  0.504284\nH 0.260455 0.000000 -0.872893\nO 3.000000 0.500000  0.000000\nH 3.758602 0.500000  0.504284\nH 3.260455 0.500000 -0.872893\n')
        assert z == 'O1\nH2 O1 0.910922\nH3 O1 0.910922 H2 107.000024\nO4 H2 2.351206 O1 132.466298 H3 -16.755013\nH5 O4 0.910922 H2 132.466298 O1 180.000000\nH6 O4 0.910922 H5 107.000024 H2 163.244987\n'

    def test_zmat_glycine(self):
        z = pymolpro.xyz_to_zmat(
            '10\ntest with glycine\nH  0.000000000000 -0.106159785893 -2.351004690069\nO  0.000000000000 -0.678075052479 -1.591046451359\nC  -0.000000000000 0.082367016245 -0.500909038651\nO  -0.000000000000 1.267791067549 -0.543752808818\nC  -0.000000000000 -0.744221268281 0.765136671851\nH  -0.870415825720 -1.402324054216 0.725864010781\nH  0.870415825720 -1.402324054216 0.725864010781\nN  0.000000000000 0.013539178369 1.987713087169\nH  0.799232626780 0.624397043016 2.007552635709\nH  -0.799232626780 0.624397043016 2.007552635709\n')
        assert z == 'C5\nN8 C5 1.438365\nC3 C5 1.511992 N8 115.069000\nH6 C5 1.091909 N8 110.370557 C3 121.752502\nH7 C5 1.091909 N8 110.370557 C3 -121.752502\nH9 N8 1.006138 C5 109.670446 C3 57.520627\nH10 N8 1.006138 C5 109.670446 C3 -57.520627\nO2 C3 1.329162 C5 111.961600 N8 -180.000000\nO4 C3 1.186198 C5 125.210000 N8 -0.000000\nH1 O2 0.951117 C3 108.138000 C5 -180.000000\n'

    def test_xyz_to_zmat_from_chemcoord(self):
        z = str(pymolpro.geometry.convert_xyz_to_chemcoordzmat(
            '6\ntest with 6 atoms\nO 0.000000 0.000000  0.000000\nH 0.758602 0.000000  0.504284\nH 0.260455 0.000000 -0.872893\nO 3.000000 0.500000  0.000000\nH 3.758602 0.500000  0.504284\nH 3.260455 0.500000 -0.872893\n'))
        self.assertMultiLineEqual(z,
                                  '  atom       b      bond    a       angle    d    dihedral\n1    O  origin  0.000000  e_z    0.000000  e_x    0.000000\n2    H       1  0.910922  e_z   56.385853  e_x   -0.000000\n3    H       1  0.910922    2  107.000024  e_x   -0.000000\n4    O       2  2.351206    1  132.466298    3  -16.755013\n5    H       4  0.910922    2  132.466298    1  180.000000\n6    H       4  0.910922    5  107.000024    2  163.244987')

    def test_from_chemcoordzmat_to_molprozmat(self):
        z = pymolpro.geometry.convert_chemcoordzmat_to_molprozmat(
            '  atom       b      bond    a       angle    d    dihedral\n1    O  origin  0.000000  e_z    0.000000  e_x    0.000000\n2    H       1  0.910922  e_z   56.385853  e_x   -0.000000\n3    H       1  0.910922    2  107.000024  e_x   -0.000000\n4    O       2  2.351206    1  132.466298    3  -16.755013\n5    H       4  0.910922    2  132.466298    1  180.000000\n6    H       4  0.910922    5  107.000024    2  163.244987')
        assert z == 'O1\nH2 O1 0.910922\nH3 O1 0.910922 H2 107.000024\nO4 H2 2.351206 O1 132.466298 H3 -16.755013\nH5 O4 0.910922 H2 132.466298 O1 180.000000\nH6 O4 0.910922 H5 107.000024 H2 163.244987\n'

    def test_from_chemcoordzmat_with_dummies_to_molprozmat(self):
        z = pymolpro.geometry.convert_chemcoordzmat_to_molprozmat(
            '  atom       b      bond    a       angle    d    dihedral\n1    O  origin  3.825100  e_z   97.575672  e_x -172.422536\n2    H       1  0.910922  e_z   13.716615  e_x -146.597056\n3    H       1  0.910922    2  107.000024  e_x  -12.723783\n4    O       2  2.351206    1  180.000000    3  -16.755013\n7    X  origin  2.604064  e_z  107.831880  e_x -175.100590\n5    H       4  0.910922    2  132.466298    7 -180.000000\n8    X  origin  2.604064  e_z  107.831880  e_x -175.100590\n6    H       4  0.910922    2  118.561123    8   18.293238\n')
        self.assertIsNone(z)


if __name__ == '__main__':
    unittest.main()

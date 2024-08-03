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


if __name__ == '__main__':
    unittest.main()

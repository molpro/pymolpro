import string
import re
import os

def xyz_to_zmat(xyz:str, algorithm='chemcoord'):
    """
    Converts from xyz to z-matrix coordinates.

    :param xyz: xyz-file
    :param algorithm: algorithm to choose, default is the chemcoord algorithm
    :return: The Z matrix as a string, or None in case of unknown algorithm,
    :rtype: str
    """

    if (algorithm == 'chemcoord'):
        return xyz_via_chemcoord_to_molprozmat(xyz)
    return None

def xyz_via_chemcoord_to_molprozmat(xyz):
    """
    Converts from xyz to z-matrix coordinates, using chemcoord-algorithm

    :param xyz: xyz-file
    :return: molprozmat , no check for errors
    """
    chemcoordzmat = convert_xyz_to_chemcoordzmat(xyz)
    # print ('z-matrix from chemcoord:')
    # print(chemcoordzmat)
    molprozmat = convert_chemcoordzmat_to_molprozmat(chemcoordzmat)
    # print('\nz-matrix for Molpro:')
    # print(molprozmat)
    return molprozmat

def convert_xyz_to_chemcoordzmat(xyzinput):
    """
    Converts from xyz to z-matrix coordinates in chemcoord-format

    :param xyz: xyz-file
    :return: chemcoordzmat , no check for errors
    """
    import chemcoord
    if os.access(xyzinput, os.R_OK):
        xyzfile = chemcoord.Cartesian.read_xyz(xyzinput, start_index = 1)
    else:
        import tempfile
        with tempfile.NamedTemporaryFile(suffix='.xyz') as f:
            f.write(xyzinput.encode('utf-8'))
            f.flush()
            xyzfile = chemcoord.Cartesian.read_xyz(f.name, start_index = 1)
    chemcoordzmat = xyzfile.get_zmat()
    return chemcoordzmat

def convert_chemcoordzmat_to_molprozmat(chemcoordzmat):
    """
    Converts z-matrix from chemcoord format to molpro format
    :param chemcoordzmat: z-matrix in chemcoord format
    :return: zmatrix , no check for errors
    assumption is that axes e_z and e_x
    are only used for the first two atoms (this may not be the case if
    the ininital z-matrix form chemcoord is further transformed,
    and dummies appear in the chemcoord z-matrix)
    """
    f = str(chemcoordzmat)
    n = 0
    for x in f.splitlines():
        n = n +1
        if (n == 1):
            continue
        atom1 = x.split()[1] + x.split()[0]
        if (n == 2):
            atomdict = {x.split()[0] : x.split()[1]}
            zmatrix = atomdict.get(x.split()[0]) + x.split()[0]+'\n'
        atomdict[x.split()[0]] = x.split()[1]
        if (n > 2):
            atom2 = atomdict.get(x.split()[2]) + (x.split()[2])
        if (n == 3):
            line2 = atom1 + " " + atom2 + " " + x.split()[3]
            zmatrix = zmatrix + line2+'\n'
        if (n > 3):
            atom3 = str(atomdict.get(x.split()[4])) + x.split()[4]
        if (n == 4):
            line3 = atom1 + " " + atom2 + " " + x.split()[3] + " " + atom3 + " " + x.split()[5]
            zmatrix = zmatrix + line3+'\n'
        atom4 = str(atomdict.get(x.split()[6])) + x.split()[6]
        if (n > 4):
            line4 = atom1 + " " + atom2 + " " + x.split()[3] + " " + atom3 + " " + x.split()[5] + " " + atom4 + " " + x.split()[7]
            zmatrix = zmatrix + line4+'\n'
    return zmatrix

if (__name__ == "__main__"):
    """
    when called from the command line
    :param xyz: xyz-file
    :param algorithm: algorithm to choose, Default is the chemcoord-algorithm
    :example: python pymolpro/geometry.py ../../tests-fuer-xyz-to-zmat-conversion/glycine.xyz
    """
    import sys
    import os
    if len(sys.argv) >1:
        if (os.path.isfile(sys.argv[1])):
            xyz_to_zmat(sys.argv[1], algorithm='chemcoord')
    else:
        pass

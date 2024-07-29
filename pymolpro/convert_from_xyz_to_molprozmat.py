import string
import re
import chemcoord

def convert_xyz_to_molprozmat(xyzinput):
    chemcoordzmat = convert_xyz_to_chemcoordzmat(xyzinput)
    print ('z-matrix from chemcoord:')
    print(chemcoordzmat)
    molprozmat = convert_chemcoordzmat_to_molprozmat(chemcoordzmat)
    print('\nz-matrix for Molpro:')
    print(molprozmat)
    return molprozmat

def convert_xyz_to_chemcoordzmat(xyzinput):
    xyzfile = chemcoord.Cartesian.read_xyz(xyzinput, start_index = 1)
    chemcoordzmat = xyzfile.get_zmat()
    return chemcoordzmat

def convert_chemcoordzmat_to_molprozmat(chemcoordzmat):
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


#!/usr/bin/env python
import math
import re
from lxml import etree
import numpy as np

namespaces = {
    'molpro-output': 'http://www.molpro.net/schema/molpro-output',
    'xsd': 'http://www.w3.org/1999/XMLSchema',
    'cml': 'http://www.xml-cml.org/schema',
    'stm': 'http://www.xml-cml.org/schema',
    'xhtml': 'http://www.w3.org/1999/xhtml',
    'xlink': 'http://www.w3.org/1999/xlink'}


def evaluateBasis(molecule, points):
    """
    Evaluate the orbital basis set on a grid

    :param molecule: lxml.etree.Element holding the molecule
    :param points: numpy :,3 array
    :return: numpy :,: array, first index basis function, second index brid point
    """
    cartesianAngularQuantumNumbers = np.array(
        [[0, 1, 0, 0, 2, 0, 0, 1, 1, 0, 3, 0, 0, 2, 2, 1, 0, 1, 0, 1, 4, 0, 0, 3, 3, 1, 0, 1, 0, 2, 2, 0, 2, 1, 1, 5, 4,
          4, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1,
          1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 1, 0, 0, 2, 0, 1, 0, 1, 0, 3, 0, 1, 0, 2, 2, 0, 1, 1, 0, 4, 0, 1, 0, 3, 3, 0, 1, 2, 0, 2, 1, 2, 1, 0, 1,
          0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0, 0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 5, 4, 3,
          2, 1, 0, 6, 5, 4, 3, 2, 1, 0],
         [0, 0, 0, 1, 0, 0, 2, 0, 1, 1, 0, 0, 3, 0, 1, 0, 1, 2, 2, 1, 0, 0, 4, 0, 1, 0, 1, 3, 3, 0, 2, 2, 1, 1, 2, 0, 0,
          1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2,
          3, 4, 5, 0, 1, 2, 3, 4, 5, 6]])
    normfac = np.array(
        [1, 1, 1, 1, 3, 3, 3, 1, 1, 1, 15, 15, 15, 3, 3, 3, 3, 3, 3, 1, 105, 105, 105, 15, 15, 15, 15, 15, 15, 9, 9, 9,
         1, 1, 1, 945, 105, 105, 45, 15, 45, 45, 9, 9, 45, 105, 15, 9, 15, 105, 945, 105, 45, 45, 105, 945])
    sphtran = []
    sphtran.append(np.array([[1]]))
    sphtran.append(np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
    sphtran.append(np.array(
        [[-0.5, -0.5, 1, 0, 0, 0], [1, -1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0], [0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1]]))

    for answer in molecule.xpath(
            'molpro-output:variables/molpro-output:variable[@name="_ANGSTROM"]/molpro-output:value',
            namespaces=namespaces):
        Angstrom = np.float64(answer.text)
    atoms = molecule.xpath('cml:molecule/cml:atomArray/cml:atom', namespaces=namespaces)

    basisSets = molecule.xpath('molpro-output:basisSet[@id="ORBITAL"]', namespaces=namespaces)
    if len(basisSets) != 1:
        raise Exception('something is wrong: there should be just one orbital basisSet')
    basisSet = basisSets[0]

    # print('Basis length =',basisSet.get('length'),' angular =',basisSet.get('angular'),' type =',basisSet.get('type'),' groups =',basisSet.get('groups'))
    basisGroups = basisSet.xpath('molpro-output:basisGroup', namespaces=namespaces)
    #    print
    # print(str(len(basisGroups))+' basisGroups:')
    # for basisGroup in basisGroups:
    #    print(basisGroup.get('id')+' '+basisGroup.get('minL')+' '+basisGroup.get('maxL')+' '+basisGroup.get('primitives')+' '+basisGroup.get('angular')+' '+basisGroup.get('contractions'))
    # print

    # now walk through basis set and evaluate it at some points
    basisAtPoints = np.empty((0, len(points)), dtype=np.float64)
    nbasis = 0
    iatom = 0
    for atom in atoms:
        nuclearCoordinates = [np.float64(atom.get('x3')) * Angstrom, np.float64(atom.get('y3')) * Angstrom,
                              np.float64(atom.get('z3')) * Angstrom]
        query = 'molpro-output:association/molpro-output:atoms[@xlink:href[contains(.,"@id=\'' + atom.get(
            'id') + '\'")]]/..'
        basisGroupAssociation = basisSet.xpath(query, namespaces=namespaces)
        if len(basisGroupAssociation) != 1:
            raise Exception('something is wrong: there should be a unique association of an atom with a basis set')
        bases = basisGroupAssociation[0].xpath('molpro-output:bases', namespaces=namespaces)
        if len(bases) != 1:
            raise Exception('something is wrong: there should be a bases node in association')
        basesString = bases[0].get('{' + namespaces['xlink'] + '}href')
        basesString = basesString[basesString.find('basisGroup['):]
        basesString = basesString[basesString.find('[') + 1:].lstrip()
        basesString = basesString[:basesString.find(']')].rstrip()
        basesString = basesString.replace('or', '').replace('\n', '').replace("'", '')
        list = basesString.split("@id=");
        for item in list:
            item = item.lstrip().rstrip()
            if item.isalnum():
                basisGroup = basisSet.xpath('molpro-output:basisGroup[@id="' + item + '"]', namespaces=namespaces)[0]
                lquant = int(basisGroup.get('minL'))
                if lquant > 5:
                    raise Exception("Sorry, I was too lazy to write this for i basis functions and higher")
                if lquant != int(basisGroup.get('maxL')):
                    raise Exception("This program cannot handle multiple-angular momentum sets")
                # print(basisGroup.get('id')+' '+basisGroup.get('minL')+' '+basisGroup.get('maxL')+' '+basisGroup.get('primitives')+' '+basisGroup.get('angular')+' '+basisGroup.get('contractions'))
                alpha = np.float64(re.sub(' +', ' ',
                                          basisGroup.xpath('molpro-output:basisExponents', namespaces=namespaces)[
                                              0].text.replace('\n', '').lstrip().rstrip()).split(" "))
                # first evaluate all the primitives for this shell on the grid
                ncompc = int(((lquant + 2) * (lquant + 1)) / 2)  # number of cartesian components
                primitivesc = np.empty((ncompc, len(alpha), len(points)))
                lqbase = (lquant * (lquant + 1) * (lquant + 2)) // 6
                for ip in range(len(points)):
                    xyz = np.subtract(points[ip], nuclearCoordinates)
                    r2 = np.dot(xyz, xyz)
                    for ia in range(len(alpha)):
                        alph = alpha[ia]
                        norm = math.sqrt(math.sqrt(2 * alph / math.pi) ** 3 * (4 * alph) ** lquant)
                        value = norm * math.exp(-alph * r2)
                        for icomp in range(ncompc):
                            k = cartesianAngularQuantumNumbers[0, icomp + lqbase]
                            l = cartesianAngularQuantumNumbers[1, icomp + lqbase]
                            m = cartesianAngularQuantumNumbers[2, icomp + lqbase]
                            primitivesc[icomp, ia, ip] = value * xyz[0] ** k * xyz[1] ** l * xyz[2] ** m / math.sqrt(
                                normfac[icomp + lqbase])

                # transformation to spherical harmonics
                if molecule.xpath('molpro-output:orbitals', namespaces=namespaces)[0].get(
                        'angular') == 'spherical':  # Molpro 2012.1 does not produced this, but cartesian only. The spherical code here is not finished, but not presently needed
                    ncomp = 2 * lquant + 1
                else:
                    ncomp = ncompc
                if ncomp < ncompc and lquant >= len(sphtran):
                    raise Exception("Spherical functions not yet coded")
                elif ncomp < ncompc:
                    primitives = np.empty((ncomp, len(alpha), len(points)))
                    # print("primitivesc",primitivesc)
                    # print("sphtran[lquant]",sphtran[lquant])
                    for ip in range(len(points)):
                        # print sphtran[lquant]
                        # print primitivesc[:,:,ip]
                        # print sphtran[lquant].shape
                        # print primitivesc[:,:,ip].shape
                        primitives[:, :, ip] = np.dot(sphtran[lquant], primitivesc[:, :, ip])
                    # print "primitives",primitives
                else:
                    primitives = primitivesc

                # next loop over contractions
                for basisContraction in basisGroup.xpath('molpro-output:basisContraction', namespaces=namespaces):
                    cc = np.float64(
                        re.sub(' +', ' ', basisContraction.text.replace('\n', '').lstrip().rstrip()).split(" "))
                    # loop over angular components in the contraction
                    basisAtPoints = np.append(basisAtPoints, np.empty((ncomp, len(points)), dtype=np.float64), axis=0)
                    for m in range(ncomp):
                        basisAtPoints[nbasis][:] = np.dot(cc, primitives[m, :, :])
                        nbasis += 1  # completed evaluating this basis function
        iatom += 1
    # print outer(basisAtPoints[:,1],basisAtPoints[:,1])
    return basisAtPoints  # end basis evaluation


def evaluateOrbitals(molecule, points, minocc=1.0, ID=None, values=False):
    """
    Evaluate the molecular orbitals on a grid

    :param molecule: lxml.etree.Element holding the molecule
    :param points: numpy :,3 array
    :param minocc: Only orbitals with at least this occupation will be returned
    :param ID: Only the orbital whose ID attribute matches this will be selected
    :param values:
    :return: array of dictionaries giving the occupation and values on the grid, or if ID is specified, a single dictionary, or if values==True, a numpy array
    """
    basisAtPoints = evaluateBasis(molecule, points)  # evaluate the basis set at the nuclei
    orbitalSets = molecule.xpath('molpro-output:orbitals', namespaces=namespaces)
    if len(orbitalSets) != 1:
        raise Exception('something is wrong: there should be just one orbital set')
    result = []
    from pymolpro import element_to_dict
    search = 'molpro-output:orbital'
    if ID: search += '[@ID="' + ID + '"]'
    for orbital in orbitalSets[0].xpath(search, namespaces=namespaces):
        occ = np.float64(orbital.get('occupation'))
        if occ >= minocc:
            mos = np.array(re.sub(' +', ' ', orbital.text.lstrip().rstrip().replace('\n', '')).split(" "))
            mopoint = np.zeros(len(points), dtype=np.float64)
            for ic in range(len(basisAtPoints)):
                for ipoint in range(len(points)):
                    mopoint[ipoint] += np.float64(mos[ic]) * basisAtPoints[ic][ipoint]
            orbdict = element_to_dict(orbital)
            orbdict['values'] = mopoint
            orbdict['occ'] = occ
            result.append(orbdict)
    assert not ID or len(result) == 1
    if values:
        return np.array(result[0]['values']) if ID else np.array([res['values'] for res in result])
    else:
        return result[0] if ID else result


def spherical_grid(radial_points, radial_weights, l):
    """
    Take a one-dimensional grid, and use it for radial integration, together with a Lebedev grid of specified order, to generate a 3-dimensional grid (cartesian coordinates, weights) for quadrature.

    :param radial_points:
    :param radial_weights:
    :param l: Largest rank of spherical harmonic for which exact quadrature is wanted.
    :return: tuple of points(numpy array [npt,3]), weights
    """
    import numgrid
    nptang = __lebedev_select(l, True)
    gridang = numgrid.get_angular_grid(nptang)
    # print("gridang:",gridang)
    nptrad = len(radial_points)
    npt3 = nptrad * nptang
    points3d = np.empty([npt3, 3])
    weights3d = np.empty(npt3)
    for j in range(nptrad):
        for k in range(nptang):
            for i in range(3):
                points3d[k + j * nptang, i] = radial_points[j] * gridang[i][k]
            weights3d[k + j * nptang] = 4 * np.pi * pow(radial_points[j], 2) * radial_weights[j] * gridang[3][k]
    # print("in spherical_grid points3d",points3d)
    return [points3d, weights3d]


def __lebedev_select(l, size=False):
    lebedev = {3: 6, 5: 14, 7: 26, 9: 38, 11: 50, 13: 74, 15: 86, 17: 110, 19: 146,
               21: 170, 23: 194, 25: 230, 27: 266, 29: 302, 31: 350, 35: 434, 41: 590, 47: 770,
               53: 974, 59: 1202, 65: 1454, 71: 1730, 77: 2030, 83: 2354, 89: 2702, 95: 3074, 101: 3470,
               107: 3890, 113: 4334, 119: 4802, 125: 5294, 131: 5810}
    for k, v in lebedev.items():
        if l <= k: return v if size else k
    return lebedev[131] if size else 131


def cubical_grid(points1d, weights1d):
    """
    Take a one-dimensional grid, and construct a 3-dimensional outer-product cubical grid

    :param points1d:
    :param weights1d:
    :return: tuple of points(numpy array [npt,3]), weights
    """
    npt1 = len(points1d)
    npt3 = pow(npt1, 3)
    points3d = np.empty([npt3, 3], np.double)
    weights3d = np.empty(npt3, np.double)
    for j in range(npt1):
        for k in range(npt1):
            for l in range(npt1):
                points3d[l + npt1 * (k + npt1 * j), 0] = points1d[j]
                points3d[l + npt1 * (k + npt1 * j), 1] = points1d[k]
                points3d[l + npt1 * (k + npt1 * j), 2] = points1d[l]
                weights3d[l + npt1 * (k + npt1 * j)] = weights1d[j] * weights1d[k] * weights1d[l]
    return [points3d, weights3d]


# example main program
if __name__ == "__main__":

    from sys import argv

    tree = etree.parse(argv[1])
    document = tree.getroot()

    for answer in document.xpath(
            '//molpro-output:variables/molpro-output:variable[@name="_ANGSTROM"]/molpro-output:value',
            namespaces=namespaces):
        Angstrom = np.float64(answer.text)

    molecules = document.xpath('//molpro-output:job/molpro-output:molecule', namespaces=namespaces)
    print(str(len(molecules)) + ' molecules:')
    for molecule in molecules:
        print('Molecule ' + molecule.get('id'))

        metadata = molecule.xpath('stm:metadataList/stm:metadata', namespaces=namespaces)
        print(str(len(metadata)) + ' metadata items:')
        for md in metadata:
            print(md.get('name') + ' = ' + md.get('content'))

        atoms = molecule.xpath('cml:molecule/cml:atomArray/cml:atom', namespaces=namespaces)
        points = np.zeros((len(atoms), 3), dtype=np.float64)
        print(str(len(atoms)) + ' atoms:')
        ids = np.empty(len((atoms)), dtype='a4')
        elementTypes = np.empty((len(atoms)), dtype='a4')
        iatom = 0
        for atom in atoms:
            print(atom.get('id') + ' ' + atom.get('elementType') + ' ' + atom.get('x3') + ' ' + atom.get(
                'y3') + ' ' + atom.get('z3'))
            ids[iatom] = atom.get('id')
            elementTypes[iatom] = atom.get('elementType')
            points[iatom, :] = [np.float64(atom.get('x3')) * Angstrom, np.float64(atom.get('y3')) * Angstrom,
                                np.float64(atom.get('z3')) * Angstrom]
            iatom += 1

        orbitalsAtPoints = evaluateOrbitals(molecule, points)
        results = np.zeros(len(points), dtype=np.float64)
        for orbital in orbitalsAtPoints:
            # print("Orbital",orbital)
            for ipoint in range(len(points)):
                results[ipoint] += orbital['values'][ipoint] ** 2 * orbital['occ']

        print("Calculated densities at the nuclei:")
        print(results)

        for answer in document.xpath('//molpro-output:property[contains(@name,"DELTA")]', namespaces=namespaces):
            name = answer.xpath('@name', namespaces=namespaces)[0]
            value = answer.xpath('@value', namespaces=namespaces)[0]
            print('From Molpro:', name, '=', value)
            closeness = 10000000000
            for ipoint in range(len(points)):
                distance = math.sqrt((results[ipoint] - np.float64(value)) ** 2)
                if distance < closeness:
                    closeness = distance
                    closest = ipoint
            print('... nearest to result ', closest, ' which differs by ', closeness)

        print('End of molecule ' + molecule.get('id'))

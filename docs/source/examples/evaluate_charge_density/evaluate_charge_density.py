from pymolpro import Project
import numpy as np
from pymolpro.grid import evaluateOrbitals

p = Project("evaluate_charge_density")
p.write_input('geometry={f;h,f,1.732};set,sewprop=0;gexpec,delta,f;gexpec,delta,h;rhf;put,xml')
p.run(wait=True)

document = p.xpath('/*')[0]
default_ns_name='molpro-output'
namespaces = {k if k is not None else default_ns_name: v for k, v in document.nsmap.items()}

for answer in document.xpath('//molpro-output:variables/molpro-output:variable[@name="_ANGSTROM"]/molpro-output:value',
                             namespaces=namespaces):
    Angstrom = np.float64(answer.text)

molecules = document.xpath('//molpro-output:job/molpro-output:molecule', namespaces=namespaces)
print(str(len(molecules)) + ' molecules:')
for molecule in molecules:
    print('Molecule ' + molecule.get('id'))

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

    import math
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
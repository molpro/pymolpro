import pymolpro
import os
import re
import lxml
import shutil
import glob

subsets = [
    'benchmarks/N2minimal',
    'Weizmann/W4-11',
    'Karton/CRBH20',
    'Bonn/GMTKN30/BHPERI',
    'Bonn/GMTKN30/DARC',
    'Bonn/GMTKN30/DC9',
    'Bonn/GMTKN30/O3ADD6',
]

directory = 'old_molpro_databases'

MOLPRO_PREFIX = ''
with open(shutil.which('molpro'), 'r') as f:
    while True:
        line = f.readline()
        print(line)
        if re.match('MOLPRO_PREFIX', line):
            MOLPRO_PREFIX = re.sub("'", '', re.sub('MOLPRO_PREFIX=', '', line)).strip()
            break

for subset in subsets:
    print("process subset", subset)
    source_dir = MOLPRO_PREFIX + '/database/sets/' + subset
    print("source_dir", source_dir, source_dir + '/' + re.sub('^.*/', '', subset))
    master_files = []
    for n in glob.glob(source_dir + '/' + '*.xml'):
        with open(n, 'r') as f:
            for line in f.readlines():
                if re.compile('<reaction').search(line):
                    master_files.append(n)
                    break
    print('master files', master_files)

    dbs = {}
    alldbs = None
    for master_file in master_files:
        # name = re.sub('/', '_', subset)
        print("subset=", subset)
        name = re.sub('/.*$', '_' + os.path.basename(master_file).replace('.xml', ''), subset.replace('Bonn/', ''))
        print("master_file", master_file, name)
        print("name", name)
        db = pymolpro.database.Database(description=name)
        root = lxml.etree.parse(master_file)
        for title in root.xpath('//db:database/@title',
                                namespaces={'db': 'http://www.molpro.net/schema/molpro-database',
                                            'xi': 'http://www.w3.org/2001/XInclude'}):
            db.description = title
        for molecule_xml in root.xpath('//xi:include/@href',
                                       namespaces={'db': 'http://www.molpro.net/schema/molpro-database',
                                                   'xi': 'http://www.w3.org/2001/XInclude'}):
            molecule_root = lxml.etree.parse(source_dir + '/' + molecule_xml)
            # print("molecule", molecule_xml, molecule_root)
            geometry = ''
            for atom in molecule_root.xpath('//cml:atom',
                                            namespaces={'molpro-output': 'http://www.molpro.net/schema/molpro-output',
                                                        'cml': 'http://www.xml-cml.org/schema'}):
                geometry += atom.get('elementType') + ' ' + atom.get('x3') + ' ' + atom.get('y3') + ' ' + atom.get(
                    'z3') + '\n'
            node = molecule_root.xpath('//molpro-output:molecule',
                                       namespaces={'molpro-output': 'http://www.molpro.net/schema/molpro-output',
                                                   'cml': 'http://www.xml-cml.org/schema'})[0]
            cmlnode = node.xpath('//cml:molecule',
                                 namespaces={'molpro-output': 'http://www.molpro.net/schema/molpro-output',
                                             'cml': 'http://www.xml-cml.org/schema'})[0]
            index = node.get('index')
            energy = float(node.get('energy')) if node.get('energy') else None
            db.add_molecule(index, geometry.strip(),
                            energy=energy,
                            InChI=node.get('InChI'),
                            SMILES=node.get('SMILES'),
                            description=node.get('title'),
                            )
            if cmlnode.get('formalCharge') and cmlnode.get('formalCharge') != '0':
                db.molecules[index]['charge'] = int(cmlnode.get('formalCharge'))
            if cmlnode.get('spinMultiplicity') and cmlnode.get('spinMultiplicity') != '1':
                db.molecules[index]['spin'] = int(cmlnode.get('spinMultiplicity')) - 1

        if name in ['benchmarks_N2minimal']:  # special case of diatomic potential curves
            for molecule in db.molecules:
                if '\n' not in db.molecules[molecule]['geometry']:
                    atom = molecule
            for molecule in db.molecules:
                if '\n' in db.molecules[molecule]['geometry']:
                    db.add_reaction(molecule, {molecule: -1, atom: 2})

        else:
            for reference_node in root.xpath('//db:database/db:reference',
                                             namespaces={'db': 'http://www.molpro.net/schema/molpro-database',
                                                         'xi': 'http://www.w3.org/2001/XInclude'}):
                if reference_node.get('doi'):
                    db.add_reference('doi', 'https://doi.org/' + reference_node.get('doi'))
            for reaction_node in root.xpath('//db:database/db:reaction',
                                            namespaces={'db': 'http://www.molpro.net/schema/molpro-database',
                                                        'xi': 'http://www.w3.org/2001/XInclude'}):
                # print('reaction')
                stoichiometry = {}
                energy = None
                for thermo_node in reaction_node.xpath('db:thermo',
                                                       namespaces={'db': 'http://www.molpro.net/schema/molpro-database',
                                                                   'xi': 'http://www.w3.org/2001/XInclude'}):
                    # print("thermo",thermo_node.text,thermo_node.get('units'))
                    energy = float(thermo_node.text) * pymolpro.database.units[thermo_node.get('units')]
                    # print('energy',energy)
                for species_node in reaction_node.xpath('db:species',
                                                        namespaces={
                                                            'db': 'http://www.molpro.net/schema/molpro-database',
                                                            'xi': 'http://www.w3.org/2001/XInclude'}):
                    species = species_node.get('index')
                    count = species_node.get('count')
                    if count is None:
                        count = 1
                    print('species', species, count)
                    if count:
                        stoichiometry[species] = int(count)
                print('stoichiometry', stoichiometry, energy)
                db.add_reaction(
                    reaction_node.get('index') if reaction_node.get('index') is not None else reaction_node.get(
                        'title'), stoichiometry, description=reaction_node.get('title'), energy=energy)
                # print("reaction energies now",db.reaction_energies)
        dbs[name] = db
        if alldbs is None:
            alldbs = db.copy()
        else:
            alldbs += db
            alldbs.add_subset(os.path.basename(master_file).replace('.xml', ''), list(db.reactions.keys()))
        print("name=", name, os.path.basename(master_file).replace('.xml', ''))
        db.dump(pymolpro.database.library_path(name))
        print(pymolpro.database.load(name))

    # global_name = re.sub('/', '_', subset)
    # print('global_name',global_name)
    # alldbs.dump(pymolpro.database.library_path(global_name))
    # print(pymolpro.database.load(global_name))

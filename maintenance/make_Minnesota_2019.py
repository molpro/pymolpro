import copy

import pymolpro
import os
import re

if not os.path.isdir('Geometries_for_Minnesota_Database_2019'):
    import requests

    resp = requests.get(
        'https://conservancy.umn.edu/bitstream/handle/11299/208752/Geometries_for_Minnesota_Database_2019.zip?sequence=6&isAllowed=y')
    with open('Geometries_for_Minnesota_Database_2019.zip', 'wb') as f: f.write(resp.content)
    import shutil

    shutil.unpack_archive('Geometries_for_Minnesota_Database_2019.zip', '.')

directory = 'Geometries_for_Minnesota_Database_2019/txt_files'
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    handle = filename[:-9]
    db = pymolpro.database.Database(description='Minnesota 2019 ' + handle)
    if handle in ['PEC4', 'MGBL193']: continue  # TODO handle these exceptional cases
    # if handle != 'S66_2.00': continue
    assert os.path.isfile(f)
    print(handle, f)
    molecule_name = ""
    with open(f, "r") as fh:
        while True:
            line = fh.readline()
            # print("first line", line)
            if not line: break
            if line.isspace(): continue
            while line[:5] != '%chk=':
                line = fh.readline()
            # print("parsed line",line)
            s = '%chk=([A-Za-z0-9_,.\\+=-]+)(\\b|\\.chk)*'
            assert (re.match(s, line))
            molecule_name = re.sub(s, r'\1', line)
            if molecule_name[-1] == '\n': molecule_name = molecule_name[:-1]
            # print("molecule_name", molecule_name[-4:])
            if molecule_name[-4:] == ".chk": molecule_name = molecule_name[:-4]
            # print("molecule_name", molecule_name)
            while True:
                line = fh.readline()
                if line[0] != '#' and not line.isspace() and line[0] != '%': break
            molecule_title = line[:-1]
            # print("molecule_title", molecule_title)
            assert fh.readline().isspace()
            line = fh.readline()
            # print("line=", line)
            pattern = r' *(-?[0-9])[, ]+([0-9])'
            assert re.match(pattern, line)
            charge = int(re.sub(pattern, r'\1', line))
            spin = int(re.sub(pattern, r'\2', line)) - 1
            charge = charge if charge != 0 else None
            spin = spin if spin != 0 else None
            geometry = ""
            while True:
                line = fh.readline()
                if line.isspace() or not line: break
                geometry += line
            # print("molecule geometry", geometry)
            db.add_molecule(molecule_name, geometry, description=molecule_title, charge=charge, spin=spin)
            while True:
                line = fh.readline()
                if not line or re.match('>+', line): break
    db.references['Geometries for Minnesota Database 2019'] = 'https://doi.org/10.13020/217y-8g32'
    db.dump(pymolpro.database.library_path('Minnesota_2019_' + handle))

    print(pymolpro.database.library('Minnesota_2019_' + handle))

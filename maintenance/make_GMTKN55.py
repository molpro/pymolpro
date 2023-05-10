import copy
import sys

import pymolpro
import os
import re
import requests
import shutil
import lxml

subsets = ['W4-11',
           'G21EA',
           'G21IP',
           'DIPCS10',
           'PA26',
           'SIE4x4',
           'ALKBDE10',
           'YBDE18',
           'AL2X6',
           'HEAVYSB11',
           'NBPRC',
           'ALK8',
           'RC21',
           'G2RC',
           'BH76RC',
           'FH51',
           'TAUT15',
           'DC13',
           'MB16-43',
           'DARC',
           'RSE43',
           'BSR36',
           'CDIE20',
           'ISO34',
           'ISOL24',
           'C60ISO',
           'PArel',
           'BH76',
           'BHPERI',
           'BHDIV10',
           'INV24',
           'BHROT27',
           'PX13',
           'WCPT18',
           'RG18',
           'ADIM6',
           'S22',
           'S66',
           'HEAVY28',
           'WATER27',
           'CARBHB12',
           'PNICO23',
           'HAL59',
           'AHB21',
           'CHB6',
           'IL16',
           'IDISP',
           'ICONF',
           'ACONF',
           'Amino20x4',
           'PCONF21',
           'MCONF',
           'SCONF',
           'UPU23',
           'BUT14DIOL', ]
# subsets = ['W4-11', 'G21EA']

directory = 'GMTKN55'


def ensure_file(path):
    global directory
    from pathlib import Path
    url_head = 'http://www.thch.uni-bonn.de/tc.old/downloads/GMTKN'
    path_ = Path(directory) / Path(path)
    if not os.path.exists(path_):
        resp = requests.get(url_head + '/' + directory + '/' + path)
        os.makedirs(os.path.realpath(os.path.join(path_, '..')), exist_ok=True)
        with open(path_, "wb") as f: f.write(resp.content)
    return path_


def html_repair(file):
    t = open(file, 'r')
    contents = t.readlines()
    t.close()
    with open(file, 'w') as t:
        for l in contents:
            l = l.replace(r'&tau;', 'tau')
            l = l.replace(r'&omega;', 'omega')
            l = l.replace(r'&aacute;', 'a')
            l = l.replace(r'&iacute;', 'i')
            l = l.replace(r'&acute;', '\'')
            l = l.replace(r'&auml;', 'ae')
            l = l.replace(r'&Rcaron;', 'R')
            l = l.replace(r'&ccaron;', 'c')
            l = re.sub(r'href=([A-Za-z0-9_-]+\.txt)>', r'href="\1">', l)
            l = l.replace(r'<br>', r'<br/>')
            l = l.replace(r'<A ', r'<a ')
            l = l.replace(r'align=right', r'align="right"')
            l = l.replace(r'Referenz<br/>', r'Referenz</p>')  # hack DIPCS10ref.html
            l = re.sub(r'(<p>.*2824-2834.<br/>)$', r'\1</p>', l)  # hack YBDE18
            l = re.sub(r'(<p>.*2012-2018.<br/>)$', r'\1</p>', l)  # hack FH51ref
            l = re.sub(r'(<p>Taken from \(this work\).<br/>)$', r'\1</p>', l)  # hack RC21ref
            l = re.sub(r'\&(downloads|lang|subsection)', r'%26\1', l)  # hack NBPRC
            l = l.replace(r'</td><td> ADIM6', r' ADIM6')  # hack ADIM6ref.html
            l = l.replace(r'</td><td> </h1>', r' </h1>')  # hack PNICO23ref.html
            l = re.sub(r'^/table>', r'</table>', l)  # hack S66ref
            l = re.sub(r'</strong>2012', r'<strong>2012', l)  # hack HAL59ref
            l = re.sub(r'Theory Comput.<em>', r'Theory Comput.</em>', l)  # hack HAL59ref
            l = re.sub(r'<td> FI pyr </td>', r'<td> FI </td><td>pyr </td>', l)  # hack HAL59ref
            l = re.sub(r'(^Taken from A L. Goerigk and S. Grimme,.*107-126.*<br/> *)$', r'<p>\1', l)  # hack PCONF21ref
            t.write(l)
    return file


kcal = 0.00159360144
for subset in subsets:
    # print("process subset",subset)
    db = pymolpro.database.Database(description='GMTKN55 ' + subset)
    db.references['GMTKN55'] = 'https://www.chemie.uni-bonn.de/grimme/de/software/gmtkn/gmtkn55'
    db.references[
        'L. Goerigk, A. Hansen. C. A. Bauer, S. Ehrlich, A. Najibi and S. Grimme in Phys. Chem. Chem. Phys., 2017'] = 'https://doi.org/10.1039/C7CP04913G'
    ensure_file(subset + 'ref.html')
    ensure_file('GMTKN.css')
    shutil.unpack_archive(ensure_file(subset + '.tar'), directory)
    for filename in os.listdir(os.path.join(directory, subset if subset != 'BH76RC' else 'BH76')):
        f = os.path.join(directory, subset, filename, 'struc.xyz')
        if not os.path.exists(f): continue
        with open(f, 'r') as fh:
            geometry = ' '.join(fh.readlines())
        try:
            with open(os.path.join(directory, subset, filename, '.UHF'), 'r') as fh:
                spin = int(' '.join(fh.readlines()))
        except:
            spin = None
        spin = spin if spin != 0 else None
        try:
            with open(os.path.join(directory, subset, filename, '.CHRG'), 'r') as fh:
                charge = int(' '.join(fh.readlines()))
        except:
            charge = None
        charge = charge if charge != 0 else None
        # print(geometry)
        db.add_molecule(filename, geometry, charge=charge, spin=spin)
    # file = ensure_file(subset + '.html')
    # html_repair(file)
    # root = lxml.etree.parse(file)
    # for root.xpath("//h3/following-sibling::ul/li/a/@href"))
    ref = lxml.etree.parse(html_repair(ensure_file(subset + 'ref.html')))
    for url in ref.xpath('//p/a/@HREF'):
        db.references['reference data'] = url
    for row in ref.xpath("//table/tr"):
        cols = row.xpath('td/text()')
        if len(cols) < 1: continue
        key = cols[0]
        ref = cols[-1]
        nvalue = (sum([1 if c.strip() else 0 for c in cols]) - 1) // 2
        offset = 1 + nvalue
        while not cols[offset].strip(): offset += 1
        # print(cols,len(cols),len(cols)//2-1,nvalue,offset)
        stoichiometry = {}
        for i in range(nvalue):
            # print("i",i)
            # print("key",str(cols[i+1]).strip())
            # print("value",cols[i+offset])
            if str(cols[i + 1]).strip():
                stoichiometry[str(cols[i + 1]).strip()] = int(cols[i + offset])
        # print("stoichiometry:",stoichiometry)
        db.add_reaction(cols[0], stoichiometry, energy=float(cols[-1]) * kcal)
    # print(db)

    if subset == 'S22':
        db.add_subset('small', ['2', '1', '8'])
        db.add_reference('Jurecka, P.; Sponer, J.; Cerny, J.; Hobza, P. Phys. Chem. Chem. Phys. 2006, 8, 1985-1993',
                         'https://doi.org/10.1039/B600027D')

    if subset == 'S66':
        db.add_reference('J. Řezáč, K. E. Riley and P. Hobza, J. Chem. Theory Comput., 2011, 7, 2427-2438',
                         'https://doi.org/10.1021/ct2002946')

    if subset in ['SIE4x4', 'BH76', 'BH76RC', 'BHPERI', 'RSE43', 'G2RC', 'O3ADD6', 'DC9', 'IDISP', 'SCONF', 'PA26']:
        db.add_reference('L. Goerigk; S. Grimme, J. Chem. Theory Comput., 2010, 6, 107-126',
                         'https://doi.org/10.1021/ct900489g')

    if subset == 'W4-11':
        db.add_reference('A. Karton, S. Daon, J. M. L. Martin and B. Ruscic, Chem. Phys. Lett., 2011, 510, 165-178',
                         'http://dx.doi.org/10.1016/j.cplett.2011.05.007')

    dbname = 'GMTKN55_' + subset

    db.dump(pymolpro.database.library_path(dbname))

    # print(pymolpro.database.library(dbname))

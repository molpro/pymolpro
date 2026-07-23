import pymolpro
import os
import re
import requests
import shutil
import lxml
import sys

subsets = ['MB08-165']
directory = 'GMTKN30'


def ensure_file(path):
    from pathlib import Path
    url_head = 'http://www.thch.uni-bonn.de/tc.old/downloads/GMTKN'
    path_ = Path(directory) / Path(path)
    if not os.path.exists(path_):
        resp = requests.get(url_head + '/' + directory + '/' + path)
        os.makedirs(os.path.realpath(os.path.join(path_, '..')), exist_ok=True)
        with open(path_, "wb") as f:
            f.write(resp.content)
    return path_


def html_repair(file):
    t = open(file, 'r')
    contents = t.readlines()
    t.close()
    with open(file, 'w') as t:
        for line in contents:
            line = line.replace(r'&tau;', 'tau')
            line = line.replace(r'&omega;', 'omega')
            line = line.replace(r'&aacute;', 'a')
            line = line.replace(r'&iacute;', 'i')
            line = line.replace(r'&acute;', '\'')
            line = line.replace(r'&auml;', 'ae')
            line = line.replace(r'&Rcaron;', 'R')
            line = line.replace(r'&ccaron;', 'c')
            line = re.sub(r'href=([A-Za-z0-9_-]+\.txt)>', r'href="\1">', line)
            line = line.replace(r'<br>', r'<br/>')
            line = line.replace(r'<A ', r'<a ')
            line = line.replace(r'align=right', r'align="right"')
            line = line.replace(r'align=center', r'align="center"')
            line = line.replace(r'</td><td> ADIM6', r' ADIM6')  # hack ADIM6ref.html
            line = re.sub(r'(^Taken from A L. Goerigk and S. Grimme,.*107-126.*<br/> *)$', r'<p>\1',
                          line)  # hack PCONF21ref
            line = line.replace(r'bz h2', r'bz </td><td> h2')
            line = line.replace(r'</td></td> <td>', r'</td>')  # hack MB08-165
            t.write(line)
    return file

def coord2xyz(fil):
    """convert Turbomole coord file to an xyz file"""
    xyz=[]
    toang=0.529177210903
    f=open(fil,'r')
    for line in f.readlines():
        data=line.split()
        if len(data)==0: continue
        if data[0][0]=='$': continue
        x,y,z=[float(x)*toang for x in data[:-1]]
        x=str(x)
        y=str(y)
        z=str(z)
        xyz.append([x,y,z,data[3]])
    f.close()
    f=open(fil+'.xyz','w')
    f.write(str(len(xyz))+'\n')
    f.write('\n')
    for i in range(len(xyz)):
        f.write("%s  %s  %s  %s\n" % (xyz[i][3].upper(),xyz[i][0],xyz[i][1],xyz[i][2]))
    f.close()

def read_xyz(fil):
    f=open(fil,'r')
    natoms=int(f.readline())
    comment=f.readline()
    geom=""
    for i in range(natoms):
        line=f.readline()
        geom+=line
    f.close()
    return natoms,geom

def read_charge_and_spin(fil):
    """read the spins of the molecules from the README file"""
    spins={}
    charges={}
    if not os.path.exists(fil):
        return charges,spins
    f=open(fil,'r')
    line=f.readline() #first line contains description
    if 'The following species are open-shell systems' in line:
        for line in f.readlines():
            data=line.split()
            name=data[0]
            spin=data[1]
            spin=spin.strip('(')
            spin=spin.strip(')')
            spin=int(spin)
            spins[name]=spin
    elif 'The following species are charged' in line:
        for line in f.readlines():
            data=line.split()
            name=data[0]
            charge=data[1]
            charge=spin.strip('(')
            charge=spin.strip(')')
            charge=int(charge)
            charges[name]=charge
    f.close()
    return charges,spins


kcal = 0.00159360144
for subset in subsets:
    # print("process subset",subset)
    db = pymolpro.database.Database(description='GMTKN30 ' + subset)
    db.references['GMTKN30'] = 'https://www.chemie.uni-bonn.de/grimme/de/software/gmtkn/gmtkn30'
    db.references[
        'L. Goerigk, and S. Grimme in Phys. J. Chem. Theory Comput., 2011'] = 'https://doi.org/10.1021/ct100466k'
    ensure_file(subset + 'ref.html')
    ensure_file('GMTKN.css')
    shutil.unpack_archive(ensure_file('strucs/'+subset + 'structures.zip'), directory)
    molecule_subset = subset if subset != 'BH76RC' else 'BH76'


    print("directory: ",directory)
    print("molecule_subset: ",molecule_subset)

    #get the charges and spins for the molecules
    charges,spins=read_charge_and_spin(os.path.join(directory, molecule_subset+'structures','README'))
    
    for filename in os.listdir(os.path.join(directory, molecule_subset+'structures')):
        if filename=='README' or '.xyz' in filename: continue        
        f = os.path.join(directory, molecule_subset+'structures', filename)
        # if subset[:4]=='BHPE': print("filename",filename,subset,f)
        if not os.path.exists(f):
            continue

        #make xyz file
        coord2xyz(f)

        #read geometry
        natoms,geom=read_xyz(f+'.xyz')

        #wavefunction info
        spin=None
        charge=None
        if filename in spins:
            spin=spins[filename]
        if filename in charges:
            charge=charges[filename]
        
        # print(geometry)
        # if subset == 'BHPERI': print('add_molecule',filename)
        db.add_molecule(filename, geometry=geom, charge=charge, spin=spin)
        
    # file = ensure_file(subset + '.html')
    # html_repair(file)
    # root = lxml.etree.parse(file)
    # for root.xpath("//h3/following-sibling::ul/li/a/@href"))
    ref = lxml.etree.parse(html_repair(ensure_file(subset + 'ref.html')))
    for url in ref.xpath('//p/a/@HREF'):
        db.references['reference data'] = url
    for row in ref.xpath("//table/tr"):
        cols = row.xpath('td/text()')
        if len(cols) < 1:
            continue
        key = cols[0]
        ref = cols[-1]
        nvalue = (sum([1 if c.strip() else 0 for c in cols]) - 1) // 2
        offset = 1 + nvalue
        while not cols[offset].strip():
            offset += 1
        # print(cols,len(cols),len(cols)//2-1,nvalue,offset)
        stoichiometry = {}
        for i in range(nvalue):
            # print("i",i)
            # print("key",str(cols[i+1]).strip())
            # print("value",cols[i+offset])
            if str(cols[i + 1]).strip():
                try:
                    key = [k for k in db.molecules if k.upper() == str(cols[i + 1]).strip().upper()][0]
                except Exception:
                    print("subset ", subset, ": cannot find ", str(cols[i + 1]).strip(), " in ", db.molecules.keys())
                stoichiometry[key] = int(cols[i + offset])
        # print("stoichiometry:",stoichiometry)
        db.add_reaction(cols[0], stoichiometry, energy=float(cols[-1]) * kcal)
    # print(db)

    if subset == 'S22':
        db.add_subset('small', ['2', '1', '8'])
        db.add_reference('Jurecka, P.; Sponer, J.; Cerny, J.; Hobza, P. Phys. Chem. Chem. Phys. 2006, 8, 1985-1993',
                         'https://doi.org/10.1039/B600027D')

    dbname = 'GMTKN30_' + subset

    db.dump(pymolpro.database.library_path(dbname))

    # print(pymolpro.database.library(dbname))

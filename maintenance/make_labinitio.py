"""
Build pymolpro databases from the datasets on https://labinitio.org/explore/chem_data_bases/
except GDB9.  Downloads go to labinitio_dl/.
"""
import csv
import os
import re
import shutil
import tarfile
import urllib.request
import zipfile

import pymolpro

root = 'labinitio_dl'
os.makedirs(root, exist_ok=True)
KCAL = 1 / 627.5094740631
KJ = 1 / 2625.499638

elements = 'H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar'.split()
Z = {e: i + 1 for i, e in enumerate(elements)}
# ground-state spin (2S) of the atoms
atom_spin = {'H': 1, 'Be': 0, 'B': 1, 'C': 2, 'N': 3, 'O': 2, 'F': 1, 'Al': 1, 'Si': 2, 'P': 3, 'S': 2, 'Cl': 1}

PAPER = 'https://drive.google.com/uc?export=download&id='
PAGE = 'https://labinitio.org/explore/chem_data_bases/'


def fetch(name):
    path = os.path.join(root, name)
    if not os.path.exists(path):
        urllib.request.urlretrieve('https://labinitio.org/databases/' + name, path)
    return path


def unpack(name):
    path = fetch(name)
    d = os.path.join(root, 'x_' + re.sub(r'\.(tar\.gz|zip)$', '', name))
    if not os.path.isdir(d):
        if name.endswith('.zip'):
            with zipfile.ZipFile(path) as z:
                z.extractall(d)
        else:
            with tarfile.open(path) as t:
                t.extractall(d)
    return d


def files_in(d, ext):
    out = {}
    for dp, dn, fn in os.walk(d):
        if '__MACOSX' in dp:
            continue
        for f in fn:
            if f.endswith(ext) and not f.startswith('.'):
                out[f[:-len(ext)]] = os.path.join(dp, f)
    return out


def read_xyz(path):
    """returns (atoms-lines, spin or None, charge or None)"""
    with open(path, encoding='utf-8', errors='replace') as fh:
        lines = fh.read().splitlines()
    n = int(lines[0].split()[0])
    header = lines[1]
    atoms = [l.strip() for l in lines[2:] if l.strip()]
    spin = charge = None
    m = re.search(r'Charge\s*=\s*(-?\d+)\s*,\s*Multiplicity\s*=\s*(\d+)', header, re.I)
    if m is None:  # some files have the charge/multiplicity line as the first atom line
        m = re.search(r'Charge\s*=\s*(-?\d+)\s*,\s*Multiplicity\s*=\s*(\d+)', ' '.join(atoms[:1]), re.I)
        if m:
            atoms = atoms[1:]
    if m is None:
        m = re.fullmatch(r'\s*(-?\d+)[ ,]+(\d+)\s*', header)
        if m is None and re.fullmatch(r'\s*(-?\d+)[ ,]+(\d+)\s*', atoms[0] if atoms else ''):
            m = re.fullmatch(r'\s*(-?\d+)[ ,]+(\d+)\s*', atoms[0])
            atoms = atoms[1:]
    if m:
        charge, spin = int(m.group(1)), int(m.group(2)) - 1
    atoms = [' '.join(a.split()[:4]) for a in atoms]
    assert len(atoms) == n, (path, n, len(atoms))
    return atoms, spin, charge


def composition(atoms):
    c = {}
    for a in atoms:
        e = a.split()[0]
        e = e[0].upper() + e[1:].lower()
        c[e] = c.get(e, 0) + 1
    return c


def add_mol(db, name, atoms, spin=None, charge=None, description=None):
    comp = composition(atoms)
    if spin is None:
        nel = sum(Z[e] * n for e, n in comp.items()) - (charge or 0)
        spin = nel % 2
    db.add_molecule(name, '\n'.join(atoms), spin=spin or None, charge=charge or None, description=description)
    return comp


def add_atoms(db, elems):
    for e in sorted(elems):
        assert e not in db.molecules, e
        db.add_molecule(e, e, energy=0.0, spin=atom_spin[e] or None)


def add_atomisation(db, name, comp, energy, description=None):
    """energy in hartree for molecule -> atoms; molecule energy is set to -energy"""
    stoi = {name: -1}
    for e, n in comp.items():
        stoi[e] = n
    db.molecule_energies[name] = -energy
    db.add_reaction(name, stoi, energy=energy, description=description)


def new_db(description, refs):
    db = pymolpro.database.Database(description=description)
    db.add_reference('Lab Initio chemical databases', PAGE)
    for k, v in refs.items():
        db.add_reference(k, v)
    return db


def save(db, key):
    pymolpro.database.library_path
    db.dump(pymolpro.database.library_path(key))
    d = pymolpro.database.load(key)
    print(key, len(d.molecules), 'molecules', len(d.reactions), 'reactions',
          os.path.getsize(pymolpro.database.library_path(key)) // 1024, 'kB')


def floats(s):
    return float(re.sub(r'\(\w\)', '', s))


# ------------------------------------------------------------------ W4-17
def parse_com(path):
    lines = open(path, encoding='utf-8', errors='replace').read().splitlines()
    i = next(k for k, l in enumerate(lines) if re.fullmatch(r'\s*-?\d+[ ,]+\d+\s*', l))
    charge, mult = [int(x) for x in re.findall(r'-?\d+', lines[i])]
    atoms = []
    for l in lines[i + 1:]:
        if not l.strip():
            break
        atoms.append(' '.join(l.split()[:4]))
    return atoms, mult - 1, charge


def w4_17_spins():
    allc = files_in(unpack('w4-17_all.zip'), '.com')
    return {k: parse_com(v) for k, v in allc.items()}


def w4_17():
    coms = w4_17_spins()
    mr = set(files_in(unpack('w4-17_mr.zip'), '.com'))
    for variant, fname, col, what in [('', 'w4-17_ref.txt', 0, 'TAEe'), ('_TAE0', 'w4-17_ref_tae0.txt', 0, 'TAE0')]:
        db = new_db('W4-17 FCI/CBS total atomisation energies (' + what + '), W4 theory or higher',
                    {'W4-17 paper': PAPER + '1ZFZpoe3GVn-7F1yPQkYMm4g5mSxyHEJv'})
        vals = {}
        for l in open(fetch(fname), encoding='utf-8', errors='replace').read().splitlines():
            m = re.match(r'^([A-Za-z0-9_\-\(\)\.,]+)\t+\s*(-?\d+\.\d+)', l)
            if m:
                vals[m.group(1)] = float(m.group(2))
        assert len(vals) == 200, len(vals)
        elems = set()
        for name, v in vals.items():
            atoms, spin, charge = coms[name]
            comp = add_mol(db, name, atoms, spin, charge)
            elems |= set(comp)
            add_atomisation(db, name, comp, v * KCAL)
        add_atoms(db, elems)
        db.add_subset('MR', sorted(k for k in vals if k in mr))
        db.add_subset('nonMR', sorted(k for k in vals if k not in mr))
        assert len(db.subsets['MR']) == 17 and len(db.subsets['nonMR']) == 183
        save(db, 'W4-17' + variant)


# ------------------------------------------------------------------ W4-11-GEOM (geometries only)
def w4_11_geom():
    geoms = files_in(unpack('w4-11-geom.tar.gz'), '.xyz')
    coms = {k.lower(): v for k, v in w4_17_spins().items()}
    db = new_db('W4-11-GEOM CCSD(T)/CBS reference geometries (no reference energies)',
                {'W4-11-GEOM paper': PAPER + '1KxMwq6j1dqnOWI_S97naZUpYoeoT_Q2A',
                 'DFT benchmark paper': PAPER + '11rnldFRSq-tCfAYVLpBs9ZMosClQgjQQ'})
    for name in sorted(geoms):
        atoms, spin, charge = read_xyz(geoms[name])
        if name.lower() in coms:
            spin = coms[name.lower()][1]
        add_mol(db, name, atoms, spin, charge)
    assert len(db.molecules) == 122, len(db.molecules)
    save(db, 'W4-11-GEOM')


# ------------------------------------------------------------------ isomerisation sets
def isomers(key, description, refs, geoms, energies, reference_of, unit):
    """energies: dict name -> relative energy; reference_of(name) -> name of the isomer it is relative to"""
    db = new_db(description, refs)
    for name in sorted(geoms):
        atoms, spin, charge = read_xyz(geoms[name])
        add_mol(db, name, atoms, spin, charge)
    for name, e in energies.items():
        ref = reference_of(name)
        if ref == name:
            continue
        db.add_reaction(name, {name: 1, ref: -1}, energy=e * unit,
                        description=name + ' relative to ' + ref)
    save(db, key)
    return db


def pah335():
    geoms = files_in(unpack('PAH335_XYZ.tar.gz'), '.xyz')
    en = {}
    for l in open(fetch('PAH335_G4MP2_Ener.txt'), encoding='utf-8', errors='replace').read().splitlines()[1:]:
        p = l.split()
        if len(p) >= 2:
            en[re.sub(r'\.log$', '', p[0])] = float(p[1])
    assert len(en) == 335 and set(en) == set(geoms)
    fam = {}
    for n, e in en.items():
        f = n.split('_')[0]
        if f not in fam or e < en[fam[f]]:
            fam[f] = n
    isomers('PAH335', 'PAH335 isomerisation energies, G4(MP2); relative to the lowest isomer of each formula',
            {'PAH335 paper': PAPER + '1ibUmqk_BWjwP1e4cQmL9SxMs5POb4a9G'}, geoms, en,
            lambda n: fam[n.split('_')[0]], KJ)


def c40():
    geoms = files_in(unpack('C40_XYZ.tar.gz'), '.xyz')
    rows = list(csv.reader(open(fetch('C40_g4mp2_ref.csv'), encoding='utf-8-sig')))[1:]
    en = {r[0]: float(r[1]) for r in rows if r[0]}
    assert set(en) == set(geoms) and len(en) == 29
    isomers('C40', 'C40 isomerisation energies (Ee), G4(MP2); relative to isomer 1',
            {'C40 paper': PAPER + '1_9vbyyfSTxDd85bbC1WemG-F9ca82aAU'}, geoms, en, lambda n: '1', KJ)


def c8h8():
    geoms = files_in(unpack('C8H8_XYZ.tar.gz'), '.xyz')
    en = {}
    started = False
    for l in open(fetch('C8H8_refs.txt'), encoding='utf-8', errors='replace').read().splitlines():
        if l.startswith('Filename'):
            started = True
            continue
        p = l.split()
        if started and len(p) == 2:
            en[p[0]] = float(p[1])
    assert set(en) == set(geoms) and len(en) == 44, (len(en), set(en) ^ set(geoms))
    isomers('C8H8', 'C8H8 isomerisation energies, CCSD(T)/CBS from W1-F12 theory; relative to isomer 1',
            {'C8H8 paper': PAPER + '1NeS72vUiXZ_sGvo2ERDX8mKN65_ppp4D',
             'Karton and Martin, Mol. Phys. 110, 2477 (2012)': 'http://dx.doi.org/10.1080/00268976.2012.698316'},
            geoms, en, lambda n: '1', KJ)


# ------------------------------------------------------------------ amino acids
def amino():
    d = unpack('Aminio_Acids_Database.tar.gz')
    geoms = files_in(d, '.xyz')
    for variant, f in [('', 'w1-f12_ref_all-e_nonrel.txt'), ('_rel', 'w1-f12_ref_all-e_rel.txt')]:
        path = [p for p in files_in(d, '.txt').items() if p[0] + '.txt' == f][0][1]
        db = new_db('Amino acids CCSD(T)/CBS total atomisation energies from W1-F12 theory, all-electron, '
                    + ('relativistic' if variant else 'nonrelativistic'),
                    {'Amino acids paper': PAPER + '1rM6kiD_pVXNiniocOl5_JtaTIYwGJjLo',
                     'Theor. Chem. Acc. 133, 1483 (2014)': 'https://doi.org/10.1007/s00214-014-1483-4'})
        vals = {}
        for l in open(path, encoding='utf-8', errors='replace').read().splitlines():
            p = l.split()
            if len(p) == 2 and p[0] in geoms:
                vals[p[0]] = float(p[1])
        assert len(vals) == 18, len(vals)
        elems = set()
        for name, v in vals.items():
            atoms, spin, charge = read_xyz(geoms[name])
            comp = add_mol(db, name, atoms, spin, charge)
            elems |= set(comp)
            add_atomisation(db, name, comp, v * KCAL)
        add_atoms(db, elems)
        save(db, 'AminoAcids' + variant)


# ------------------------------------------------------------------ WCPT27
def wcpt27():
    geoms = files_in(unpack('Water_Cat_XYZ.tar.gz'), '.xyz')
    uncat = [153.80, 196.70, 293.30, 198.80, 287.10, 339.90, 130.30, 167.80, 256.20]
    w1 = [22.60, 56.40, 158.70, 36.30, 154.10, 249.50, 20.80, 61.40, 149.20]
    w2 = [-26.70, 11.20, 99.20, -17.10, 87.50, 207.00, -19.40, 26.30, 85.50]
    # cross-check against the reference file
    txt = open(fetch('Water_Cat_Refs.txt'), encoding='utf-8', errors='replace').read()
    nums = [float(x) for x in re.findall(r'^\s*(?:\S+)?\s*\d\s+(-?\d+\.\d+)\s*$', txt, re.M)]
    assert nums == uncat + w1 + w2, nums
    db = new_db('WCPT27 water-catalysed proton-transfer barrier heights, CCSD(T)/CBS from W1 or W2.2 theory',
                {'WCPT27 paper': PAPER + '1dWd-1tmMru7AXbMOGx_QixtuhoObtkKe',
                 'Karton, O\'Reilly and Radom, J. Phys. Chem. A 116, 4211 (2012)': 'http://dx.doi.org/10.1021/jp301499y'})
    for name in sorted(geoms):
        atoms, spin, charge = read_xyz(geoms[name])
        add_mol(db, name, atoms, spin, charge)
    for suffix, vals, sub, nw in [('', uncat, 'uncat', 0), ('_1w', w1, '1H2O', 1), ('_2w', w2, '2H2O', 2)]:
        names = []
        for n in range(1, 10):
            stoi = {'ts%d%s' % (n, suffix): 1, 'prod%d' % n: -1}
            if nw:
                stoi['h2o'] = -nw
            r = '%s_%d' % (sub, n)
            db.add_reaction(r, stoi, energy=vals[n - 1] * KJ, description=sub + ' reaction %d barrier height' % n)
            names.append(r)
        db.add_subset(sub, names)
    save(db, 'WCPT27')


# ------------------------------------------------------------------ DBOC
def dboc(key, tarname, csvname, description, refs, col):
    geoms = files_in(unpack(tarname), '.xyz')
    db = new_db(description, refs)
    rows = list(csv.reader(open(fetch(csvname), encoding='utf-8-sig')))
    elems = set()
    n = 0
    for r in rows:
        if len(r) > col and r[0].endswith(('.xyz', '.out')):
            try:
                v = float(r[col])
            except ValueError:
                continue
            name = r[0].rsplit('.', 1)[0]
            atoms, spin, charge = read_xyz(geoms[name])
            comp = add_mol(db, name, atoms, spin, charge)
            elems |= set(comp)
            add_atomisation(db, name, comp, v * KJ)
            n += 1
    assert n == len(geoms), (n, len(geoms))
    add_atoms(db, elems)
    save(db, key)


# ------------------------------------------------------------------ BH28
def bh28():
    geoms = files_in(unpack('bh28.zip'), '.xyz')
    rows = list(csv.reader(open(fetch('bh28_ref.csv'), encoding='utf-8-sig')))[1:]
    refs, subsets, typ = {}, {}, None
    for row in rows:
        if len(row) < 3 or not row[1].strip().isdigit():
            continue
        typ = row[0].strip() or typ
        key = typ + '_' + row[1].strip()
        refs[key] = float(row[2])
        subsets.setdefault(typ, []).append(key)
    assert len(refs) == 28
    db = new_db('BH28 CCSDT(Q)/CBS barrier heights from W3-F12 theory',
                {'BH28 paper': PAPER + '17GoLx5_ld_UNz1fl4kUqtL3HZNxrmXOm'})
    # the xyz files carry no charge/spin; all species are neutral closed-shell by electron-count parity
    for name in sorted(geoms):
        atoms, spin, charge = read_xyz(geoms[name])
        add_mol(db, name, atoms, spin, charge)
    for key, barrier in refs.items():
        stoi = {key + '_ts': 1, key + '_min': -1}
        if key.startswith('CADBH'):
            stoi['CADBH_c2h4'] = -1
        db.add_reaction(key, stoi, energy=barrier * KJ, description=key.split('_')[0] + ' reaction %s barrier height' % key.split('_')[1])
    for typ, keys in subsets.items():
        db.add_subset(typ, keys)
    save(db, 'BH28')


if __name__ == '__main__':
    import sys
    todo = sys.argv[1:] or ['w4_17', 'w4_11_geom', 'pah335', 'c40', 'c8h8', 'amino', 'wcpt27', 'dboc', 'bh28']
    for t in todo:
        if t == 'dboc':
            dboc('DBOC200', 'XYZ_DBOC200_SI.tar.gz', 'DBOC200_SI.csv',
                 'DBOC200 diagonal Born-Oppenheimer corrections to the total atomisation energy, CCSD/CBS (best estimate)',
                 {'DBOC200 paper': PAPER + '1Tiyej6th4CWH4x41DQS313YiJD3wWRod'}, 5)
            dboc('DBOC12', 'XYZ_DBOC12_SI.tar.gz', 'DBOC12_SI.csv',
                 'DBOC12 independent test set: diagonal Born-Oppenheimer corrections to the total atomisation energy, CCSD/cc-pVDZ',
                 {'DBOC200 paper': PAPER + '1Tiyej6th4CWH4x41DQS313YiJD3wWRod'}, 4)
        else:
            globals()[t]()

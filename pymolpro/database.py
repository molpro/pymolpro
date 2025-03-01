import copy

import pymolpro
from pymolpro import resolve_geometry
import json
import os.path
import re

__all__ = ['Database', 'load', 'run', 'analyse', 'basis_extrapolate', 'units']


class Database:
    """
    Database of molecular structures and reactions. The class contains and supports the following data.

    * Descriptions of molecules: geometry, charge, spin state.
    * Descriptions of reactions between the molecules: stochiometric factors.
    * Optionally, computed data: energies of molecules and/or reactions,s
      the parameters used to obtain them,
      and references to persistent data that can be used to further query or restart.
    * Description and external URLs that serve to further define the provenance of the data.
    * Specification of recommended subsets of the reactions

    :param list molecules: Initial molecules to be added with default options to :py:meth:`add_molecule()`.
    :param list reactions: Initial reactions to be added with default options to :py:meth:`add_reaction()`.
    :param str description: Text describing the database
    :param str method: A string specifying the ansatz that was used to compute the energies.
      In the Molpro context, this will be a valid fragment of Molpro input.
    :param str basis: A string specifying the orbital basis set that was used to compute the energies.
    :param str options: A string specifying any additional options used to compute energies.
    :param str project_directory: The path of the directory where support files that were generated
      in calculating energies can be found.

    """

    def __init__(self, molecules={}, reactions={}, description=None, method=None, basis=None, options=None,
                 project_directory=None):
        self.molecules = {}  #: Dictionary of molecules
        self.reactions = {}  #: Dictionary of reactions involving the :py:data:`molecules` together with stoichiometric factors
        for key, value in molecules.items():
            self.add_molecule(key, value)
        for key, value in reactions.items():
            self.add_reaction(key, value)
        self.preamble = ""  #: any Molpro commands that should be executed before geometry specification. Typically `angstrom` could be specified if the geometry specification is in Z-matrix format with numerical values that would, by default, be interpreted as Bohr.
        self.description = "" if description is None else description  #: Text describing the database
        self.references = {}  #: A dictionary of external references to the data. The keys should be a short-form string that you want printed, eg author, year, and the values URLs that lead to the resource.
        self.subsets = {}  #: A dictionary defining named subsets of the database reactions
        self.molecule_energies = {}  #: A dictionary with molecule keys giving the molecular energy/Hartree reference values associated with the Database. The dictionary could be empty or only partly filled.
        self.reaction_energies = {}  #: A dictionary with reaction keys giving the reaction energy/Hartree reference values associated with the Database. The dictionary could be empty or only partly filled. A database might have either, both or none of :py:data:`molecule_energies` or :py:data:`reaction_energies`.
        self.method = None  #: A string specifying the ansatz used to compute the energies. In the Molpro context, this will be a valid fragment of Molpro input.
        self.basis = None  #: A string specifying the orbital basis-set used to compute the energies.
        self.options = None  #: Any additional options used to compute reference values.
        self.project_directory = None  #: A string giving the path of the directory where support files generated in calculating energies can be found.
        self.projects = {}  #: A dictionary with molecule handles pointing to filesystem project bundles for each Molpro job that has been run.
        self.failed = {}  #: Subset of :py:data:`projects` corresponding to jobs that did not complete successfully.

    def copy(self):
        db = Database()
        db.molecules = copy.deepcopy(self.molecules)
        db.reactions = copy.deepcopy(self.reactions)
        db.preamble = str(self.preamble)
        return db

    def __eq__(self, other):
        if other is None and self is not None:
            return False
        return all([getattr(self, attr) == getattr(other, attr) for attr in
                    ['molecules', 'reactions', 'molecule_energies', 'reaction_energies']])

    def __len__(self):
        return len(self.reactions)

    def __iadd__(self, other):
        for mkey, molecule in other.molecules.items():
            self.add_molecule(mkey, molecule['geometry'])
            if mkey in other.molecule_energies:
                self.molecule_energies[mkey] = other.molecule_energies[mkey]
            for key in ['description', 'InChI', 'SMILES', 'spin', 'charge']:
                if key in molecule:
                    self.molecules[mkey][key] = molecule[key]
        for mkey, reaction in other.reactions.items():
            self.add_reaction(mkey, reaction['stoichiometry'])
            if mkey in other.reaction_energies:
                self.reaction_energies[mkey] = other.reaction_energies[mkey]
            for key in ['description']:
                if key in reaction:
                    self.reactions[mkey][key] = reaction[key]
        return self

    def __add__(self, other):
        result = self.copy()
        result += other
        return result

    def add_molecule(self, name, geometry, energy=None, description=None, InChI=None,
                     SMILES=None, spin=None, charge=None, preamble=None):
        r"""
        Add a molecule to the database.  The minimal information that is stored is the geometry, but information from each of the optional arguments, if given, is also stored in the :py:data:`molecules` dictionary.

        :param str name: The key for the molecule in :py:data:`molecules`.
        :param str geometry: The geometry. Any format recognised by Molpro can be used. This includes xyz, with or without the two header lines, or Z matrix, and lines can be separated either with newline or `;`. The geometry can be specified either as a string, or a filename or url reference, in which case the contents of the reference are resolved now.
        :param float energy:  The reference value for the energy of the molecule in Hartree units
        :param str description: Descriptive text
        :param str InChI: `InChI <https://www.inchi-trust.org>`_ string describing the molecule
        :param str SMILES: `SMILES <http://opensmiles.org/opensmiles.html>`_ string describing the molecule
        :param int spin: The spin multiplicity minus one
        :param int charge: Electrical charge of molecule
        :param str preamble: Additional Molpro input to be added before geometry
        :return: The added molecule
        :rtype: dict
        """
        _name = name.strip()
        self.molecules[_name] = {
            'geometry': resolve_geometry(geometry),
        }
        self.molecules[_name]['description'] = description if description is not None else _name
        if energy is not None:
            self.molecule_energies[_name] = energy
        if spin is not None:
            self.molecules[_name]['spin'] = spin
        if charge is not None:
            self.molecules[_name]['charge'] = charge
        if InChI is not None:
            self.molecules[_name]['InChI'] = InChI
        if SMILES is not None:
            self.molecules[_name]['SMILES'] = SMILES
        if preamble is not None:
            self.molecules[_name]['preamble'] = preamble
        return self.molecules[_name]

    def add_reaction(self, name, stoichiometry, energy=None, description=None):
        r"""
        Add a reaction to the database.  The minimal information that is stored is the stoichiometry, which references existing molecules in the database, but information from each of the optional arguments, if given, is also stored in the :py:data:`reactions` dictionary.

        :param name:  The key for the reaction in :py:data:`reactions`.
        :param stoichiometry: A dictionary describing the stoichiometry of the reaction. Each key should be a key in :py:data:`molecules`, and the value is an integer giving the number of equivalents of the molecule in the reaction, with the sign convention of positive for products, negative for reactants.
        :param energy: The reference value for the energy change of the reaction in Hartree units. If not given, and if all molecules in the reaction have an energy, it will be computed.
        :param description:  Descriptive text
        :return:  The added reaction
        :rtype: dict

        """
        _name = name.strip()
        if energy:
            __reference_energy = energy
        else:
            try:
                __reference_energy = 0.0
                for reagent, stoi in stoichiometry.items():
                    __reference_energy += stoi * self.molecule_energies[reagent]
            except Exception:
                __reference_energy = None
        self.reactions[_name] = {
            'stoichiometry': Stoichiometry(stoichiometry),
        }
        if __reference_energy is not None:
            self.reaction_energies[_name] = __reference_energy
        if description is not None:
            self.reactions[_name]['description'] = description
        return self.reactions[_name]

    def add_subset(self, subset_name, subset):
        self.subsets[subset_name.strip()] = subset if type(subset) is list else [subset]
        assert all([reaction in self.reactions for reaction in self.subsets[subset_name.strip()]])

    def add_reference(self, key, url):
        self.references[key.strip()] = url.strip()

    def subset(self, subset=None, open_shell=True, max_atoms=None, max_electrons=None):
        """
        Extract a subset of this database as a new database

        :param subset: Either a key in the :py:data:`subsets` or a list of keys in :py:data:`reactions`. If not given, all the reactions in the database
        :param open_shell: Whether to include reactions involving open-shell molecules, or prune to closed-shell only
        :param max_atoms: Exclude reactions with one or more molecules with more than this number of atoms
        :param max_electrons: Exclude reactions with one or more molecules with more than this number of electrons
        :return: The subset
        :rtype: Database
        """
        subset_list = subset if type(subset) is list else self.reactions.keys() if subset is None else self.subsets[
            subset]
        if not open_shell:
            subset_list = [kr for kr in subset_list if not any(
                ['spin' in self.molecules[k] and self.molecules[k]['spin'] != 0 for k in
                 self.reactions[kr]['stoichiometry']])]
        if max_atoms:
            subset_list = [kr for kr in subset_list if not any(
                [len(elements(self.molecules[k]['geometry'])) > max_atoms for k in
                 self.reactions[kr]['stoichiometry']])]
        if max_electrons:
            subset_list = [kr for kr in subset_list if not any(
                [(electrons(self.molecules[k]['geometry'])) > max_electrons for k in
                 self.reactions[kr]['stoichiometry']])]
        db = Database(description=self.description)
        for reaction in subset_list:
            db.reactions[reaction] = self.reactions[reaction]
            if reaction in self.reaction_energies:
                db.reaction_energies[reaction] = self.reaction_energies[reaction]
        for molecule in self.molecules:
            if any([molecule in reaction['stoichiometry'] for reaction in db.reactions.values()]):
                db.molecules[molecule] = self.molecules[molecule]
                if molecule in self.molecule_energies:
                    db.molecule_energies[molecule] = self.molecule_energies[molecule]
        db.preamble = str(self.preamble)
        if subset:
            db.description += " (subset " + str(subset) + ")"
        if not open_shell:
            db.description += " (closed shell only)"
        if max_atoms:
            db.description += " (maximum number of atoms " + str(max_atoms) + ")"
        if max_electrons:
            db.description += " (maximum number of electrons " + str(max_electrons) + ")"
        return db

    def dump(self, filename=None):
        r"""
        Dump the database in json format.

        :param str filename: If specified, the destination of the result
        :return: If :py:data:`filename` is not specified, the json is returned as a string
        :rtype: str
        """
        if filename is not None:
            with open(filename, "w") as f_:
                json.dump(self, f_, default=vars)
        else:
            return json.dumps(self, default=vars)

    def load(self, source):
        r"""
        Load the database from a json dump in either the library, a file or a string

        :param str source: Source of dump
        :return: The database
        :rtype: Database
        """
        __j = None
        from json.decoder import JSONDecodeError
        try:
            for name in [source, source + '.json', library_path(source)]:
                if os.path.isfile(name):
                    with open(name, "r") as f_:
                        __j = json.load(f_)
                    break
            if not __j:
                __j = json.loads(source)
        except JSONDecodeError as e:
            raise ValueError(
                'Cannot resolve "' + source + '" as a library key, library file name, or library-defining json string') from e

        self.molecules = __j['molecules']
        self.reactions = __j['reactions']
        for reaction in self.reactions.values():
            if reaction['stoichiometry'] is not None:
                reaction['stoichiometry'] = Stoichiometry(reaction['stoichiometry'])
        if 'molecule_energies' in __j:
            self.molecule_energies = __j['molecule_energies']
        if 'reaction_energies' in __j:
            self.reaction_energies = __j['reaction_energies']
        if 'references' in __j:
            self.references = __j['references']
        if 'preamble' in __j:
            self.preamble = __j['preamble']
        if 'description' in __j:
            self.description = __j['description']
        else:
            self.description = "Molecular and reaction database"
        if 'subsets' in __j:
            self.subsets = __j['subsets']
        return self

    def calculate_reaction_energies(self, check=False):
        if len(self.failed) == 0:
            self.reaction_energies = {}
            for reaction_name, reaction in self.reactions.items():
                if not check:
                    self.reaction_energies[reaction_name] = 0.0
                    for reagent, stoichiometry in reaction['stoichiometry'].items():
                        if not self.molecule_energies[reagent]:
                            raise ValueError(
                                "Missing database molecule energy for " + reagent + ", method=" + str(
                                    self.method) + ", basis=" + self.basis + ", project directory=" + self.project_directory + str(
                                    self.molecule_energies))
                        self.reaction_energies[reaction_name] += stoichiometry * self.molecule_energies[reagent]

    def __str__(self, rst=False, geometry=True, title=None):
        result = ''
        if title:
            result = _header(title, c='-', rst=rst)
            if self.description and self.description.replace('_', ' ') != title.replace('_', ' '):
                result += self.description.replace('_', ' ') + '\n\n'
        else:
            result = _header("Database" if self.description == "" or self.description is None else self.description,
                             c='-', rst=rst) + '\n'
        if len(self.references) > 0:
            result += _header('References', rst=rst)
            for key, value in self.references.items():
                result += ('| ' if rst else '') + key + ': ' + value + '\n'
            result += '\n'
        if len(self.molecules) > 0:
            result += _header('Molecules', rst)
            for name, molecule in self.molecules.items():
                result += name + (': ' + str(molecule) if geometry else '') + (', energy = ' + str(
                    self.molecule_energies[name]) if name in self.molecule_energies else "") + '\n'
            result += '\n'
        if len(self.reactions) > 0:
            result += _header('Reactions', rst)
            for name, reaction in self.reactions.items():
                if rst:
                    result += '| '
                result += name + ': ' + str(reaction['stoichiometry']) + (
                    ' (' + reaction['description'] + ')' if 'description' in reaction else '') + (', energy = ' + str(
                    self.reaction_energies[name]) if name in self.reaction_energies else "") + '\n'
            result += '\n'
        if len(self.subsets) > 0:
            result += _header('Subsets of reactions', rst)
            for name, subset in self.subsets.items():
                result += name + ': ' + str(subset) + '\n'
            result += '\n'
        if self.method is not None and self.method != "":
            result += _header('Method', rst) + str(self.method) + '\n'
        if self.preamble is not None and self.preamble != "":
            result += _header('Preamble', rst) + str(self.preamble) + '\n'
        return result


def _header(text, rst=True, c='^'):
    result = text.replace('_', ' ')
    if rst:
        result += '\n' + ''.join([c for t in text]) + '\n'
    else:
        result += ':\n'
    return result


def load(source):
    r"""
    Construct a :py:class:`Database` from a json dump in either the library, a file or a string

    :param str source: Source of dump
    :return: The database
    :rtype: Database
    """
    return Database().load(source)


def library_path(key):
    import os.path
    return os.path.realpath(os.path.join(__file__, '..', 'share', 'database', key.replace('.json', '') + '.json'))


def library(expression=None):
    r"""
    Obtain a, possibly filtered, list of databases in the library

    :param expression: Regular expression to limit the databases found
    :return:
    :rtype list:
    """
    import os.path
    import pathlib
    result = []
    library = os.path.realpath(os.path.join(__file__, '..', 'share', 'database'))
    for file in os.listdir(library):
        if pathlib.Path(file).suffix != '.json':
            continue
        if expression and not re.match(expression, pathlib.Path(file).stem):
            continue
        result.append(pathlib.Path(file).stem)
    return result


def run(db, ansatz=None, method="hf", basis="cc-pVTZ", location=".", parallel=None, backend="local",
        clean=False, initial="", check=False, check_energy=True, **kwargs):
    r"""
    Construct and run a Molpro job for each molecule in a :py:class:`Database`,
    and compute reaction energies.

    :param Database db:  The database that defines molecules and reactions.
    :param str method: The computational method for constructed input. Anything accepted as Molpro input, including parameters and directives, can be given.  If the method needs a preceding Hartree-Fock calculation, this is prepended automatically.
    :param str basis: The orbital basis set for constructed input. Anything that can appear after `basis=` in Molpro input is accepted.
    :param str ansatz: String of the form method/basis//geometry_method/geometry_basis or method/basis which is parsed to give the same effect as the
    method and basis parameters. If geometry_method/geometry_basis is specified, the calculation will be preceded by a geometry optimisation at that level of theory.
    :param str func: This should be one of

        * `energy` for a single geometry
        * `opt` for a geometry optimisation

    :param str extrapolate: If specified, carry out basis-set extrapolation. Anything that can appear after `extrapolate,basis=` in Molpro input is accepted.
    :param str location: The filesystem directory in which projects will be constructed.
    :param int parallel: The number of simultaneous jobs to be launched. The default is the number of cores on the local machine.
    :param str backend: The sjef backend to be used for running jobs.
    :param bool clean: Whether to destroy the project bundles on successful completion. This should not normally be done, since later invocations of :py:meth:`run()` will use cached results when possible. If there are errors, this parameter is ignored.
    :param str initial: Any valid molpro input to be placed before the geometry specification.
    :param bool check: Whether to check for status of jobs instead of running them.
    :param bool check_energy: Whether to throw an exception if any job did not set the Molpro ENERGY variable
    :param kwargs: Any other options to pass to :py:meth:`project.Project.run()`, including `func`, `extrapolate`, `preamble`, `postamble`.
    :return: A new database which is a copy of :py:data:`db` but with the new results overwriting any old ones
    :rtype: Database
    """
    if parallel is None:
        from multiprocessing import cpu_count
        __parallel = cpu_count()
    else:
        __parallel = parallel
    from shutil import rmtree
    import hashlib
    from multiprocessing.dummy import Pool
    from operator import methodcaller
    from pymolpro import Project
    import os
    if type(db) is str:
        db = load(db)
    newdb = db.copy()
    newdb.project_directory = os.path.realpath(
        os.path.join(location,
                     hashlib.sha256(
                         ((str(ansatz) if ansatz is not None else str(method) + str(basis)) + str(initial) +
                          str(tuple(sorted(kwargs.items())))).encode('utf-8')).hexdigest()[-8:]))
    if not os.path.exists(newdb.project_directory):
        os.makedirs(newdb.project_directory)
    newdb.projects = {}
    for molecule_name, molecule in db.molecules.items():
        method_ = method if type(method) is str else method[1] if 'spin' in molecule and int(molecule['spin']) > 0 else \
            method[0] # TODO implement this for ansatz too
        initial_ = initial + '\n' + db.preamble
        if 'preamble' in molecule:
            initial_ += '\n' + molecule['preamble']
        if re.match('^(\n+;+ +)*$', initial_):
            initial_ = None
        try:
            newdb.projects[molecule_name] = Project(molecule_name, geometry=molecule['geometry'],
                                                    ansatz=ansatz,
                                                    method=method_,
                                                    basis=basis,
                                                    location=newdb.project_directory,
                                                    initial=initial_,
                                                    spin=molecule['spin'] if 'spin' in molecule else None,
                                                    charge=molecule['charge'] if 'charge' in molecule else None,
                                                    **kwargs)
        except Exception:
            raise FileNotFoundError(
                "pymolpro project " + molecule_name + " in directory " + newdb.project_directory + " cannot be opened and might be corrupt")
    if check:
        for k, p in newdb.projects.items():
            print("Project", k, "status:", p.status)
    else:
        with Pool(processes=__parallel) as pool:
            pool.map(methodcaller('run', backend=backend, wait=True), newdb.projects.values(), 1)

    newdb.failed = {}
    for molecule in db.molecules:
        project = newdb.projects[molecule]
        if check:
            print("checking for failure of", molecule, project.filename())
        if check:
            print("status", project.status)
        if check:
            print("no_errors", pymolpro.no_errors([project]))
        if project.status != 'completed' or not pymolpro.no_errors([project]):
            if check:
                print("failed")
            newdb.failed[molecule] = project
        else:
            if check:
                print("not failed")
    if check:
        print("after storing failed")

    newdb.molecule_energies = {}
    if check:
        print('constructing newdb and getting run results')
    for molecule_name in db.molecules:
        try:
            newdb.molecule_energies[molecule_name] = newdb.projects[molecule_name].variable('energy')
            if check_energy and newdb.molecule_energies[molecule_name] is None:
                raise ValueError('ENERGY variable is empty')
            if type(newdb.molecule_energies[molecule_name]) is type([]):
                newdb.molecule_energies[molecule_name] = newdb.molecule_energies[molecule_name][0]

        except Exception:
            if not check:
                raise ValueError(
                    "Failure to get value of ENERGY variable from " + newdb.projects[molecule_name].filename("xml"))
        if 'func' in kwargs and kwargs['func'][:3].lower() == 'opt':
            try:
                Angstrom = 1.88972612462577
                newdb.molecules[molecule_name]['geometry'] = '\n'.join(
                    [atom['elementType'] + ' ' + ' '.join([str(c / Angstrom) for c in atom['xyz']]) for atom in
                     (newdb.projects[molecule_name].geometry())])
            except Exception:
                if not check:
                    raise ValueError(
                        "Failure to get geometry from " + newdb.projects[molecule_name].filename("xml"))
            if check:
                print('got geometry', newdb.molecules[molecule_name]['geometry'])
    if check:
        print("after getting molecule_energies")
    newdb.method = method
    newdb.basis = basis
    newdb.options = sorted(kwargs.items())
    newdb.calculate_reaction_energies(check)
    if clean:
        newdb.projects = {}
        rmtree(newdb.project_directory)
        newdb.project_directory = None

    return newdb


import collections


class Units(collections.abc.Mapping):
    def __init__(self, d):
        self._d = d
        self._s = dict((k.lower(), k) for k in d)

    def __contains__(self, k):
        return k.lower() in self._s

    def __len__(self):
        return len(self._s)

    def __iter__(self):
        return iter(self._s)

    def __getitem__(self, k):
        if k is None:
            return 1.0
        if type(k) is float or type(k) is int:
            return float(k)
        return self._d[self._s[k.lower()]]

    def keys(self):
        return self._d.keys()

    def items(self):
        return self._d.items()

    def __str__(self):
        if len(self) == 0:
            return ""
        s = "Defined units:\n"
        for k, v in self.items():
            s += "1 " + k + " = " + str(v) + " a.u.\n"
        return s

    def actual_key_case(self, k):
        return self._s.get(k.lower())


units = Units({
    'kJ/mol': 1 / 2625.49963948,
    'kcal/mol': 4.184 / 2625.49963948,
    'eV': 1 / 27.211386249880,
    'meV': 1 / 27.211386249880 / 1000.0,
    'mEh': 1 / 1000.0,
    'mH': 1 / 1000.0,
    'uEh': 1 / 1000.0 / 1000.0,
    'uH': 1 / 1000.0 / 1000.0,
    'cm-1': 1 / 219474.6313632,
    'K': 1 / 315776.177845143,
})  #: dictionary of units, giving their values in atomic units


def analyse(databases, reference_database=None, unit=None, **kwargs):
    r"""
    Analyse and format the results in one or more databases

    :param list(Database)|Database databases:
    :param Database reference_database:
    :param unit: Either a string or a float specifying desired units for the output. In the case of  a string, it should match (case insensitive) one of the keys of :py:data:`units`; if an explicit value is given, it should be the value of the desired units in atomic units, ie the results will be divided by the value.
    :return: A dictionary with keys, where possible,

        * `molecule energies`
        * `molecule energy deviations`
        * `molecule statistics`
        * `molecule violin`
        * `reaction energies`
        * `reaction energy deviations`
        * `reaction statistics`
        * `reaction violin`

    The energy deviations are, for each molecule or reaction, the difference between the energies in each element of :py:data:`databases` and the value in :py:data:`reference_database`.
    The values are pandas Dataframe objects, whose columns correspond to the elements of :py:data:`databases`. In the case of the molecule or reaction energies or energy deviations, the rows correspond to the individual molecule or reaction.
    In the case of `statistics`, the rows are

        * `MSD`: the mean of the deviations
        * `STDEVD`: the standard deviation of the deviations
        * `MAD`: the mean of the absolute value of the deviations
        * `MAXD`: the maximum absolute value of the deviations
        * `RMSD`: the root mean square deviation

    :rtype: dict
    """
    import statistics
    from math import sqrt
    databases_ = list(databases) if type(databases) is list else [databases]
    output = {}
    for typ in ['reaction', 'molecule']:
        results = []
        for database in databases_:
            if len(getattr(database, typ + '_energies')) == 0:
                continue
            results.append({})
            results[-1]['method'] = database.method
            results[-1]['basis'] = database.basis
            results[-1]['unit'] = unit if unit else 'Hartree'
            if len(getattr(database, typ + '_energies')) == 0:
                continue
            results[-1][typ + ' energies'] = {key: value / units[unit] for key, value in
                                              getattr(database, typ + '_energies').items()}
            if reference_database is not None and len(getattr(reference_database, typ + '_energies')) > 0:
                results[-1][typ + ' energy deviations'] = {
                    key: value - getattr(reference_database, typ + '_energies')[key] / units[unit]
                    for
                    key, value in
                    results[-1][typ + ' energies'].items()}
                results[-1][typ + ' statistics'] = {
                    'MAD': statistics.mean([abs(v) for v in results[-1][typ + ' energy deviations'].values()]),
                    # mean absolute deviation
                    'MAXD': max([abs(v) for v in results[-1][typ + ' energy deviations'].values()]),
                    # maximum absolute deviation
                    'RMSD': sqrt(statistics.mean([v * v for v in results[-1][typ + ' energy deviations'].values()])),
                    # root mean square deviation
                    'MSD': statistics.mean(results[-1][typ + ' energy deviations'].values()),  # mean of the deviations
                    'STDEVD': statistics.stdev(results[-1][typ + ' energy deviations'].values()) if len(
                        results[-1][typ + ' energy deviations']) > 1 else 0.0,
                    # standard deviation of the deviations
                }
        output['results'] = results
        for table in [typ + ' energies', typ + ' energy deviations', typ + ' statistics']:
            if results and all([table in result for result in results]):
                output[table] = __compare_database_runs_format_table(results, table)
    if reference_database is not None:
        output['molecule violin plot'] = violin_plot(output, reactions=False,
                                                     reference_method=reference_database.method,
                                                     **kwargs)
        output['reaction violin plot'] = violin_plot(output, reactions=True, reference_method=reference_database.method,
                                                     **kwargs)
    return output


def __compare_database_runs_format_table(results, dataset):
    import numpy as np
    import pandas as pd
    output = pd.DataFrame(np.array([list(result[dataset].values()) for result in results]).transpose(),
                          index=(list(results[0][dataset].keys())))
    column_headers = [
        ((result['method'] if type(result['method']) == str else result['method'][-1]).upper() if 'method' in result and
                                                                                                  result[
                                                                                                      'method'] else "")
        + ("/" + result['basis'] if 'basis' in result and result['basis'] else "")
        for result in results
    ]
    if all(column_headers):
        output.columns = column_headers
    output.style.set_table_attributes("style='display:inline'").set_caption(dataset)
    return output


def violin_plot(analysis, reactions=True, omitted_methods=[], reference_method=None, zero_hline=True, title=None):
    r"""
    Construct a violin plot of deviations from a collection of databases

    :param dict analysis: result of a call to database.analyse().
    :param bool reactions: Whether to analyse reaction energies rather than molecule energies
    :param list omitted_methods:
    :param str reference_method:
    :param bool zero_hline: Whether to draw the y=0 x axis
    :param str title:
    :return:
    """
    # print("analysis keys",analysis.keys())
    if reference_method:
        ref_meth = reference_method
    else:
        for k in ['reaction energy deviations', 'molecule energy deviations']:
            if k in analysis:
                ref_meth = analysis[k].keys()[-1]
    if not ref_meth: return None
    # ref_meth = list(analysis['molecule energy deviations'].keys())[-1] if reference_method is None else reference_method
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError or ImportError:
        return None
    fig, pane = plt.subplots(nrows=1, ncols=1, sharey=True, figsize=(6, 6))
    deviations = 'reaction energy deviations' if reactions else 'molecule energy deviations'
    if deviations not in analysis:
        return None
    methods_pruned = [method for method in analysis[deviations] if
                      method not in omitted_methods and method != ref_meth]
    if len(methods_pruned) == 0:
        return None
    data = analysis[deviations][methods_pruned].to_numpy()
    if data.size == 0:
        return None
    pane.violinplot(data, showmeans=True, showextrema=True, vert=True, bw_method='silverman')
    pane.set_xticks(range(1, len(methods_pruned) + 1), labels=methods_pruned, rotation=-90)
    pane.set_title(title if title else analysis['results'][-1]['basis'])
    pane.set_ylabel('Error relative to ' + str(ref_meth) + ' / ' + analysis['results'][-1]['unit'])
    if zero_hline:
        plt.axhline(color='gray', linestyle=':', linewidth=0.5)
    plt.close()
    return fig


def basis_extrapolate(results, hf_results, x=None):
    assert len(hf_results) == len(hf_results)
    if x == None:
        assert len(results) > 1
        xs = []
        newresults = []
        for result in results:
            letter = re.sub(r'.*V([DTQ5-9])Z.*', r'\1', result.basis.upper())
            if len(letter) == 1:
                xs.append(int(letter.replace('D', '2').replace('T', '3').replace('Q', '4')))
                newresults.append(result)
        return basis_extrapolate(newresults,
                                 [list(hf_results)[list(results).index(newresult)] for newresult in newresults], xs)

    assert len(x) == len(results)
    assert len(x) == len(hf_results)
    if len(results) > 2:
        for first in range(len(x) - 1):
            try:
                partner = x.index(x[first] + 1)
                return basis_extrapolate(
                    [list(results)[first], list(results)[partner]],
                    [list(hf_results)[first], list(hf_results)[partner]],
                    [x[first], x[partner]]
                ) + basis_extrapolate(
                    [list(results)[i] for i in range(len(results)) if i != first],
                    [list(hf_results)[i] for i in range(len(results)) if i != first],
                    [x[i] for i in range(len(results)) if i != first]
                )
            except Exception:
                pass
        return []
    if len(x) < 2 or x[1] != x[0] + 1:
        return []
    assert len(results) == 2
    newdb = results[0].copy()
    for molecule_name in results[0].molecule_energies:
        newdb.molecule_energies[molecule_name] = __extrapolate_single(
            [result.molecule_energies[molecule_name] for result in results],
            [result.molecule_energies[molecule_name] for result in hf_results],
            x
        )
    for reaction_name in results[0].reaction_energies:
        newdb.reaction_energies[reaction_name] = __extrapolate_single(
            [result.reaction_energies[reaction_name] for result in results],
            [result.reaction_energies[reaction_name] for result in hf_results],
            x
        )
    newdb.method = results[0].method
    newdb.basis = re.sub('(.*[Vv])[dtqDTQ567]([Zz].*)', f'\\1[{str(min(x)) + str(max(x))}]\\2',
                         results[0].basis) if re.match('.*[Vv][dtqDTQ567][Zz].*',
                                                       results[0].basis) is not None else '[' + str(min(x)) + str(
        max(x)) + ']'
    return [newdb]


def __extrapolate_single(energies, hf_energies, x):
    return (hf_energies[0] if x[0] > x[1] else hf_energies[1]) + (
            x[0] * x[0] * x[0] * (energies[0] - hf_energies[0])
            - x[1] * x[1] * x[1] * (energies[1] - hf_energies[1])
    ) / (x[0] * x[0] * x[0] - x[1] * x[1] * x[1])


class Stoichiometry(dict):
    def __str__(self):
        return self.__format_reagents(-1) + " -> " + self.__format_reagents(1)

    def __format_reagents(self, direction):
        reagents = ""
        for k, v in self.items():
            if v * direction > 0:
                reagents += ' + ' + (str(v * direction) if v * direction != 1 else "") + k
        return reagents[3:] if len(reagents) > 0 else ""


def elements(geo):
    geom = geo if type(geo) is str else geo['geometry']
    lines = geom.replace(';', '\n').strip().split('\n')
    try:
        natom = len(lines) - 2 if (int(lines[0]) == len(lines) - 2) else len(lines)
    except Exception:
        natom = len(lines)
    return [line.replace(',', ' ').strip().split(' ')[0] for line in lines[-natom:]]


def electrons(geom):
    from pymolpro.project import periodic_table
    return sum([[symbol.upper() for symbol in periodic_table].index(el.upper()) + 1 for el in (elements(geom))])


import pymolpro
from pymolpro import resolve_geometry
import json

__all__ = ['Database', 'library', 'run', 'compare']


class Database:
    """
    Database of molecular structures and reactions

    :param list molecules: Initial molecules to be added with default options to :py:meth:`add_molecule()`.
    :param list reactions: Initial reactions to be added with default options to :py:meth:`add_reaction()`.
    :param str description: Text describing the database

    """

    def __init__(self, molecules={}, reactions={}, description="None"):
        self.molecules = {}  #: Dictionary of molecules
        self.reactions = {}  #: Dictionary of reactions involving the :py:data:`molecules` together with stoichiometric factors
        for key, value in molecules.items():
            self.add_molecule(key, value)
        for key, value in reactions.items():
            self.add_reaction(key, value)
        self.preamble = ""  #: any Molpro commands that should be executed before geometry specification. Typically `angstrom` could be specified if the geometry specification is in Z-matrix format with numerical values that would, by default, be interpreted as Bohr.
        self.description = "" if description is None else description  #: Text describing the database
        self.references = {}  #: A dictionary of external references to the data. The keys should be a short-form string that you want printed, eg author, year, and the values URLs that lead to the resource.

    def __len__(self):
        return len(self.reactions)

    def add_molecule(self, name, geometry, fragment_lengths=None, reference_energy=None, description=None, InChI=None,
                     SMILES=None, spin=None, charge=None):
        r"""
        Add a molecule to the database.  The minimal information that is stored is the geometry, but information from each of the optional arguments, if given, is also stored in the :py:data:`molecules` dictionary.

        :param str name: The key for the molecule in :py:data:`molecules`.
        :param str geometry: The geometry. Any format recognised by Molpro can be used. This includes xyz, with or without the two header lines, or Z matrix, and lines can be separated either with newline or `;`. The geometry can be specified either as a string, or a filename or url reference, in which case the contents of the reference are resolved now.
        :param list fragment_lengths:  For a molecule that is to be considered as a supramolecular complex, the lengths of each of the fragments. The last value can be omitted.
        :param float reference_energy:  The reference value for the energy of the molecule in Hartree units
        :param str description: Descriptive text
        :param str InChI: `InChI <https://www.inchi-trust.org>`_ string describing the molecule
        :param str SMILES: `SMILES <http://opensmiles.org/opensmiles.html>`_ string describing the molecule
        :param int spin: The spin multiplicity minus one
        :param int charge: Electrical charge of molecule
        :return: The added molecule
        :rtype: dict
        """
        self.molecules[name] = {
            'geometry': resolve_geometry(geometry),
        }
        self.molecules[name]['description'] = description if description is not None else name
        if reference_energy is not None: self.molecules[name]['reference energy'] = reference_energy
        if fragment_lengths is not None: self.molecules[name]['fragment lengths'] = fragment_lengths
        if spin is not None: self.molecules[name]['spin'] = spin
        if charge is not None: self.molecules[name]['charge'] = charge
        if InChI is not None: self.molecules[name]['InChI'] = InChI
        if SMILES is not None: self.molecules[name]['SMILES'] = SMILES
        return self.molecules[name]

    def add_reaction(self, name, stoichiometry, reference_energy=None, description=None):
        r"""
        Add a reaction to the database.  The minimal information that is stored is the stoichiometry, which references existing molecules in the database, but information from each of the optional arguments, if given, is also stored in the :py:data:`reactions` dictionary.

        :param name:  The key for the reaction in :py:data:`reactions`.
        :param stoichiometry: A dictionary describing the stoichiometry of the reaction. Each key should be a key in :py:data:`molecules`, and the value is an integer giving the number of equivalents of the molecule in the reaction, with the sign convention of positive for products, negative for reactants.
        :param reference_energy: The reference value for the energy change of the reaction in Hartree units. If not given, and if all molecules in the reaction have a reference energy, it will be computed.
        :param description:  Descriptive text
        :return:  The added reaction
        :rtype: dict

        """
        if reference_energy:
            __reference_energy = reference_energy
        else:
            try:
                __reference_energy = 0.0
                for reagent, stoi in stoichiometry.items():
                    __reference_energy += stoi * self.molecules[reagent]['reference energy']
            except:
                __reference_energy = None
        self.reactions[name] = {
            'stoichiometry': stoichiometry,
            'reference energy': __reference_energy,
        }
        if __reference_energy is not None: self.reactions[name]['reference energy'] = __reference_energy
        if description is not None: self.reactions[name]['description'] = description
        return self.reactions[name]

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

    def load(self, filename=None, string=""):
        r"""
        Load the database from a json dump in either a file or a string

        :param str filename: Source of dump
        :param str string: Alternate source of dump if :py:data:`filename` is not given
        :return: The database
        :rtype: Database
        """
        if filename is not None:
            with open(filename, "r") as f_:
                __j = json.load(f_)
        else:
            __j = json.loads(string)
        self.molecules = __j['molecules']
        self.reactions = __j['reactions']
        if 'references' in __j:
            self.references = __j['references']
        if 'preamble' in __j:
            self.preamble = __j['preamble']
        if 'description' in __j:
            self.description = __j['description']
        else:
            self.description = "Molecular and reaction database"
        return self

    def reference_results(self):
        r"""
        The reference values stored in the database.

        :return: A dictionary containing, where possible, entries `molecule energies` and `reaction energies` constructed from the `reference energy` fields of each :py:data:`molecules` and :py:data:`reactions` entry.
        :rtype: dict:
        """
        __results = {}
        if all(['reference energy' in molecule for molecule in self.molecules.values()]):
            __results["molecule energies"] = {key: value['reference energy'] for key, value in self.molecules.items()}
        if all(['reference energy' in reaction for reaction in self.reactions.values()]):
            __results["reaction energies"] = {key: value['reference energy'] for key, value in self.reactions.items()}
        return __results

    def __str__(self):
        result = "Database" if self.description == "" or self.description is None else self.description
        if len(self.references) > 0:
            result += '\n\nReferences:\n' + str(self.references)
        result += '\n\nMolecules:\n' + str(self.molecules)
        result += '\n\nReactions:\n' + str(self.reactions)
        if self.preamble is not None and self.preamble != "":
            result += '\n\nPreamble:\n' + str(self.preamble)
        return result


def library(key):
    r"""
    Construct a :py:class:`Database` from the library
    :param str key: File name, without trailing .json, of the library database in /share/database
    :return: The database
    :rtype: Database
    """
    import os.path
    db = Database()
    db.load(os.path.realpath(os.path.join(__file__, '..', '..', 'share', 'database', key + '.json')))
    return db


def run(db, method="hf", basis="cc-pVTZ", location=".", parallel=None, backend="local",
        clean=False, initial="", **kwargs):
    r"""
    Construct and run a Molpro job for each molecule in a :py:class:`Database`,
    and compute reaction energies.

    :param Database db:  The database that defines molecules and reacitions
    :param str method: The computational method for constructed input. Anything accepted as Molpro input, including parameters and directives, can be given.  If the method needs a preceding Hartree-Fock calculation, this is prepended automatically.
    :param str basis: The orbital basis set for constructed input. Anything that can appear after `basis=` in Molpro input is accepted.
    :param str func: This should be one of

        * `energy` for a single geometry
        * `opt` for a geometry optimisation

    :param str extrapolate: If specified, carry out basis-set extrapolation. Anything that can appear after `extrapolate,basis=` in Molpro input is accepted.
    :param str location: The filesystem directory in which projects will be constructed.
    :param int parallel: The number of simultaneous jobs to be launched. The default is the number of cores on the local machine.
    :param str backend: The sjef backend to be used for running jobs.
    :param bool clean: Whether to destroy the project bundles on successful completion. This should not normally be done, since later invocations of :py:meth:`run()` will use cached results when possible. If there are errors, this parameter is ignored.
    :param str initial: Any valid molpro input to be placed before the geometry specification.
    :param kwargs: Any other options to pass to :py:meth:`project.Project.run()`, including `func`, `extrapolate`, `preamble`, `postamble`, `initial`.
    :return: A dictionary containing the results, containing the following.

        * `molecule energies` A dictionary with molecule handles pointing to the contents of the Molpro `energy` variable, which could be a scalar or a list.
        * `reaction energies` A dictionary with molecule handles pointing to the evaluated reaction energy changes.
        * `project directory`
        * `method`
        * `basis`
        * `options` Extra options that were used in constructing the input.
        * `failed` Normally empty or absent, but in the case of execution errors a dictionary with molecule handles pointing to filesystem project bundles for each job that failed.
    :rtype: dict
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
    if type(db) == str: db = library(db)
    project_dir_ = os.path.realpath(
        os.path.join(location, method.upper() + "_" + basis + "_" + hashlib.sha256(
            str(tuple(sorted(kwargs.items()))).encode('utf-8')).hexdigest()[-8:]))
    if not os.path.exists(project_dir_):
        os.makedirs(project_dir_)
    projects = {}
    for molecule_name, molecule in db.molecules.items():
        projects[molecule_name] = Project(molecule_name, geometry=molecule['geometry'], method=method, basis=basis,
                                          location=project_dir_, initial=initial + "; " + db.preamble,
                                          **kwargs)
    with Pool(processes=__parallel) as pool:
        pool.map(methodcaller('run', backend=backend, wait=True), projects.values(), 1)

    failed_jobs = {}
    for molecule in db.molecules:
        project = projects[molecule]
        if project.status != 'completed' or not pymolpro.no_errors([project]):
            failed_jobs[molecule] = project

    molecule_energies = {}
    for molecule_name in db.molecules:
        molecule_energies[molecule_name] = projects[molecule_name].variable('energy')
    result = {
        "project directory": project_dir_,
        "method": method,
        "basis": basis,
        "options": sorted(kwargs.items()),
        "molecule energies": molecule_energies,
    }
    if len(failed_jobs) > 0:
        result['failed'] = failed_jobs
        return result

    result['reaction energies'] = {}
    for reaction_name, reaction in db.reactions.items():
        result['reaction energies'][reaction_name] = 0.0
        for reagent, stoichiometry in reaction['stoichiometry'].items():
            result['reaction energies'][reaction_name] += stoichiometry * molecule_energies[reagent]

    if clean: rmtree(project_dir_)

    return result


def compare(results, reference_result, reactions=False, molecules=False):
    import statistics
    results_ = list(results) if type(results) == list else [results]
    output = {}
    for typ in ['reaction', 'molecule']:
        for result in results_:
            result[typ + ' energy errors'] = {key: value - reference_result[typ + ' energies'][key] for key, value in
                                              result[typ + ' energies'].items()}
            result[typ + ' statistics'] = {
                'mean': statistics.mean(result[typ + ' energy errors'].values()),
                'stdev': statistics.stdev(result[typ + ' energy errors'].values()),
                'meanabs': statistics.mean([abs(v) for v in result[typ + ' energy errors'].values()]),
                'maxabs': max([abs(v) for v in result[typ + ' energy errors'].values()]),
            }
        for table in [typ + ' energies', typ + ' energy errors', typ + ' statistics']:
            output[table] = __compare_database_runs_format_table(results_, table)
    return output


def __compare_database_runs_format_table(results, dataset):
    import numpy as np
    import re
    import pandas as pd
    table = [list(result[dataset].values()) for result in results]
    row_labels = list(results[0][dataset].keys())
    output = pd.DataFrame(np.array(table).transpose(), index=row_labels)
    output.columns = pd.MultiIndex.from_arrays([
        [result['method'].upper() for result in results],
        [result['basis'] for result in results],
        # [re.compile('^.*_').sub('', (result['project directory'] if 'project_directory' in list(result.keys()) else '')) for result in results],
    ])
    output.style.set_table_attributes("style='display:inline'").set_caption(dataset)
    return output


def basis_extrapolate(results, hf_results, x):
    import re
    molecule_energies = {}
    for molecule_name in results[0]['molecule energies']:
        molecule_energies[molecule_name] = __extrapolate_single(
            [result['molecule energies'][molecule_name] for result in results],
            [result['molecule energies'][molecule_name] for result in hf_results],
            x
        )
    reaction_energies = {}
    for reaction_name in results[0]['reaction energies']:
        reaction_energies[reaction_name] = __extrapolate_single(
            [result['reaction energies'][reaction_name] for result in results],
            [result['reaction energies'][reaction_name] for result in hf_results],
            x
        )
    return {
        "method": results[0]['method'],
        "basis": re.sub('(.*[Vv])[DTQ567]([Zz])', f'\\1[{str(min(x)) + str(max(x))}]\\2',
                        results[0]['basis']) if re.match('.*[Vv][DTQ567][Zz]',
                                                         results[0]['basis']) is not None else '[' + str(min(x)) + str(
            max(x)) + ']',
        "molecule energies": molecule_energies,
        "reaction energies": reaction_energies,
    }


def __extrapolate_single(energies, hf_energies, x):
    return (hf_energies[0] if x[0] > x[1] else hf_energies[1]) + (
            x[0] * x[0] * x[0] * (energies[0] - hf_energies[0])
            - x[1] * x[1] * x[1] * (energies[1] - hf_energies[1])
    ) / (x[0] * x[0] * x[0] - x[1] * x[1] * x[1])

from pymolpro import resolve_geometry
import json

__all__ = ['Database', 'library', 'run', 'compare']


class Database:
    """
    Database of molecular structures and reactions
    """

    def __init__(self, molecules={}, reactions={}):
        self.molecules = {}
        self.reactions = {}
        for key, value in molecules.items():
            self.add_molecule(key, value)
        for key, value in reactions.items():
            self.add_reaction(key, value)
        self.preamble = ""

    def __len__(self):
        return len(self.reactions)

    def add_molecule(self, name, geometry, fragment_lengths=None, reference_energy=None, description=None, InChi=None,
                     SMILES=None):
        self.molecules[name] = {
            'geometry': resolve_geometry(geometry),
        }
        self.molecules[name]['description'] = description if description is not None else name
        if reference_energy is not None: self.molecules[name]['reference energy'] = reference_energy
        if fragment_lengths is not None: self.molecules[name]['fragment lengths'] = fragment_lengths
        if InChi is not None: self.molecules[name]['InChi'] = InChi
        if SMILES is not None: self.molecules[name]['SMILES'] = SMILES

    def add_reaction(self, name, stoichiometry, reference_energy=None, description=None):
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

    def dump(self, filename=None):
        if filename is not None:
            with open(filename, "w") as f_:
                json.dump(self, f_, default=vars)
        else:
            return json.dumps(self, default=vars)

    def load(self, filename=None, string=""):
        if filename is not None:
            with open(filename, "r") as f_:
                __j = json.load(f_)
        else:
            __j = json.loads(string)
        self.molecules = __j['molecules']
        self.reactions = __j['reactions']
        if 'preamble' in __j:
            self.preamble = __j['preamble']

    def reference_results(self):
        __results = {}
        if all(['reference energy' in molecule for molecule in self.molecules.values()]):
            __results["molecule energies"] = {key: value['reference energy'] for key, value in self.molecules.items()}
        if all(['reference energy' in reaction for reaction in self.reactions.values()]):
            __results["reaction energies"] = {key: value['reference energy'] for key, value in self.reactions.items()}
        return __results


def library(key):
    import os.path
    db = Database()
    db.load(os.path.realpath(os.path.join(__file__, '..', '..', 'share', 'database', key + '.json')))
    return db


from multiprocessing import cpu_count


def run(db, method="hf", basis="cc-pVTZ", location=".", parallel=cpu_count(), backend="local",
        clean=False, initial="", **kwargs):
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
    with Pool(processes=parallel) as pool:
        pool.map(methodcaller('run', backend=backend, wait=True), projects.values(), 1)

    molecule_energies = {}
    for molecule_name in db.molecules:
        molecule_energies[molecule_name] = projects[molecule_name].variable('energy')
    reaction_energies = {}
    for reaction_name, reaction in db.reactions.items():
        reaction_energies[reaction_name] = 0.0
        for reagent, stoichiometry in reaction['stoichiometry'].items():
            reaction_energies[reaction_name] += stoichiometry * molecule_energies[reagent]

    if clean: rmtree(project_dir_)

    return {
        "project directory": project_dir_,
        "method": method,
        "basis": basis,
        "options": sorted(kwargs.items()),
        "molecule energies": molecule_energies,
        "reaction energies": reaction_energies,
    }


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
        "basis": '[' + str(min(x)) + str(max(x)) + ']',
        "molecule energies": molecule_energies,
        "reaction energies": reaction_energies,
    }


def __extrapolate_single(energies, hf_energies, x):
    return (hf_energies[0] if x[0] > x[1] else hf_energies[1]) + (
            x[0] * x[0] * x[0] * (energies[0] - hf_energies[0])
            - x[1] * x[1] * x[1] * (energies[1] - hf_energies[1])
    ) / (x[0] * x[0] * x[0] - x[1] * x[1] * x[1])

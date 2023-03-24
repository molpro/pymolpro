from pymolpro import resolve_geometry
import json
import os.path


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
        pass

    def add_molecule(self, name, geometry, fragment_lengths=[], description=None):
        self.molecules[name] = {
            'description': description,
            'geometry': resolve_geometry(geometry),
            'fragment_lengths': fragment_lengths,
        }

    def add_reaction(self, name, stoichiometry, reference_energy=None, description=None):
        self.reactions[name] = {
            'description': description,
            'stoichiometry': stoichiometry,
            'reference_energy': reference_energy,
        }

    def dump(self, filename=None):
        if filename is not None:
            with open(filename, "w") as f_:
                json.dump(self, f_, default=vars)
        else:
            return json.dumps(self, default=vars)

    def load(self, filename=None, string=""):
        if filename is not None:
            with open(filename, "r") as f_:
                j_ = json.load(f_)
        else:
            j_ = json.loads(string)
        self.molecules = j_['molecules']
        self.reactions = j_['reactions']


def library_database(key):
    db = Database()
    db.load(os.path.realpath(os.path.join(__file__, '..', '..', 'share', 'database', key + '.json')))
    return db

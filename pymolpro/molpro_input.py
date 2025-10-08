#!env python
from __future__ import annotations
import logging
from pymolpro.defbas import periodic_table
import os
import pathlib
import re
from collections import UserDict, OrderedDict
from copy import deepcopy
import json
from dataclasses import dataclass

import jsonschema
import pymolpro

with open((pathlib.Path(__file__).parent / 'molpro_input.json').as_posix(), 'r') as f:
    schema = json.load(f)

_logger = logging.getLogger(__name__)
_symmetry_commands = {
    'automatic': '',
    'none': 'symmetry,nosym'
}
assert (set(_symmetry_commands.keys()).issubset(set(schema['properties']['symmetry']['enum'])))


def symmetry_commands():
    return _symmetry_commands


_symmetry_command_aliases = {
    'nosym': 'symmetry,nosym',
}

_hamiltonians = {
    'AE': {'text': 'All Electron', 'basis_string': ''},
    'PP': {'text': 'Pseudopotential', 'basis_string': '-PP'},
    'DK': {'text': 'Douglas-Kroll-Hess', 'basis_string': '-DK'},
    'DK3': {'text': 'Douglas-Kroll-Hess 3', 'basis_string': '-DK3'},
}
assert (set(_hamiltonians.keys()).issubset(set(schema['properties']['hamiltonian']['enum'])))


def hamiltonians():
    return _hamiltonians


_local_orbital_types = {
    'ibo': {'text': 'Intrinsic Bond', 'command': 'ibba'},
    'pipek': {'text': 'Pipek-Mezey', 'command': 'locali,pipek'},
    'nbo': {'text': 'NBO', 'command': 'nbo'},
    'boys': {'text': 'Boys', 'command': 'locali'},
}


def local_orbital_types():
    r"""
    Return the supported local orbital types.
    """
    return _local_orbital_types


_parameter_commands = {
    'parameter': 'gparam',
    'threshold': 'gthresh',
    'print': 'gprint',
}

_job_types = {
    'SP': 'Single geometry',
    'OPT': 'Geometry optimization',
    'FREQ': 'Hessian',
    'OPT+FREQ': 'Geometry and hessian',
    'INTERACT': 'Non-covalent complex',
}
assert (set(_job_types.keys()).issubset(set([a['const'] for a in schema['properties']['job_type']['anyOf']])))


def job_types():
    r"""
    Return the supported job types.
    rtype: dict
    """
    return _job_types


_job_type_aliases = {
    '{optg}': 'optg',
    '{freq}': 'frequencies',
    'freq': 'frequencies',
}
_orientation_options = {
    'mass': 'mass',
    'charge': 'charge',
    'none': 'noorient'
}
assert (set(_orientation_options.keys()).issubset(set(schema['properties']['orientation']['enum'])))


def orientation_options():
    return _orientation_options


properties = {
    'Dipole moment': 'gexpec,dm',
    'Quadrupole moment': 'gexpec,qm',
    'Second moment': 'gexpec,sm',
    'Kinetic energy': 'gexpec,ekin',
    'Cowan-Griffin': 'gexpec,rel',
    'Mass-velocity': 'gexpec,massv',
    'Darwin': 'gexpec,darw',
}

_initial_orbital_methods = ['HF', 'KS']

_supported_methods = None
_procedures_registry = None


def supported_methods():
    r"""
    Returns a list of supported methods.
    """
    from pymolpro.registry import procedures_registry
    global _supported_methods
    if _supported_methods is not None: return _supported_methods

    _procedures_registry = procedures_registry()
    _supported_methods = []
    for keyfound in _procedures_registry.keys():
        if _procedures_registry[keyfound]['class'] == 'PROG':
            _supported_methods.append(_procedures_registry[keyfound]['name'])

    return _supported_methods


def procedures_registry():
    r"""
    Returns a dictionary with procedure names as keys and procedures as values.
    """
    global _procedures_registry
    if _procedures_registry is not None: return _procedures_registry
    try:
        _procedures_registry = procedures_registry()
        if not _procedures_registry:
            raise ValueError
    except Exception as e:
        _procedures_registry = {}
    return _procedures_registry


@dataclass
class JobStep:
    command: str
    options: list
    directives: list

    def __init__(self, card: str):
        self.load(card)

    def load(self, card: str):
        lines = card.replace(';', '\n').split('\n')
        fields = lines[0].split(',')
        self.command = fields[0]
        # self.options = OrderedDict()
        self.options = fields[1:]
        self.directives = lines[1:]

    def dump(self, braces=True):
        card = ('{' if braces else '') + self.command
        for option in self.options:
            card += ',' + option
        for directive in self.directives:
            card += '\n' + directive
        if braces:
            card += '}'
        return card


_default_job_type_commands = {k: v['default'] for k, v in schema['properties']['job_type_commands']['items'].items()}


def _split_quote_protected_string(string, delimiter):
    result = ['']
    quoted = False
    for s in string:
        if not quoted and s == delimiter:
            result.append('')
        else:
            if quoted and s in '\'"':
                quoted = False
            elif not quoted and s in '\'"':
                quoted = True
            result[-1] = result[-1] + s
    return result


def _load_permissive_json(string):
    if type(string) != str:
        return string
    string_with_quoted_parts_removed = re.sub(r"'[^']*'", '', re.sub(r'"[^"]*"', '', string))
    if ';' in string_with_quoted_parts_removed or '\n' in string_with_quoted_parts_removed:
        return string
    _string = string.strip()
    if _string[0] == '{' and string[-1] == '}': _string = _string[1:-1]
    try:
        result = json.loads('{' + _string + '}')
        jsonschema.validate(result, schema)
    except json.decoder.JSONDecodeError:
        result = json.loads(_convert_keyval_to_json(_string))
        jsonschema.validate(result, schema)
    return result


def _convert_keyval_to_json(string):
    split = {}
    result = ''
    for keyval in _split_quote_protected_string(string, ','):
        key, value = re.split(r'[=:]', keyval, 1)
        key = re.sub(r'["\']', '', key).strip()
        value = re.sub(r'["\']', '', value).strip()
        if key == 'basis' and 'default' not in value:
            value = 'default=' + value
        if key == 'basis':
            value_split = value.split(',')
            value = '{'
            for index, field in enumerate(value_split):
                field_split = field.split('=')
                if index == 1:
                    value = value + '"elements":{'
                value = value + '"' + field_split[0] + '":"' + (field_split[1] if len(field_split) > 1 else '') + '",'
                if index == len(value_split) - 1 and index != 0:
                    value = value[:-1] + '},'
            value = value[:-1] + '}'
        elif value == 'True':
            value = True
        elif value == 'False':
            value = False
        else:
            value = '"' + value + '"'
        result += (',' if result else '{') + '"' + key + '":' + value
    return result + '}'


class InputSpecification(UserDict):
    """ A declarative specification of a Molpro input. Almost any simple input - no loops or logic - that requests a single quantum-mechanical method, possibly with geometry optimisation using a different method, can be represented. The class can be instantiated from a procedural Molpro input, or from a JSON string or dictionary.    If both an input and a specification are provided, the specification is ignored.

    The class is a dictionary with some extra methods, with all the data contained in the dictionary, which is constrained to conform with the molpro_input JSON schema, https://www.molpro.net/schema/molpro_input.json ."""

    _hartree_fock_methods = ['RHF', 'RKS', 'UHF', 'UKS', 'LDF-RHF', 'LDF-UHF']

    @classmethod
    def default_instance(cls) -> 'InputSpecification':
        """Factory method for creating a default instance generated from the schema."""
        self = cls()
        for k, v in schema['properties'].items():
            if 'default' in v:
                self[k] = v['default']
            elif k == 'basis':
                self[k] = {'default': v['properties']['default']['default']}
            elif k == 'job_type_commands':
                self[k] = {k2: v2['default'] for k2, v2 in v['items'].items()}
        self.validate()
        return self

    @property
    def with_defaults(self) -> InputSpecification:
        """A copy of this specification with all defaults filled in from the schema."""
        result = self.default_instance()
        for k, v in self.items():
            if k in result and type(result[k]) is str and type(v) is list and len(v) == 1:
                result[k] = v[0]
            elif k in result and type(result[k]) is dict and type(v) is dict:
                for k2, v2 in v.items():
                    result[k][k2] = v2
            else:
                result[k] = v
        if 'spin' not in self:
            result['spin'] = (self.open_shell_electrons) % 2 - 2
        for key in ['method', 'basis']:
            if 'geometry_' + key not in result and key in result:
                result['geometry_' + key] = result[key]
        return result

    @property
    def without_defaults(self) -> InputSpecification:
        """A copy of this specification with any values equal to the schema default removed."""
        result = type(self)()
        di = self.default_instance()
        for k, v in self.items():
            if k in di and type(di[k]) is str and type(v) is list and len(v) == 1:
                if v[0] != di[k]:
                    result[k] = v[0]
            else:
                if k not in di or v != di[k]:
                    result[k] = v
        return result

    def _ensure_orbital_method(self):
        """
        Ensure that the 'method' field starts with a valid orbital method.
        If not present, set to default. If invalid, prepend 'rhf' or 'uhf'.
        """
        if 'method' not in self or not self['method']:
            self['method'] = self.with_defaults['method']
            return
        methods = self['method'].split(';') if isinstance(self['method'], str) else self['method']
        if not methods:
            self['method'] = self.with_defaults['method']
            return
        if methods[0].lower() not in ['hf', 'rhf', 'uhf', 'ks', 'rks', 'uks', 'avas']:
            self['method'] = ['hf' if methods[0].lower()[0] != 'u' else 'uhf'] + methods

    @property
    def job_steps(self) -> list[JobStep]:
        r"""
        Returns the job steps in this input
        """
        job_steps = []
        defaulted = self.with_defaults
        if 'method' in defaulted:
            if type(defaulted['method']) is list:
                for method_step in defaulted['method']:
                    job_steps.append(JobStep(method_step))
            else:
                job_steps.append(JobStep(defaulted['method']))
        for step in defaulted['job_type_commands'][defaulted['job_type']]:
            job_steps.append(JobStep(step))
        return job_steps

    def set_job_step(self, job_step: JobStep, index: int):
        if index < len(self['method']):
            self['method'][index] = job_step.dump(braces=False)
        else:
            if 'job_type_commands' not in self:
                self['job_type_commands'] = []
            for k in range(len(self['job_type_commands']), index):
                self['job_type_commands'].append(_default_job_type_commands[self['job_type']][k])
            if index < len(self['job_type_commands']):
                self['job_type_commands'][index] = job_step.dump(braces=False)
            else:
                self['job_type_commands'].append(job_step.dump(braces=False))

    def __init__(self, input:str=None, allowed_methods:list[str]=[], debug:bool=False, specification:dict=None, directory:str=None):
        """

        :param input: Either a text string or a file name containing Molpro procedural input.
        :param allowed_methods: A list of allowed methods to be appended to the defaults.
        :param specification: Initial data to be included in the specification.
        :param directory: The directory where any auxiliary files are located. If not specified, the current working directory will be used.
        """
        super(InputSpecification, self).__init__()
        self.allowed_methods = list(set(allowed_methods).union(set(supported_methods())))
        self.directory = directory
        self.procname = 'ansatz'
        # print('self.allowed_methods',self.allowed_methods)
        # print(specification, type(specification))
        if specification is not None:
            if type(specification) is str:
                _specification = _load_permissive_json(specification)
            else:
                _specification = specification
            # print('Specification:', _specification)
            for k in _specification:
                self[k] = deepcopy(_specification[k])
            self._ensure_orbital_method()
            # print('InputSpecification() input=',input,'specification=',specification, '_specification=',_specification)
        if input is not None:
            # self.parse(input, debug=debug)
            try:
                self.parse(input, debug=debug)
            except Exception as e:
                print('Warning: InputSpecification.parse() has thrown an exception', e,
                      '\nIf unexpected, please report, with a copy of the input, at https://github.com/molpro/iMolpro/issues/new')
                self.clear()
        if False and self.data:
            for field in ['hamiltonian']:
                if field not in self:
                    self[field] = self.with_defaults[field]

    def validate(self):
        """ Validate the specification according to the schema."""
        jsonschema.validate(instance=json.loads(json.dumps(dict(self))), schema=schema)
        pass

    @property
    def ansatz(self) -> str:
        """A string of the form method/basis, or method/basis//geometry_method/geometry_basis summarising the ansatz where possible, otherwise an empty string."""
        _method = self.with_defaults['method']
        _basis = self.with_defaults['basis']['default']
        if 'elements' in self.with_defaults['basis']:
            return ''
        if type(_method) is list:
            _method = _method[-1]
        # if _method[:3].lower() == 'df-':
        #     _method = _method[3:]
        if _method[:3].lower() != 'df-' and self.with_defaults['density_fitting']:
            _method = 'DF-' + _method
        _method = re.sub('[ru]?ks,', '', _method, flags=re.IGNORECASE)
        result = _method.upper() + '/' + _basis
        _geometry_method = self['geometry_method'] if 'geometry_method' in self else self.with_defaults['method']
        _geometry_basis = self['geometry_basis'] if 'geometry_basis' in self else self.with_defaults['basis']
        if self.with_defaults['job_type'][
            :3] == 'OPT':  # and (_geometry_method != self.with_defaults['method'] or _geometry_basis != self.with_defaults['basis']):
            if 'elements' in _geometry_basis:
                return ''
            _geometry_basis = _geometry_basis['default']
            if type(_geometry_method) is list:
                _geometry_method = _geometry_method[-1]
            if _geometry_method[:3].lower() == 'df-':
                _geometry_method = _geometry_method[3:]
            _geometry_method = re.sub('[ru]?ks,', '', _geometry_method, flags=re.IGNORECASE)
            result += '//' + _geometry_method.upper() + '/' + _geometry_basis
        return result

    def parse(self, input: str, debug=False) -> 'InputSpecification':
        r"""
        Take a procedural molpro input, and logically parse it. If it is impossible to do so, an empty specification is returned.

        :param input: Either text that is the input, or a file name containing it.
        """
        if (os.path.isfile(input) and os.access(input, os.R_OK)):
            with open(input, 'r') as f:
                return self.parse(f.read())

        if debug:
            print('InputSpecification.parse() input=', input)

        # print('allowed_methods', self.allowed_methods)
        precursor_methods = ['LOCALI', 'CASSCF', 'OCC', 'CORE', 'CLOSED', 'FROZEN', 'WF',
                             'LOCAL', 'DFIT',
                             'DIRECT', 'EXPLICIT', 'THRESH', 'GTHRESH', 'PRINT', 'GRID']
        df_prefixes = ['', 'DF-']
        postscripts = ['PUT', 'TABLE', 'NOORBITALS', 'NOBASIS']  # FIXME not very satisfactory

        self.clear()
        variables = {}
        geometry_active = False
        self.procname = 'ansatz'
        canonicalised_input_ = re.sub('basis\n(.*)\n *end', r'basis={\1}', input,
                                      flags=re.MULTILINE | re.IGNORECASE | re.DOTALL)
        canonicalised_input_ = re.sub('basis={\n', r'basis={', canonicalised_input_,
                                      flags=re.MULTILINE | re.IGNORECASE | re.DOTALL)
        old_input_ = ''
        count = 100
        while (canonicalised_input_ != old_input_ and count):
            count -= 1
            old_input_ = canonicalised_input_
            canonicalised_input_ = re.sub('basis={([^}]+[^,}])\n([^}]+=[^}]+)}', r'basis={\1,\2}', canonicalised_input_,
                                          flags=re.DOTALL | re.IGNORECASE)
        if not re.match('.*basis={ *s[pdfghi]* *[,}].*', canonicalised_input_, flags=re.DOTALL | re.IGNORECASE):
            canonicalised_input_ = re.sub('basis={ *([^}]*)\n*}', r'basis, \1', canonicalised_input_,
                                          flags=re.DOTALL | re.IGNORECASE)
        canonicalised_input_ = canonicalised_input_.replace('{FREQ}', '{frequencies\nthermo}')  # hack for gmolpro

        if debug:
            print('canonicalised_input_=', canonicalised_input_)

        # parse and protect {....}
        line_end_protected_ = '±'
        for i in range(len(canonicalised_input_)):
            if canonicalised_input_[i] == '{':
                for j in range(i + 1, len(canonicalised_input_)):
                    if canonicalised_input_[j] == '}':
                        canonicalised_input_ = canonicalised_input_[:j] + '}\n' + canonicalised_input_[j + 1:];
                        break
                    elif canonicalised_input_[j] in ';\n':
                        canonicalised_input_ = canonicalised_input_[:j] + line_end_protected_ + canonicalised_input_[
                            j + 1:]
        canonicalised_input_ = canonicalised_input_.replace(';', '\n').replace(line_end_protected_, ';')
        methods_still_possible = True
        methods_started = False
        self['job_type'] = self.with_defaults['job_type']
        for line in canonicalised_input_.split('\n'):
            line = re.sub('basis *,', 'basis=', line, flags=re.IGNORECASE)
            line = re.sub('basis=$,', 'basis=cc-pVDZ-PP', line, flags=re.IGNORECASE)
            group = line.strip()
            if not re.match('.*basis={ *s[pdfghi]* *[,}].*', line, flags=re.DOTALL | re.IGNORECASE):
                line = group.split(line_end_protected_)[0].replace('{', '').strip()
            command = re.sub('[;, !].*$', '', line, flags=re.IGNORECASE).replace('}', '').lower()
            for df_prefix in df_prefixes:
                if command == df_prefix.lower() + 'hf': command = df_prefix.lower() + 'rhf'
                if command == df_prefix.lower() + 'ks': command = df_prefix.lower() + 'rks'
                if command == df_prefix.lower() + 'ldf-ks': command = df_prefix.lower() + 'ldf-rks'
            # print('command', command,'line',line,'group',group)
            for m in _initial_orbital_methods:
                if m.lower() in command.lower() and not any([s + m.lower() in command.lower() for s in ['r', 'u']]):
                    loc = command.lower().index(m.lower())
                    command = re.sub(m.lower(), 'r' + m.lower(), command, flags=re.IGNORECASE)
                    line = re.sub(m.lower(), 'r' + m.lower(), line, flags=re.IGNORECASE)
            if re.match('^orient *, *', line, re.IGNORECASE):
                line = re.sub('^orient *, *', '', line, flags=re.IGNORECASE)
                for orientation_option in _orientation_options.keys():
                    if (line.lower() == _orientation_options[orientation_option].lower()):
                        self['orientation'] = orientation_option
                        break
            elif command.lower() == 'angstrom':
                self['angstrom'] = True
            elif command.lower() == 'proc':
                self.procname = re.sub('^ *proc *,* *', '', line, flags=re.IGNORECASE)
                # _default_job_type_commands[list(_default_job_type_commands.keys())[0]][0]['command'] = self.procname #TODO do this differently
            elif command.lower() == 'endproc':
                pass
            elif ((command.lower() == 'nosym') or (re.match('^symmetry *, *', line, re.IGNORECASE))):
                line = re.sub('^symmetry *, *', '', line, flags=re.IGNORECASE)
                line = "symmetry," + line
                for symmetry_command in _symmetry_commands.keys():
                    if (line.lower() == _symmetry_commands[symmetry_command]):
                        self['symmetry'] = symmetry_command
                        break
            elif re.match('^dkho *=.*', command, re.IGNORECASE):
                self['hamiltonian'] = re.sub('^dkho *= *', 'DK', command, flags=re.IGNORECASE).replace('DK1', 'DK')
            elif line.lower() in properties.values():
                if 'properties' not in self: self['properties'] = []
                self['properties'] += [k for k, v in properties.items() if line.lower() == v]
            elif line.lower().strip().replace('}', '').replace('{', '') in [_local_orbital_types[k]['command'] for k in
                                                                            _local_orbital_types.keys()]:
                for k in _local_orbital_types:
                    if line.lower().strip().replace('}', '').replace('{', '') == _local_orbital_types[k]['command']:
                        if 'orbitals' not in self: self['orbitals'] = []
                        self['orbitals'].append(k)
            elif re.match('^geometry *= *{', group, re.IGNORECASE):
                # print('geometry matched')
                if methods_started: self.data.clear(); return self  # input too complex
                if 'geometry' in self: self.data.clear(); return self  # input too complex
                self['geometry'] = re.sub(';', '\n',
                                          re.sub('^geometry *= *{ *\n*', '', group + '\n', flags=re.IGNORECASE)).strip()
                if '}' in self['geometry']:
                    self['geometry'] = re.sub('}.*$', '', self['geometry']).strip()
                else:
                    geometry_active = True
                # print('self[geometry]',self['geometry'])
            elif geometry_active:
                assert "should not be here" != ""
                self['geometry'] += re.sub(' *[}!].*$', '', line)
                self['geometry'] = self['geometry'].rstrip(' \n') + '\n'
                geometry_active = not re.match('.*}.*', line)
            elif re.match('^geometry *=', line, re.IGNORECASE):
                if methods_started: self.data.clear(); return self  # input too complex
                if 'geometry' in self: self.data.clear(); return self  # input too complex
                self['geometry'] = re.sub('geometry *= *', '', line, flags=re.IGNORECASE)
                self['geometry'] = re.sub(' *!.*', '', self['geometry'])
                # self['geometry_external'] = True
            elif command == 'basis':
                raise ValueError('** warning should not happen basis', line)
            elif re.match('^basis *= *[^{]', line, re.IGNORECASE):
                if methods_started: self.data.clear(); return self  # input too complex
                self['basis'] = {'default': (re.sub(',.*', '', re.sub(' *basis *= *{*(default=)*', '',
                                                                      group.replace('{', '').replace('}', ''),
                                                                      flags=re.IGNORECASE)))}
                fields = line.replace('}', '').split(',')
                self['basis']['elements'] = {}
                for field in fields[1:]:
                    ff = field.split('=')
                    if ff[0].strip(' ')[0] != '!':
                        if len(ff) < 2: self.data.clear(); return self
                        self['basis']['elements'][ff[0][0].upper() + ff[0][1:].lower()] = ff[1].strip('\n ')
                # print('made basis specification',self)
            elif re.match('^basis *=', line, re.IGNORECASE):
                # raise ValueError('unparseable basis', line)
                self.data.clear();
                return self
                pass
            elif re.match('(set,)?[a-z][a-z0-9_]* *=.*$', line, flags=re.IGNORECASE):
                line = re.sub(' *!.*$', '', re.sub('set *,', '', line, flags=re.IGNORECASE)).strip()
                while (
                        newline := re.sub(r'(\[[0-9!]+),', r'\1!',
                                          line)) != line: line = newline  # protect eg occ=[3,1,1]
                fields = _split_comma(line)
                for field in fields:
                    key = re.sub(' *=.*$', '', field)
                    value = re.sub('.*= *', '', field)
                    # print('field, key=', key, 'value=', value)
                    if key == 'charge':
                        self['charge'] = int(value)
                    elif key == 'spin':
                        self['spin'] = int(value)
                    else:
                        variables[key] = value.replace('!', ',')  # unprotect
            elif command in _parameter_commands.values():
                spec_field = [k for k, v in _parameter_commands.items() if v == command][0]
                fields = re.sub('^ *' + command.lower() + ' *,*', '', line.strip().lower(), flags=re.IGNORECASE).split(
                    ',')
                self[spec_field] = {
                    field.split('=')[0].strip().lower(): field.split('=')[1].strip().lower() if len(
                        field.split('=')) > 1 else '' for field in fields}
                if '' in self[spec_field]: del self[spec_field]['']

            elif command == 'core':
                self['core_correlation'] = (line + ',').split(',')[1].lower()
            elif command in [re.sub(r'[,;].*$', '', v[0]) for v in _default_job_type_commands.values() if v]:
                methods_still_possible = False
                if command[:4] == 'freq' and self['job_type'] == 'OPT':
                    self['job_type'] = 'OPT+FREQ'
                    self['job_type_commands'][self['job_type']] = self['job_type_commands'].pop('OPT')
                else:
                    self['job_type'] = \
                        [k for k, v in _default_job_type_commands.items() if
                         v and re.sub(r'[,;].*$', '', v[0]) == command][
                            0]
                    self['job_type_commands'] = {}
                    self['job_type_commands'][self['job_type']] = []
                self['job_type_commands'][self['job_type']].append(
                    line.replace('}', '').replace('{', '').replace(',proc=' + self.procname, ''))
                # print('job_type', self['job_type'])
                # print('job_type_commands', self['job_type_commands'])
            elif methods_still_possible and any([re.fullmatch('{?' + df_prefix + re.escape(method), command,
                                                              flags=re.IGNORECASE) for
                                                 df_prefix
                                                 in df_prefixes
                                                 for method in self.allowed_methods]):
                step = {}
                method_with_options = re.sub('^{', '', re.sub('}$', '', group))
                # print('matched for method; line=',line,'command=',command,'group=',group)
                method_ = command
                if command[:3] == 'df-':
                    self['density_fitting'] = True
                    method_ = command[3:]
                elif command[:4] == 'pno-' or command[:4] == 'ldf-':
                    self['density_fitting'] = True
                elif False and 'density_fitting' in self and self['density_fitting'] and not any(
                        [step_['command'] == command for job_type in _default_job_type_commands for step_ in
                         _default_job_type_commands[job_type]]):  # what was this all about?
                    self.data.clear()
                    return self
                # method_options = (re.sub(';.*$', '', line.lower()).replace('}', '') + ',').split(',', 1)[1]

                # method_options_ = _split_comma(method_options.strip(', \n'))
                # if method_options_ and method_options_[-1] == '': method_options_ = method_options_[:-2]
                # print('method_options_',method_options_)
                step = JobStep(line.replace(command, method_).replace('}', '').lower())
                # step['command'] = method_
                # if method_options_:
                #     step['options'] = method_options_
                # TODO parsing of extras from following directives
                # print('group before directives',group)
                # directives = group.replace('}', '').split(';')[1:]
                # print('directives', directives)
                # print('intial step', step)
                # for directive in directives:
                #     cmd, opt = (directive + ',').split(',', 1)
                #     print('cmd',cmd,'opt',opt)
                # opts = {m1.split('=')[0].strip(): (m1.split('=')[1].strip() if len(m1.split('=')) > 1 else '') for
                #         m1 in opt.rstrip(',').split(',')}
                # if '' in opts: del opts['']
                # if 'directives' not in step: step['directives'] = []
                # opts = opt.rstrip(',').split(',')
                # if opts and opts[-1] == '': opts = opts[:-2]
                # d = {'command': cmd}
                # if opts: d['options'] = opts
                # step['directives'].append(d)
                if 'method' not in self:
                    self['method'] = []
                # method_ = command
                # print('step.dump()',step.dump())
                self['method'].append(step.dump().replace('{', '').replace('}', ''))
                methods_started = True
                # print('self[method]', self['method'])
            elif not methods_still_possible and line:
                if 'epilogue' not in self:
                    self['epilogue'] = []
                self['epilogue'].append(line.replace('}', '').lower())

        if 'job_type_commands' in self and self['job_type_commands'].get('job_type', '') == \
                schema['properties']['job_type_commands']['items'][self['job_type']]['default']:
            self.pop('job_type_commands')
        # if 'method' not in self and 'precursor_methods' in self:
        #     parse_method(self, self['precursor_methods'][-1])
        #     self['precursor_methods'].pop()
        if variables:
            self['variables'] = variables
        if 'hamiltonian' not in self:
            self['hamiltonian'] = self.basis_hamiltonian
        # self.regularise_procedure_references()

        # spin_ = self.open_shell_electrons
        # print('initial spin_',spin_)
        # if 'variables' in self and 'spin' in self['variables'] and int(self['variables']['spin'])%2 == spin_%2:
        #     spin_ = self['variables']['spin']
        # if 'variables' not in self: self['variables'] = {}
        # self['variables']['spin'] = spin_

        # print('before deduce_job_type', self)
        # self.deduce_job_type()
        # print('after deduce_job_type', self)
        # print('final self',self)
        self.validate()
        return self

    # def regularise_procedure_references(self):
    #     for step in self['steps']:
    #         job_type_commands = [job_type_ste['command'] for job_type_step in _default_job_type_commands.values() for
    #                              job_type_ste
    #                              in job_type_step if
    #                              'command' in job_type_ste and job_type_ste['command'] != self.procname] + [
    #                                 'frequencies']
    #         if 'command' in step.keys() and step['command'] in job_type_commands:
    #             step['options'] = [option for option in step['options'] if
    #                                'proc' not in option] if 'options' in step else []
    #             step['options'].append('proc=' + self.procname)

    def molpro_input(self) -> str:
        r"""
        Create the procedural Molpro input from the declarative specification
        """
        _input = ''
        defaulted_spec = self.with_defaults
        if 'prologue' in self:
            if type(self['prologue']) is list:
                for prologue in self['prologue']:
                    _input += prologue + '\n'
            else:
                _input += self['prologue'] + '\n'

        if 'orientation' in self:
            _input += 'orient,' + _orientation_options[self['orientation']] + '\n'

        if 'symmetry' in self:
            _input += _symmetry_commands[self['symmetry']] + '\n'

        if 'angstrom' in self and self['angstrom']:
            _input += 'angstrom' + '\n'

        if 'geometry' in self:
            _geometry = self['geometry'].strip()
            if (os.path.isfile(_geometry) and os.access(_geometry, os.R_OK)) or re.match(r'.*\.(xyz|h5)$',
                                                                                         _geometry) or (
                    _geometry[0] == '{' and _geometry[-1] == '}'):
                _input += 'geometry=' + _geometry + '\n'
            else:
                _input += 'geometry=' + '{' + _geometry + '}' + '\n'

        _job_type = defaulted_spec['job_type']
        if _job_type[:3] == 'OPT':
            if 'geometry_basis' in self:
                _input += self._input_from_basis(self['geometry_basis'])
            elif 'basis' in self:
                _input += self._input_from_basis(self['basis'])
        else:
            if 'basis' in self:
                _input += self._input_from_basis(self['basis'])

        if 'hamiltonian' in self and self['hamiltonian'][:2] == 'DK':
            _input += 'dkho=' + self['hamiltonian'][2] if len(self['hamiltonian']) > 2 else '1' + '\n'

        if 'charge' in self:
            _input += 'charge=' + str(self['charge']) + '\n'

        if 'spin' in self:
            # print('spin=' + self['spin'])
            _input += 'spin=' + str(self['spin']) + '\n'

        if 'variables' in self:
            for k, v in self['variables'].items():
                if v != '':
                    _input += k + '=' + v + '\n'

        if 'properties' in self:
            for p in self['properties']:
                _input += properties[p] + '\n'  # TODO review. It looks wrong

        for typ, command in _parameter_commands.items():
            if typ in self and len(self[typ]) > 0:
                _input += command
                for k, v in self[typ].items():
                    _input += ',' + k.lower() + ('=' + str(v) if str(v) != '' else '')
                _input += '\n'

        if 'core_correlation' in self:
            _input += 'core,' + self['core_correlation'] + '\n'

        _job_type_commands = defaulted_spec['job_type_commands'][_job_type]
        if len(_job_type_commands) > 0:
            _input += '\nproc ' + self.procname + '\n'
        _method = self.with_defaults['geometry_method'] if 'geometry_method' in self else self.with_defaults['method']
        _input += self._input_from_method(_method)

        if len(_job_type_commands) > 0:
            _input += 'endproc\n\n'
            for step in _job_type_commands:
                _step = JobStep(step)
                for option in _step.options:
                    if re.match(r'^proc=.*$', option):
                        _index = step.options.index(option)
                try:
                    _step.options[_index] = 'proc=' + self.procname
                except NameError:
                    _step.options.append('proc=' + self.procname)
                _input += _step.dump() + '\n'

        if _job_type[:3] == 'OPT' and (self.with_defaults['method'] != self.with_defaults['geometry_method']
                                       or self.with_defaults['geometry_basis'] != self.with_defaults['basis']
        ):
            if 'geometry_basis' in self and 'basis' in self and self.with_defaults['geometry_basis'] != \
                    self.with_defaults['basis']:
                _input += self._input_from_basis(self['basis'])
            _input += self._input_from_method(self.with_defaults['method'])

        if 'extrapolate' in self:
            _input += 'extrapolate,' + self['extrapolate'] + '\n'

        if 'orbitals' in self:
            for k in self['orbitals']:
                if _local_orbital_types[k]['command'].strip(): _input += '{' + _local_orbital_types[k][
                    'command'] + '}\n'

        if 'epilogue' in self:
            if type(self['epilogue']) is list:
                for epilogue in self['epilogue']:
                    _input += epilogue + '\n'
            else:
                _input += self['epilogue'] + '\n'
        return _input.rstrip('\n') + '\n'

    def _input_from_basis(self, basis) -> str:
        _input = 'basis=' + basis['default']
        if 'elements' in basis:
            for e, b in basis['elements'].items():
                _input += ',' + e + '=' + b
        _input += '\n'
        return _input

    def _input_from_method(self, method) -> str:
        _input = ''
        if type(method) is str:
            method = [method]
        for step in method:
            _input += '{'
            if 'density_fitting' in self and self['density_fitting'] and step.lower()[:4] != 'pno-' and \
                    step.lower()[:4] != 'ldf-':
                _input += 'df-'
            _input += step
            if False:  # TODO don't lose this logic
                if 'options' in step:
                    numerical_options = []
                    for option in step['options']:
                        if re.match('[0123456789]+', option):
                            numerical_options.append(option)
                    other_options = step['options'][len(numerical_options):]
                    other_options.sort()
                    step['options'] = numerical_options + other_options
                    for option in step['options']:
                        _input += ',' + str(option)
                if 'directives' in step:
                    for directive in step['directives']:
                        _input += ';' + directive['command']
                        if 'options' in directive:
                            for option in directive['options']:
                                _input += ',' + str(option)
            _input += '}\n'
        return _input

    def set_job_type(self, new_job_type):
        if self['job_type'] == new_job_type: return
        self['job_type'] = new_job_type
        if 'steps' not in self: self['steps'] = []
        old_len = len(self['steps'])
        new_steps = [
                        step for step in self['steps']
                        if
                        step['command'].lower() not in [s['command'].lower() for j in _default_job_type_commands for s
                                                        in
                                                        _default_job_type_commands[j]]
                    ] + _default_job_type_commands[new_job_type]
        for new_step in new_steps:
            for step in self['steps']:
                if step['command'].lower() == new_step['command'].lower() and 'options' in step:
                    new_step['options'] = step['options']
        self['steps'] = new_steps

    @property
    def method(self) -> str:
        r"""
        The single method implemented by the job, represented as the Molpro input command name that implements it.
        """
        methods = []
        defaulted_spec = self.with_defaults
        if 'method' in defaulted_spec:
            methods = defaulted_spec['method']
        if type(methods) is str:
            return methods.replace('\n', ';').split(';')[0].split(',')[0]
        else:
            return methods[-1].replace('\n', ';').split(';')[0].split(',')[0]

    @method.setter
    def method(self, method: str):
        if method is None or method == '' or (self.method is not None and method.lower() == self.method.lower()): return
        if method.lower() not in [m.lower() for m in self._hartree_fock_methods]:
            self['method'] = ['rhf' if method[0].lower() != 'u' else 'uhf', method.lower()]
        else:
            self['method'] = method.lower()

    @property
    def method_options(self) -> str:
        r"""The options for a single-method job
        """
        _command = self.method.split(',')[0].split(';')[0]
        if 'method' not in self:
            return []
        for step in self['method'] if type(self['method']) is list else [self['method']]:
            js = JobStep(step)
            if js.command == _command:
                return js.options
        return []

    @method_options.setter
    def method_options(self, options):
        if 'method' not in self:
            self['method'] = self.with_defaults['method']
        _command = self.method.split(',')[0].split(';')[0]
        if not isinstance(self['method'], list):
            self['method'] = [self['method']]
        for index, step in enumerate(self['method']):
            js = JobStep(step)
            if js.command == _command:
                js.options = [k + '=' + v for k, v in options.items()] if type(options) is not list else options
                self['method'][index] = js.dump().replace('{', '').replace('}', '')
        if len(self['method']) == 1:
            self['method'] = self['method'][0]

    @property
    def basis_quality(self) -> int:
        r"""The cardinal number of the basis set used in the job. If not a correlation-consistent basis set, returns 0."""
        quality_letters = {2: 'D', 3: 'T', 4: 'Q', 5: '5', 6: '6', 7: '7'}
        if 'basis' in self:
            bases = [self['basis']['default']]
            if 'elements' in self['basis']: bases += self['basis']['elements'].values()
            qualities = []
            for basis in bases:
                quality = 0
                for q, l in quality_letters.items():
                    if re.match(r'.*V\(?.*' + l, basis, flags=re.IGNORECASE): quality = q
                qualities.append(quality)
            if all(quality == qualities[0] for quality in qualities):
                return qualities[0]
        return 0

    @property
    def basis_hamiltonian(self) -> str:
        """The hamiltonian for which the orbital basis set is designed"""
        result = 'AE'
        for v, k in _hamiltonians.items():
            if k and 'basis' in self and 'default' in self['basis'] and k['basis_string'] in \
                    self['basis']['default']: result = v
        if 'variables' in self and 'dkho' in self['variables']:
            result = 'DK' + str(self['variables']['dkho']) if str(
                self['variables']['dkho']) != '1' else 'DK'
        return result

    _last_density_functional = 'LDA'

    @property
    def density_functional(self) -> str:
        """For a Kohn-Sham calculation, the name of the density functional used"""
        if 'method' not in self or not self['method']:
            return None
        methods = self['method'] if type(self['method']) is list else [self['method']]
        js = JobStep(methods[0])
        if js.command not in ['ks', 'rks', 'uks']:
            return None
        if not js.options:
            return self._last_density_functional
        self._last_density_functional = js.options[0].upper()
        return js.options[0].upper()

    @density_functional.setter
    def density_functional(self, density_functional):
        if 'method' not in self or not self['method']:
            return
        methods = self['method'] if type(self['method']) is list else [self['method']]
        js = JobStep(methods[0])
        if js.command not in ['ks', 'rks', 'uks']:
            return
        if js.options and '=' not in js.options[0]:
            js.options[0] = density_functional
        else:
            js.options = [density_functional] + js.options
        methods[0] = js.dump().replace('{', '').replace('}', '')
        self['method'] = methods if len(methods) > 1 else methods[0]

    @property
    def open_shell_electrons(self) -> int:
        r"""
        The number of open-shell electrons in the molecule's normal state.  This will typically be 0 or 1, but for some special cases (eg atoms) might be higher.
        """
        if 'spin' in self:
            return self['spin']
        # TODO set up a cache if input has not changed and geometry file has not changed
        if 'geometry' not in self: return 0
        try:
            with open(pathlib.Path(self.directory if self.directory is not None else '.') / self['geometry'],
                      'r') as f:
                geometry = ''.join(f.readlines())
        except:
            geometry = self['geometry']
        # print('geometry',geometry)
        line_number = 0
        start_line = 1
        total_nuclear_charge = 0
        for line in geometry.replace(';', '\n').split('\n'):
            line_number += 1
            if line.strip().isdigit() and line_number == 1: start_line = 3
            if line_number >= start_line and line.strip():
                word = re.sub(r'\s+', ',', line.strip()).split(',')[0]
                word = re.sub(r'\d.*$', '', word[0].upper() + word[1:].lower())
                try:
                    atomic_number = periodic_table.index(word) + 1
                except ValueError:
                    atomic_number = 0
                total_nuclear_charge += atomic_number
        charge = self['charge'] if 'charge' in self else 0
        total_electrons = total_nuclear_charge - charge
        electrons = total_electrons % 2
        # implementing default spin > 1 is tricky because of handling of input files that do not contain spin specification
        # if atomic_number == total_nuclear_charge:
        #     if total_electrons in [6, 8, 14, 16, 32, 34, 50, 52, 82, 84]: electrons = 2
        #     if total_electrons in [7, 15, 33, 51, 83]: electrons = 3
        # print('Electrons: ' + str(electrons))
        return electrons

    def polish(self):
        """Ensure that the job specification is valid and complete"""
        self._clean_coupled_cluster_property_input()

    def _clean_coupled_cluster_property_input(self):
        for step_index, step in enumerate(self.job_steps):
            # print('clean_coupled_cluster_property_input, step',step)
            if step.command.lower()[:4] in ['ccsd', 'bccd', 'qcisd']:
                for index, directive in enumerate(step.directives):
                    _directive = JobStep(directive)
                    if _directive.command.lower() == 'expec':
                        operator = _directive.options[0].lower().replace('expec,', '')
                        property = [k for k, v in properties.items() if v == 'gexpec,' + operator][0]
                        if 'properties' in self and property not in \
                                self['properties']:
                            # step.directives.remove(_directive)
                            step.directives.pop(index)
                        self.set_job_step(step, step_index)
                if 'properties' in self:
                    for property in self['properties']:
                        cmd = properties[property]
                        operator = cmd.lower().replace('gexpec,', '').strip()
                        directive = 'expec,' + operator
                        if directive not in step.directives:
                            step.directives.append(directive)
                            self.set_job_step(step, step_index)


def canonicalise(input):
    result = re.sub('\n}', '}',
                    re.sub(' *= *', '=',
                           re.sub('{\n', r'{',
                                  re.sub('\n+', '\n',
                                         re.sub(' *, *', ',',
                                                input.replace(';',
                                                              '\n')))))).rstrip(
        '\n ').lstrip(
        '\n ') + '\n'
    result = re.sub(',+}', '}', result)
    # push variable assignments below geometry=file.xyz to hack compatibility with gmolpro guided
    # print('before hack', result)
    # hack for gmolpro geomtyp:
    old_result = ''
    while (old_result != result):
        old_result = result
        result = re.sub('(\\w+=\\w+)\n(orient,mass)', '\\2\n\\1', result, flags=re.MULTILINE | re.IGNORECASE)
    old_result = ''
    while (old_result != result):
        old_result = result
        result = re.sub('(\\w+=\\w+)\n(nosym)', '\\2\n\\1', result, flags=re.MULTILINE | re.IGNORECASE)
    old_result = ''
    while (old_result != result):
        old_result = result
        result = re.sub('(\\w+=\\w+)\n(geometry=[\\w.{}]*)', '\\2\n\\1', result, flags=re.MULTILINE | re.IGNORECASE)
    old_result = ''
    while (old_result != result):
        old_result = result
        result = re.sub('(\\w+=\\w+)\n(basis={[^\n]*})', '\\2\n\\1', result,
                        flags=re.MULTILINE | re.IGNORECASE | re.DOTALL)
    # print('after 1st hack', result)
    result = re.sub('(dkho=\\d)\n(geomtyp=xyz)', '\\2\n\\1', result, flags=re.MULTILINE | re.IGNORECASE)
    # hack for gmolpro-style frequencies:
    # print('after 2nd hack', result)
    result = result.replace('{FREQ}', '{frequencies\nthermo}')
    result = re.sub('basis={\n', 'basis={', result, flags=re.IGNORECASE | re.DOTALL)
    # print('after 3rd hack', result)
    new_result = ''
    in_group = False
    for line in re.sub('set[, ]', '', result.strip(), flags=re.IGNORECASE).split('\n'):

        if not in_group:
            in_group = '{' in line
        # transform out alternate formats of basis
        line = re.sub('basis *, *', 'basis=', line.rstrip(' ,'), flags=re.IGNORECASE)
        line = re.sub('basis= *{(.*)} *(!.*)?$', r'basis=\1 \2', line, flags=re.IGNORECASE)
        line = re.sub('basis= *default *= *', r'basis=', line, flags=re.IGNORECASE).lower()
        line = re.sub(' *!.*$', '', line)
        for cmd in ['hf', 'ks']:
            for bra in ['', '{']:
                line = re.sub('^ *' + bra + ' *' + cmd, bra + 'r' + cmd, line,
                              flags=re.IGNORECASE)  # TODO unify with following

        # transform out alternate spin markers
        # for m in initial_orbital_methods:
        #     line = re.sub('r' + m, m, line, flags=re.IGNORECASE)
        # transform in alternate spin markers
        for m in _initial_orbital_methods:
            line = re.sub('^{' + m, '{r' + m.lower(), line, flags=re.IGNORECASE)

        if line.lower().strip() in _job_type_aliases.keys(): line = _job_type_aliases[line.lower().strip()]
        if line.lower().strip() in _symmetry_command_aliases.keys():
            line = _symmetry_command_aliases[line.lower().strip()]
        line = line.replace('!', '&&&&&')  # protect trailing comments
        while (newline := re.sub(r'(\[[0-9!]+),', r'\1!', line)) != line: line = newline  # protect eg occ=[3,1,1]
        if re.match(r'[a-z][a-z0-9_]* *= *\[?[!a-z0-9_. ]*\]? *,', line, flags=re.IGNORECASE):
            line = line.replace(',', '\n')
        line = re.sub(' *}', '}', line)
        line = re.sub('{ *', '{', line)
        line = line.replace('!', ',').strip() + '\n'  # unprotect
        line = line.replace('&&&&&', '!').strip() + '\n'  # unprotect
        # print('line before bracketing',line, in_group)
        if line.strip() and line.strip()[0] != '{' and not re.match(r'^ *\w+ *=', line) and not in_group and not any(
                [v in line for v in _parameter_commands.values()]):
            comment_split = line.split('!')
            line = '{' + comment_split[0].strip() + '}'  # + (comment_split[1] if len(comment_split) > 1 else '')
        # print('line after bracketing',line)
        in_group = in_group and not '}' in line
        if line.strip('\n') != '':
            new_result += line.strip('\n ') + '\n'
    return new_result.strip('\n ') + '\n'


def equivalent(input1, input2, debug=False):
    if isinstance(input1, InputSpecification): return equivalent(input1.molpro_input(), input2, debug)
    if isinstance(input2, InputSpecification): return equivalent(input1, input2.molpro_input(), debug)
    if debug:
        _logger.debug('equivalent: input1=', input1)
        _logger.debug('equivalent: input2=', input2)
        _logger.debug('equivalent: canonicalise(input1)=', canonicalise(input1))
        _logger.debug('equivalent: canonicalise(input2)=', canonicalise(input2))
        _logger.debug('will return this', canonicalise(input1).lower() == canonicalise(input2).lower())
    return canonicalise(input1).lower() == canonicalise(input2).lower()


def _split_comma(string):
    string_ = string
    string__ = re.sub(r'(\[.*),(.*\])', r'\1@@@\2', string_)
    while string__ != string_:
        string_ = string__
        string__ = re.sub(r'(\[.*),(.*\])', r'\1@@@\2', string_)
    return [field.replace('@@@', ',') for field in (string__.split(','))]

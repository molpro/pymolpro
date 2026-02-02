import atexit
import copy

import math
import os
try:
    import pwd
except ImportError:
    pass
import pathlib
import tempfile

import pysjef
import subprocess
import re
import json
import numpy as np
import shutil

from pandas.core.window.doc import kwargs_scipy
from scipy.special import factorial2
import inspect

import pymolpro
from .molpro_input import schema, InputSpecification

from .elements import periodic_table
import logging
logger = logging.getLogger(__name__)

def no_errors(projects, ignore_warning=True):
    """
    Checks that none of the projects have any errors. Projects can by running.

    :param projects: list of projects to check for errors
    :param ignore_warning: Whether to count warnings as well as errors.
    :return: True/False whether any projects have errors
    """
    return not any([p.xpath('//error[@type != "Warning"]' if ignore_warning else '//error') for p in projects])


def element_to_dict(node, attributes=True):
    """
    Convert an lxml.etree node tree into a dict.

    :param node: A node in an lxml tree
    :type node: lxml.etree.ElementTree
    :param attributes: whether to include attributes in the result
    :type attributes: bool, optional
    :return: A dictionary representing the tree
    :rtype: dict
    """
    result = {}
    if attributes:
        for item in node.attrib.items():
            key, result[key] = item

    for element in node.iterchildren():
        # Remove namespace prefix
        key = element.tag.split('}')[1] if '}' in element.tag else element.tag

        # Process element as tree element if the inner XML contains non-whitespace content
        if element.text and element.text.strip():
            value = element.text
        else:
            value = element_to_dict(element)
        if key in result:
            if type(result[key]) is list:
                result[key].append(value)
            else:
                result[key] = [result[key], value]
        else:
            result[key] = value
    return result


def resolve_geometry(geometry):
    import re
    if re.match(re.compile(
            r'^(?:http|ftp|file)s?://'  # http:// or https://
            r'(?:(?:[A-Z0-9](?:[A-Z0-9-]{0,61}[A-Z0-9])?\.)+(?:[A-Z]{2,6}\.?|[A-Z0-9-]{2,}\.?)|'  # domain...
            r'localhost|'  # localhost...
            r'\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})'  # ...or ip
            r'(?::\d+)?'  # optional port
            r'(?:/?|[/?]\S+)$', re.IGNORECASE), geometry) is not None:
        import urllib.request
        return urllib.request.urlopen(geometry).read()
    else:
        try:
            with open(geometry, "r") as file:
                return file.read()
        except Exception:
            return geometry


def method_from_commands(commands):
    method = re.sub(r'^(df-)?[ur]?hf[\s;]*$', 'HF', commands.strip(), flags=re.IGNORECASE)
    method = re.sub(r'^(df-)?[ur]?hf\s*; *', '', method.strip(), flags=re.IGNORECASE)
    method = re.sub(r'^(df-)?[ur]?ks\s*,', '', method, flags=re.IGNORECASE)
    return method


def molpro_xyz_powers(l):
    components = ['',
                  'x', 'y', 'z',
                  'xx', 'yy', 'zz', 'xy', 'xz', 'yz',
                  'xxx', 'yyy', 'zzz', 'xxy', 'xxz', 'xyy', 'yyz', 'xzz', 'yzz', 'xyz',
                  'xxxx', 'yyyy', 'zzzz', 'xxxy', 'xxxz', 'xyyy', 'yyyz', 'xzzz', 'yzzz', 'xxyy', 'xxzz', 'yyzz',
                  'xxyz', 'xyyz', 'xyzz']
    if l < 5:
        result = [[components[l * (l + 1) * (l + 2) // 6 + m].count('xyz'[xyz]) for xyz in range(3)] for m in
                  range((l + 1) * (l + 2) // 2)]
    else:
        result = [[p[xyz] for xyz in range(3)] for p in lexical_xyz_powers(l)]
    return result


def lexical_xyz_powers(l):
    powers = [l, 0, 0]
    while True:
        powers[2] = l - powers[1] - powers[0]
        yield powers
        if powers[1] > 0:
            powers[1] -= 1
        elif powers[0] == 0:
            return
        else:
            powers[0] -= 1
            powers[1] = l - powers[0]


class Project(pysjef.project.Project):
    r"""
    A :py:class:`Project` holds all the data associated with a single Molpro job. This includes input, output, any auxiliary files, and information about
    the state of the job. :py:class:`Project` acts as a reference to a
    `sjef <https://molpro.github.io/sjef>`_
    Project with some added functionality.
    All of the data is stored in the project bundle on the file system, with the consequence that it is safe to construct and recreate multiple
    :py:class:`Project` objects all mapping the same underlying project.

    The underlying functionality of `sjef` includes job submission and monitoring on local or remote machines, and structured searching of Molpro's xml output stream.
    This class provides additional convenience operations:

    * Preparation of simple Molpro input from provided geometry, method, basis set and other options.

    * Molpro-specific error checking.

    * Computed properties, energies and geometries.

    * Computed orbitals, including their evaluation on a grid.

    * Values of Molpro variables.


    """

    def __init__(self, name: str = None, input: str | dict | None = None, specification: str | dict | None = None, ansatz: str | None = None,
                 files: list[str] | None = None, **kwargs):
        r"""

        :param name: The base filename of the filesystem bundle carrying the project. If the bundle does not yet exist, it is created.
        :param input: General specification of the input. If it looks like JSON or is a dictionary, it is treated as if it had been passed as the specification parameter. If it has the form method/basis or method/basis//geometry_method/geometry_basis, it is treated as if it had been passed as the ansatz parameter. Otherwise, it is treated as the desired contents of the Molpro input file. If input is specified, all other arguments are ignored.
        :param specification: Either a dictionary or a JSON string conforming to the JSON schema https://www.molpro.net/schema/molpro_input.json
        :param ansatz: A string of the form method/basis//geometry_method/geometry_basis or method/basis which is parsed to give the same effect as the method and basis parameters. If geometry_method/geometry_basis is specified, the calculation will be preceded by a geometry optimisation at that level of theory.
        :param files: External files to be copied into the project directory. If one of these is a molpro-output xml file, it is used to construct the input file.
        :param kwargs: Any of the top-level keywords in the JSON schema https://www.molpro.net/schema/molpro_input.json, or any of the arguments accepted by the parent sjef.Project class constructor.
"""
        # print('Project(): name',name,'input',input,'specification',specification,'ansatz',ansatz,'files',files,'kwargs',kwargs)
        possible_arguments_matching_schema = [k for k in schema['properties'] if
                                              k not in inspect.signature(self.__init__).parameters]
        # print('possible_arguments_matching_schema',possible_arguments_matching_schema)
        # print('remaining arguments',{k: v for k, v in kwargs.items() if
        #                             k not in possible_arguments_matching_schema + ['suffix']})
        if not hasattr(self, '__initialized'):
            try:
                kwargs_ = {k: v for k, v in kwargs.items() if k not in possible_arguments_matching_schema + ['suffix']}
                if not name:
                    from pysjef import __version__ as pysjef_version
                    from packaging.version import Version
                    _name = self._anonymous_name(input, specification, ansatz, **kwargs)
                    if Version(pysjef_version) >= Version("1.42.1"):
                        kwargs_['record_as_recent'] = False
                    kwargs_['location'] = (pathlib.Path(tempfile.gettempdir()) / 'pymolpro_projects').as_posix()
                    os.makedirs(kwargs_['location'], exist_ok=True)
                    atexit.register(self._unconditionally_destroy)
                else:
                    _name = name

                super().__init__(name=_name, suffix='molpro', **kwargs_)
                self.__initialized = True
            except Exception:
                raise FileNotFoundError("Cannot open project " + name)

        self.local_molpro_root_ = None

        self.input(input, specification, ansatz, **kwargs)

        if files is not None and isinstance(files, list) and len(files) > 0:
            project_directory = pathlib.Path(self.filename('', '', -1))
            run_directory = project_directory
            project_name = str(project_directory.stem)
            rundir = False
            for file in files:
                if (os.path.isfile(file)):
                    path = pathlib.Path(file)
                    if path.suffix in ['.out', '.xml']:
                        if not rundir:
                            self.run_directory_new()
                            rundir = True
                        shutil.copyfile(file, self.filename(path.suffix[1:]))
                        if path.suffix == '.xml':
                            input_from_output = self.input_from_output()
                            self.write_input('\n'.join(input_from_output))
                    elif path.suffix in ['.inp', '.xyz'] and path.stem != 'optimised':
                        shutil.copyfile(file, pathlib.Path(self.filename('', run=-1)) / path.name)
                    else:
                        if not rundir:
                            self.run_directory_new()
                            rundir = True
                        shutil.copyfile(file, pathlib.Path(self.filename('')) / path.name)

    def _unconditionally_destroy(self):
        try:
            self.erase()
            shutil.rmtree(self.filename(), ignore_errors=True)
        except:
            pass


    def input(self, input: str | dict = None, specification: str | dict = None, ansatz: str = None, **kwargs):
        r"""
        Defines the input for running Molpro.

        :param input:General specification of the input. If it looks like JSON or is a dictionary, it is treated as if it had been passed as the specification parameter. If it has the form method/basis or method/basis//geometry_method/geometry_basis, it is treated as if it had been passed as the ansatz parameter. Otherwise, it is treated as the desired contents of the Molpro input file. If input is specified, all other arguments are ignored.
        :param specification: Either a dictionary or a JSON string conforming to the JSON schema https://www.molpro.net/schema/molpro_input.json
        :param ansatz: A string of the form method/basis//geometry_method/geometry_basis or method/basis which is parsed to give the same effect as the method and basis parameters. If geometry_method/geometry_basis is specified, the calculation will be preceded by a geometry optimisation at that level of theory.
        :param kwargs: Any of the top-level keywords in the JSON schema https://www.molpro.net/schema/molpro_input.json
        """

        if input is not None and (isinstance(input, dict) or re.match(r'^{.*}$', input)):
            self.input(specification=input, **kwargs)
            return
        elif input is not None and re.match(r'^[\w][\w\d/]*[\w\d]$', input):
            self.input(ansatz=input, **kwargs)
            return
        elif input is not None:
            self.write_input(input)
            return

        if specification is not None and type(specification) is str:
            try:
                self.input(specification=json.loads(specification), **kwargs)
                return
            except:
                pass

        if specification is not None and isinstance(specification, dict):
            self.input_specification = InputSpecification(specification=specification)
            self.write_input(self.input_specification.molpro_input())
            self.property_set({'input_specification': json.dumps(dict(self.input_specification))})

        if ansatz is not None:
            _parsed_ansatz = self.parse_ansatz(ansatz)
            _kwargs = {k: v for k, v in kwargs.items()}
            for k, v in _parsed_ansatz.items():
                if k != 'ansatz' and v is not None:  _kwargs[k] = v
            if 'geometry_method' in _parsed_ansatz and _parsed_ansatz[
                'geometry_method'] is not None and 'job_type' not in _kwargs:
                _kwargs['job_type'] = 'OPT'
            # print('ansatz processing', _kwargs)
            self.input(**_kwargs)
            return

        if specification is not None and isinstance(specification, dict):
            _input = copy.deepcopy(specification)
        else:
            _input = {}
        for key in [k for k in kwargs if k in schema['properties']]:
            if kwargs[key] is not None:
                _input[key] = kwargs[key]
        for key in ['basis', 'geometry_basis']:
            if key in _input and type(_input[key]) is str:
                _input[key] = {"default": _input[key]}
        if _input:
            self.input_specification = InputSpecification(specification=_input)
            self.write_input(self.input_specification.molpro_input())
            self.property_set({'input_specification': json.dumps(dict(self.input_specification))})

    # def set_method(self,method,basis="cc-pVTZ",geometry_method=None, geometry_basis=None):
    #     pass
    def commandify_method(self, method):
        if method.strip().upper() in self.registry('DFUNC').keys():
            return 'df-ks,' + method
        __method = method.lower().strip()
        if re.match(r'^(df-)?[ur]?(hf|ks).*', __method):
            return method
        if __method[-2:] != 'hf' and __method[1:3] != 'hf' and __method[:2] != 'hf' and __method[
                                                                                        -2:] != 'ks' and 'ks,' not in __method and 'ks ' not in __method:
            if __method[:2] == 'df':
                if __method[:4] == 'df-u':
                    __method = 'df-uhf; df-' + __method[4:]
                else:
                    __method = 'df-hf; ' + __method
            else:
                if __method[:1] == 'u':
                    __method = 'uhf; ' + __method[1:]
                else:
                    __method = 'hf; ' + __method
        return __method


    @property
    def ansatz(self):
        return self.input_specification.ansatz

    def parse_ansatz(self, ansatz):
        parsed = {}
        parts = ansatz.split('//')
        # if len(parts) == 1:
        #     parts.append(parts[0])
        result = []
        for part in parts:
            bits = part.split('/')
            # print('part',part,'bits',bits)
            if bits[0].lower()[:3] == 'df-':
                parsed['density_fitting'] = True
                bits[0] = bits[0][3:]
            # if bits[0].strip().lower() == 'hf':
            #     bits[0] = 'df-hf'
            # elif bits[0].strip().lower() == 'mp2':
            #     bits[0] = 'df-hf;df-mp2'
            elif bits[0].strip().upper() in self.registry('DFUNC').keys():
                bits[0] = 'ks,' + bits[0]
            # elif bits[0].strip().lower()[:3] == 'df-' and not re.match(bits[0].strip().lower()[3:], '[ur]?hf'):
            #     bits[0] = 'hf; ' + bits[0]
            # else:
            #     bits[0] = 'hf; ' + bits[0]
            result += bits
        parsed['ansatz'] = ansatz
        parsed['method'] = result[0]
        parsed['basis'] = result[1] if len(result) > 1 else None
        parsed['geometry_method'] = result[2] if len(result) > 2 else None
        parsed['geometry_basis'] = result[3] if len(result) > 3 else None
        return parsed

    def errors(self, ignore_warning=True):
        '''
        Return all error nodes

        :return: list of error nodes
        :rtype: lxml.etree
        '''
        errors = self.select('//error')
        if ignore_warning:
            errors = [x for x in errors if
                      not ('type' in x.attributes and x.attributes['type'].lower() == 'warning')]
        return errors

    def properties_old(self, name='Energy', principal=True, *, value=False, **kwargs, ):
        """
        Obtain selected properties from the job output

        :param name: name of property
        :param principal: principal property
        :param value: return by value
        :param kwargs: any other attribute selectors for the select function
        :return: list of properties

        Shorthand for `p.select('//properties[name= {}, principal ={}, ... ]')`
        """
        string = f'name={name}'
        if principal:
            string += ', principal=true'
        for key, val in kwargs.items():
            string += f',{key}={val}'
        selector = f'//property[{string}]'
        if value:
            selector += '.value'
        return self.select(selector)

    def properties(self, *args, **kwargs):
        """
        Obtain selected properties from the job output

        :param args: any trailing XPath qualifier, for example ``'[@StateSymmetry = "1" or @StateSymmetry = "4"]'``. For restrictions that should be combined with ``and``, it's simpler to instead use a kwarg argument, eg ``StateSymmetry=1``
        :param kwargs: any attribute selectors for the select function, including the following
        :param name: name of property
        :param principal: restrict to principal property
        :param preamble: XPath expression to prepend the search for the property node
        :param command: restrict to properties contained in a jobstep with this command
        :param dict: return property value, scalar float or list of floats as appropriate, if ``False`` or omitted. Otherwise, a dictionary containing all the data in the property node is returned
        :return: list of properties
        """

        query = 'property'
        preamble = '//'
        dict = False
        for v in args:
            query += v
        for k, v in kwargs.items():
            if k == 'command':
                preamble = '//jobstep[@command="' + v + '"]/'
            elif k == 'preamble':
                preamble = v
            elif k == 'dict':
                dict = v
            elif k == 'principal':
                query += '[@principal="true"]' if v else '[@principal="false"]'
            else:
                query += '[@' + k + '="' + v + '"]'
        dictarray = [element_to_dict(node) for node in self.xpath(preamble + query)]
        for e in dictarray:
            r = [float(v) for v in e['value'].split()]
            e['value'] = r[0] if len(r) == 1 else r
        if dict:
            return dictarray
        else:
            return [e['value'] for e in dictarray]

    def energies(self, *args, **kwargs):
        """
        Wrapper for :py:meth:`properties()` that restricts to energy values

        :param args:
        :param kwargs:
        :return:
        """
        return self.properties(
            "[ contains( translate(@name, 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'energy' ) ]",
            *args, **kwargs)

    def property(self, *args, **kwargs):
        """
        Wrapper for :py:meth:`properties()` that returns a single value as a scalar. An exception is thrown if more than one value is found.

        :param args:
        :param kwargs:
        :return:
        """
        result = self.properties(*args, **kwargs)
        if len(result) == 1:
            return result[0]
        elif len(result) == 0:
            return None
        else:
            raise ValueError('property() has found more than one result; use properties() instead')

    def energy(self, *args, **kwargs):
        """
        Wrapper for :py:meth:`energies()` that returns a single value as a scalar. An exception is thrown if more than one value is found.

        :param args:
        :param kwargs:
        :return:
        """
        result = self.energies(*args, **kwargs)
        if len(result) == 1:
            return result[0]
        elif len(result) == 0:
            return None
        else:
            raise ValueError('energy() has found more than one result; use energies() instead')

    def geometry(self, preamble='//', instance=-1):
        """
        Obtain the geometry from the job output

        :param preamble: Prepend this to the xpath search expression to locate the geometry
        :param instance: In the case of multiple geometries in the output stream, which occurence
        :return: dictionary holding the geometry. Coordinates are in bohr
        """
        return self.geometries(preamble=preamble)[instance]

    def geometries(self, preamble='//', parent_node=None):
        """
         Obtain all geometries from the job output

         :param preamble: Prepend this to the xpath search expression to locate the geometry
         :return: list of dictionaries holding the geometry. Coordinates are in bohr
         """
        search = preamble + "*/cml:atomArray"
        nodes = self.xpath(search, parent_node)
        Angstrom = 1.88972612462577
        geoms = []
        for node in nodes:
            atoms = []
            for atom in node.xpath("cml:atom", namespaces={'cml': "http://www.xml-cml.org/schema"}):
                atoms.append({
                    'id': atom.xpath("@id")[0],
                    'elementType': atom.xpath("@elementType")[0],
                    'xyz': [Angstrom * float(atom.xpath("@x3")[0]), Angstrom * float(atom.xpath("@y3")[0]),
                            Angstrom * float(atom.xpath("@z3")[0])]
                })
            geoms.append(atoms)
        return geoms

    def xyz(self, preamble='//', instance=-1, title=''):
        """
        Obtain a geometry in xyz file format

        :param instance: In the case of multiple geometries in the output stream, which occurence
        :param preamble: Prepend this to the xpath search expression to locate the geometry
        :param title:  A title to be injected into the second line
        :return:  A string containing the xyz representation of the geometry
        """
        geom = self.geometry(preamble=preamble, instance=instance)
        result = str(len(geom)) + '\n' + title + '\n'
        for atom in geom:
            result += atom['elementType']
            for i in range(3):
                result += ' ' + str(atom['xyz'][i])
            result += '\n'
        return result

    def xyzs(self, preamble='//', title=''):
        """
        Obtain a set of geometries in xyz file format

        :param preamble: Prepend this to the xpath search expression to locate the geometry
        :param title:  A title to be injected into the second line
        :return:  A list of strings containing the xyz representation of the geometry
        """
        geoms = self.geometries(preamble=preamble)
        result = []
        for geom in geoms:
            result.append(str(len(geom)) + '\n' + title + '\n')
            for atom in geom:
                result[-1] += atom['elementType']
                for i in range(3):
                    result[-1] += ' ' + str(atom['xyz'][i])
                result[-1] += '\n'
        return result

    def orbitals_to_trexio(self, filename=None, instance=-1, overwrite=True):
        r"""
        Create a TrexIO dump containing the geometry and orbitals.

        :param filename: Name of the Trexio file to be placed in the run directory
        :type filename: str
        :param instance: Which instance of the molecule node in the xml file
        :type instance: int
        :param overwrite: Overwrite existing TrexIO file
        :type overwrite: bool
        :return: file, label
        :rtype: str, str
        """
        import trexio
        Angstrom = 1.88972612462577
        molecule_node = self.xpath('//*/molecule')[instance]
        molecule_info = self.molecule(instance)
        orbitalSets = molecule_info['orbitalSets']
        if len(orbitalSets) < 1:
            raise Exception('No orbital sets found')

        file = self.filename('h5', filename) if filename is not None else self.filename('h5')
        if instance != -1:
            file = file.replace('.h5', '_' + str(instance) + '.h5')
        if overwrite:
            pathlib.Path(file).unlink(missing_ok=True)
        with trexio.File(file, mode='w', back_end=trexio.TREXIO_HDF5) as f:
            atoms = molecule_info['geometry']
            trexio.write_nucleus_num(f, len(atoms))
            trexio.write_nucleus_label(f, [atom['elementType'] for atom in atoms])
            trexio.write_nucleus_coord(f, [[c for c in atom['xyz']] for atom in atoms])
            trexio.write_nucleus_charge(f, [float(periodic_table.index(atom['elementType']) + 1) for atom in atoms])

            basisSets = self.xpath('basisSet[@id="ORBITAL"]', molecule_node)
            if len(basisSets) != 1:
                raise Exception('something is wrong: there should be just one orbital basisSet')
            basisSet = basisSets[0]

            prim_num = 0
            shell_num = 0
            nucleus_index = []
            shell_ang_mom = []
            shell_factor = []
            shell_index = []
            exponent = []
            coefficient = []
            prim_factor = []
            ao_shell = []
            ao_num = 0
            map_aos = []
            ao_normalization = []
            for atom_index, atom in enumerate(atoms):
                query = 'association/atoms[@xlink:href[contains(.,"@id=\'' + atom.get(
                    'id') + '\'")]]/..'
                basisGroupAssociation = self.xpath(query, basisSet)
                if len(basisGroupAssociation) != 1:
                    raise Exception(
                        'something is wrong: there should be a unique association of an atom with a basis set')
                bases = self.xpath('bases', basisGroupAssociation[0])
                if len(bases) != 1:
                    raise Exception('something is wrong: there should be a bases node in association')
                basesString = bases[0].get('{http://www.w3.org/1999/xlink}href')
                basesString = basesString[basesString.find('basisGroup['):]
                basesString = basesString[basesString.find('[') + 1:].lstrip()
                basesString = basesString[:basesString.find(']')].rstrip()
                basesString = basesString.replace('or', '').replace('\n', '').replace("'", '')
                list = basesString.split("@id=")
                for item in list:
                    item = item.lstrip().rstrip()
                    if item.isalnum():
                        basisGroup = self.xpath('basisGroup[@id="' + item + '"]', basisSet)[0]
                        lquant = int(basisGroup.get('minL'))
                        if lquant != int(basisGroup.get('maxL')):
                            raise Exception("This program cannot handle multiple-angular momentum sets")
                        alpha = np.float64(re.sub(' +', ' ',
                                                  self.xpath('basisExponents', basisGroup)[
                                                      0].text.replace('\n', '').lstrip().rstrip()).split(" "))
                        for basisContraction in self.xpath('basisContraction', basisGroup):
                            nucleus_index.append(atom_index)
                            shell_ang_mom.append(lquant)
                            shell_factor.append(1.0)
                            cc = np.float64(
                                re.sub(' +', ' ', basisContraction.text.replace('\n', '').lstrip().rstrip()).split(" "))
                            ao_num_base = int(ao_num)
                            molpro_powers = molpro_xyz_powers(lquant)
                            for i, powers in enumerate(lexical_xyz_powers(lquant)):
                                ao_shell.append(shell_num)
                                for j, mp in enumerate(molpro_powers):
                                    if powers == mp:
                                        map_aos.append(ao_num_base + j)
                                nfac = 1
                                for k in range(3):
                                    if powers[k] > 0:
                                        nfac = nfac * round(factorial2(2 * powers[k] - 1))
                                ao_normalization.append(1.0 / np.sqrt(nfac))
                                ao_num += 1
                            for i in range(len(cc)):
                                shell_index.append(shell_num)
                                exponent.append(alpha[i])
                                coefficient.append(cc[i])
                                normalization = (2.0 * alpha[i] / np.pi) ** 0.75 * np.sqrt(2.0 ** lquant) * (
                                        np.sqrt(2.0 * alpha[i]) ** lquant)
                                prim_factor.append(normalization)
                                prim_num += 1
                            shell_num += 1

            trexio.write_metadata_author_num(f,1)
            try:
                trexio.write_metadata_author(f,[pwd.getpwuid(os.getuid()).pw_gecos])
            except NameError:
                pass
            trexio.write_metadata_code_num(f,1)
            trexio.write_metadata_code(f,["pymolpro "+pymolpro.__version__])

            trexio.write_basis_type(f, "Gaussian")
            trexio.write_basis_prim_num(f, int(prim_num))
            trexio.write_basis_shell_num(f, int(shell_num))
            trexio.write_basis_nucleus_index(f, nucleus_index)
            trexio.write_basis_shell_ang_mom(f, shell_ang_mom)
            trexio.write_basis_shell_factor(f, shell_factor)
            trexio.write_basis_shell_index(f, shell_index)
            trexio.write_basis_r_power(f, [0 for _ in range(shell_num)])
            trexio.write_basis_exponent(f, exponent)
            trexio.write_basis_coefficient(f, coefficient)
            trexio.write_basis_prim_factor(f, prim_factor)
            trexio.write_ao_cartesian(f, 1)
            trexio.write_ao_num(f, int(ao_num))
            trexio.write_ao_shell(f, ao_shell)
            trexio.write_ao_normalization(f, ao_normalization)

            mo_num = 0
            mo_energies = []
            mo_occupations = []
            mo_spins = []
            mo_coefficients = []
            for orbital_instance in range(
                    min(2, len(orbitalSets))):  # forces that we just get the UHF alpha and beta, not natural
                spin = orbitalSets[orbital_instance]['spin']
                orbitals = orbitalSets[orbital_instance]['orbitals']
                for orb in orbitals:
                    mo_num += 1
                    if hasattr(orb, 'energy'): mo_energies.append(orb.energy)
                    mo_spins.append(1 if spin == 'beta' else 0)
                    if hasattr(orb, 'occupation'):
                        mo_occupations.append(orb.occupation)
                    mo_coefficients += [orb.coefficients[map_aos[i]] for i, c in enumerate(orb.coefficients)]
            trexio.write_mo_num(f, mo_num)
            trexio.write_mo_spin(f, mo_spins)
            if len(mo_energies) > 0: trexio.write_mo_energy(f, mo_energies)
            if len(mo_occupations) > 0: trexio.write_mo_occupation(f, mo_occupations)
            trexio.write_mo_coefficient(f, mo_coefficients)
            label = molecule_info['orbitalSets'][0]['method']
            if len(orbitalSets) == 1: label += '/' + molecule_info['orbitalSets'][0]['type']
            trexio.write_mo_type(f, label)

            nelec = int(orbitalSets[orbital_instance]['state_nelec'])
            ms2 = int(orbitalSets[orbital_instance]['state_ms2'])
            trexio.write_electron_num(f, nelec)
            trexio.write_electron_up_num(f, (nelec + ms2) // 2)
            trexio.write_electron_dn_num(f, (nelec - ms2) // 2)

            f.close()
        return file, label

    def orbitals_to_molden(self, filename=None, instance=-1, minocc=1.0, ID=None):
        Angstrom = 1.88972612462577
        molecule_node = self.xpath('//*/molecule')[instance]
        molecule_info = self.molecule(instance)
        orbitalSets = molecule_info['orbitalSets']
        if len(orbitalSets) < 1:
            raise Exception('No orbital sets found')

        file = self.filename('molden', filename) if filename is not None else self.filename('molden')
        if instance != -1:
            file = file.replace('.molden', '_' + str(instance) + '.molden')
        with open(file, 'w') as f:
            f.write('[Molden Format]\n')
            f.write('[Atoms] Angs\n')
            atoms = molecule_info['geometry']
            atomindex = 0
            for atom in atoms:
                atomindex += 1
                f.write(atom['elementType'] + ' ' + str(atomindex) + ' ' + str(
                    periodic_table.index(atom['elementType']) + 1) + ' ' + ' '.join(
                    [str(c / Angstrom) for c in atom['xyz']]) + '\n')

            f.write('[GTO]\n')

            basisSets = self.xpath('basisSet[@id="ORBITAL"]', molecule_node)
            if len(basisSets) != 1:
                raise Exception('something is wrong: there should be just one orbital basisSet')
            basisSet = basisSets[0]

            count = 0
            for atom in atoms:
                count += 1
                f.write(str(count) + ' 0\n')
                query = 'association/atoms[@xlink:href[contains(.,"@id=\'' + atom.get(
                    'id') + '\'")]]/..'
                basisGroupAssociation = self.xpath(query, basisSet)
                if len(basisGroupAssociation) != 1:
                    raise Exception(
                        'something is wrong: there should be a unique association of an atom with a basis set')
                bases = self.xpath('bases', basisGroupAssociation[0])
                if len(bases) != 1:
                    raise Exception('something is wrong: there should be a bases node in association')
                basesString = bases[0].get('{http://www.w3.org/1999/xlink}href')
                basesString = basesString[basesString.find('basisGroup['):]
                basesString = basesString[basesString.find('[') + 1:].lstrip()
                basesString = basesString[:basesString.find(']')].rstrip()
                basesString = basesString.replace('or', '').replace('\n', '').replace("'", '')
                list = basesString.split("@id=")
                for item in list:
                    item = item.lstrip().rstrip()
                    if item.isalnum():
                        basisGroup = self.xpath('basisGroup[@id="' + item + '"]', basisSet)[0]
                        lquant = int(basisGroup.get('minL'))
                        if lquant > 5:
                            raise Exception("Sorry, I was too lazy to write this for i basis functions and higher")
                        if lquant != int(basisGroup.get('maxL')):
                            raise Exception("This program cannot handle multiple-angular momentum sets")
                        alpha = np.float64(re.sub(' +', ' ',
                                                  self.xpath('basisExponents', basisGroup)[
                                                      0].text.replace('\n', '').lstrip().rstrip()).split(" "))
                        for basisContraction in self.xpath('basisContraction', basisGroup):
                            cc = np.float64(
                                re.sub(' +', ' ', basisContraction.text.replace('\n', '').lstrip().rstrip()).split(" "))
                            f.write('spdfghiklmnopqrst'[lquant] + ' ' + str(len(cc)) + ' 1.00\n')
                            for i in range(len(cc)):
                                f.write(str(alpha[i]) + ' ' + str(cc[i]) + '\n')
                f.write('\n')

            for orbital_instance in range(
                    min(2, len(orbitalSets))):  # forces that we just get the UHF alpha and beta, not natural
                spin = orbitalSets[orbital_instance]['spin']
                orbitals = orbitalSets[orbital_instance]['orbitals']
                for orb in orbitals:
                    f.write('[MO]\nSym=' + orb.ID + '\n')
                    if hasattr(orb, 'energy'): f.write('Ene= ' + str(orb.energy) + '\n')
                    f.write('Spin= ' + spin + '\n')
                    if hasattr(orb, 'occupation'): f.write('Occup= ' + str(orb.occupation) + '\n')
                    count = 0
                    for c in orb.coefficients:
                        count += 1
                        f.write(str(count) + ' ' + str(c) + '\n')
        label = molecule_info['orbitalSets'][0]['method']
        if len(orbitalSets) == 1: label += '/' + molecule_info['orbitalSets'][0]['type']
        return file, label

    def orbitals(self, instance=-1, minocc=1.0, ID=None, orbital_instance=-1):
        """
        Obtain some or all of the orbitals in the job output

        :param instance: Which molecule object to get orbitals for
        :param minocc: Only orbitals with at least this occupation will be returned
        :param ID: Only the orbital whose ID attribute matches this will be selected
        :param orbital_instance: Which set of orbitals
        :return: a list of Orbital objects
        """
        # molecule = self.molecule(instance=instance)
        # if 'orbitals' not in molecule:
        #     raise Exception('No orbital set found')
        return self.orbitals_old(instance=instance, ID=ID, minocc=minocc, orbital_instance=orbital_instance)

    def molecule(self, instance=-1):
        """
        Obtain a molecule in the job output

        :param instance: Which molecule object to get
        :return: a dictionary containing all the information available
        """
        from pymolpro import Orbital
        molecule = {}
        molecule_node = self.xpath('//*/molecule')[instance]
        molecule['id'] = molecule_node.get('id')
        molecule['method'] = molecule_node.get('method')
        if molecule_node.get('energy') not in ['', None]:
            molecule['energy'] = molecule_node.get('energy')
        orbitalSets = self.xpath('orbitals', molecule_node)
        molecule['orbitalSets'] = []
        for orbitalSet in orbitalSets:
            molecule['orbitalSets'].append({})
            for attribute in orbitalSet.keys():
                molecule['orbitalSets'][-1][attribute] = orbitalSet.get(attribute)
            molecule['orbitalSets'][-1]['orbitals'] = []
            for orbital in self.xpath('orbital', orbitalSet):
                molecule['orbitalSets'][-1]['orbitals'].append(Orbital(orbital, self.filename()))
        geometry_nodes = self.xpath('//*/molecule[@id="' + molecule['id'] + '"]/cml:molecule/cml:atomArray',
                                    molecule_node)
        if len(geometry_nodes) != 1:
            raise IndexError('Molecule node does not contain exactly one geometry node')
        Angstrom = 1.88972612462577
        atoms = []
        for atom in geometry_nodes[0].xpath("cml:atom", namespaces={'cml': "http://www.xml-cml.org/schema"}):
            atoms.append({
                'id': atom.xpath("@id")[0],
                'elementType': atom.xpath("@elementType")[0],
                'xyz': [Angstrom * float(atom.xpath("@x3")[0]), Angstrom * float(atom.xpath("@y3")[0]),
                        Angstrom * float(atom.xpath("@z3")[0])]
            })
        molecule['geometry'] = atoms
        return molecule

    def orbitals_old(self, instance=-1, minocc=1.0, ID=None, orbital_instance=-1):
        """
        Obtain some or all of the orbitals in the job output

        :param instance: Which molecule object to get orbitals for
        :param minocc: Only orbitals with at least this occupation will be returned
        :param ID: Only the orbital whose ID attribute matches this will be selected
        :param orbital_instance: Which set of orbitals
        :return: a list of Orbital objects
        """
        molecule = self.xpath('//*/molecule')[instance]
        orbitalSets = self.xpath('orbitals', molecule)
        result = []
        search = 'orbital'
        if ID:
            search += '[@ID="' + ID + '"]'
        from pymolpro import Orbital
        for orbital in self.xpath(search, orbitalSets[orbital_instance]):
            if float(orbital.get('occupation')) >= minocc:
                result.append(Orbital(orbital, self.filename()))
        assert not ID or len(result) == 1
        return result

    def evaluateOrbitals(self, points, instance=-1, minocc=1.0e-10, ID=None, values=False):
        """
        Evaluate molecular orbitals on a grid of points

        :param points: List of geometries specified as 3-list of values in bohr, or numpy [:,3] array
        :param instance: Which set of orbitals
        :param minocc: Only orbitals with at least this occupation will be returned
        :param ID: Only the orbital whose ID attribute matches this will be selected
        :param values:
        :return: array of dictionaries giving the occupation and values on the grid, or if ID is specified, a single dictionary, or if values==True, a numpy array
        """
        if type(points) is list:
            return self.evaluateOrbitals(np.array(points), instance=instance, minocc=minocc,
                                         ID=ID, values=values)
        molecule = self.xpath('//*/molecule')[instance]
        import pymolpro.grid
        return pymolpro.grid.evaluateOrbitals(molecule, points, minocc=minocc, ID=ID, values=values, directory=self.filename())

    def pairs(self, instance=-1):
        """
        Obtain some or all of the correlation pairs in the job output

        :param instance: Which set of pairs
        :return: a list of Pair objects
        """
        # jobsteps = self.xpath('//*/jobstep[pair]') #TODO use this instead of following when pysjef.xpath() supports default namespace in []
        jobsteps = self.xpath('//*/jobstep[@commandset="CCSD"]')
        if len(jobsteps) == 0:
            raise Exception('No orbital pairs found')
        if (instance >= 0 and len(jobsteps) <= instance) or len(jobsteps) < abs(instance):
            raise Exception('Not enough pair-containing jobsteps found')
        from pymolpro import Pair
        return [Pair(pair) for pair in self.xpath('pair', jobsteps[instance])]

    def singles(self, instance=-1):
        """
        Obtain some or all of the correlation singles in the job output

        :param instance: Which set of singles
        :return: a list of Single objects
        """
        # jobsteps = self.xpath('//*/jobstep[single]') #TODO use this instead of following when pysjef.xpath() supports default namespace in []
        jobsteps = self.xpath('//*/jobstep[@commandset="CCSD"]')
        if len(jobsteps) == 0:
            raise Exception('No orbital singles found')
        if (instance >= 0 and len(jobsteps) <= instance) or len(jobsteps) < abs(instance):
            raise Exception('Not enough single-containing jobsteps found')
        from pymolpro import Single
        return [Single(single) for single in self.xpath('single', jobsteps[instance])]

    def variable(self, name, instance=-1, list=False, dict=False):
        """
        Return the value of a variable in the output xml stream

        :param name: The name of the variable
        :param instance: index of occurence in output
        :param list: Whether to force returning a list. If true, a list is always returned; otherwise if the result is a scalar, a scalar is returned, and if no match is found, `None` is returned
        :param dict: return property value, scalar float or list of floats as appropriate, if ``False`` or omitted. Otherwise, a dictionary containing all the data in the property node is returned
        :return:
        """
        matches = self.xpath('//variables/variable[@name="' + name.upper() + '"]')
        if len(matches) == 0:
            matches = self.xpath('//variables/variable[@name="_' + name.upper() + '"]')
        if len(matches) == 0:
            matches = self.xpath('//variables/variable[@name="!' + name.upper() + '"]')
        if len(matches) == 0 or len(matches) <= instance or len(matches) < abs(instance):
            if list:
                return []
            return None
        value_nodes = matches[instance].xpath('molpro-output:value', namespaces={
            'molpro-output': 'http://www.molpro.net/schema/molpro-output'})
        values = [node.text for node in value_nodes]
        type = 'xsd:string'
        if len(matches[instance].xpath('@type')) > 0:
            type = matches[instance].xpath('@type')[0]
        length = 1
        if len(matches[instance].xpath('@length')) > 0:
            length = int(matches[instance].xpath('@length')[0])
        if type == 'xsd:double':
            values = [float(k) for k in values] if list or length > 1 else float(values[0])
        if dict:
            result = element_to_dict(matches[instance])
            del result['type']
            del result['length']
            result['value'] = values
            return result
        else:
            return values

    def variables(self, instance=-1):
        """
        Return a list of all defined variables

        :param instance: index of occurence of variables node in output
        :return:
        """
        matches = self.xpath('//variables')
        if len(matches) == 0:
            return []
        return [node.xpath('@name')[0] for node in (matches[instance].xpath('molpro-output:variable', namespaces={
            'molpro-output': 'http://www.molpro.net/schema/molpro-output'}))]

    def input_from_output(self, instance=-1):
        """
        Return a list of all defined input in the output stream

        :param instance: index of occurence of input node in output
        :return:
        """
        matches = self.xpath('//input')
        if len(matches) == 0:
            return []
        instance_ = matches[instance]
        return [node.text for node in (instance_.xpath('molpro-output:p', namespaces={
            'molpro-output': 'http://www.molpro.net/schema/molpro-output'}))]

    def run_local_molpro(self, options: list):
        logger.debug(f"run_local_molpro with options: {options}")
        logger.debug(f"PATH: {os.environ.get('PATH')}")
        return subprocess.run(['molpro', *options], capture_output=True, text=True)
        def get_first(string):
            result = ''
            quoted = False
            for i in range(len(string)):
                if string[i] == "'":
                    quoted = not quoted
                elif string[i] == ' ' and not quoted:
                    return result
                else:
                    result = result + string[i]

        command_ = (self.backend_get('local', 'run_command')).strip(' ')
        command_ = re.sub('mpiexec', 'mpiexec -n 1', re.sub('  *', ' ', re.sub('{.*?}', '', command_))).strip(' ')
        split_options = command_.split(' ')[1:] + options
        run = subprocess.run([shutil.which(command_.split(' ')[0])] + split_options, capture_output=True, shell=False)
        return run

    def registry(self, set=None):
        """
        Access the registry belonging to Molpro in the `local` backend.

        :param set: If present, obtain as a dictionary the specified registry set. If absent, obtain a list of all sets
        :return:
        """
        if set:
            if not hasattr(self, 'registry_cache'):
                self.registry_cache = {}
            if set not in self.registry_cache:
                try:
                    run = self.run_local_molpro(['--registry', set])
                    if not run.stdout:
                        return None
                    l0 = str(run.stdout)
                    l1 = l0[l0.find(':\\n') + 1:].rstrip("'").replace('\\n', '')
                    # print('l1',l1)
                    line = l1.replace('{', '').strip('\\n').strip('"').split('}')
                    # print('l',l)
                    self.registry_cache[set] = {}
                    for li in line:
                        # print('li',li)
                        name = re.sub('["\'].*', '', re.sub('.*name *= *["\']', '', li))
                        if name:
                            json_ = '{' + re.sub('"type" *: * (.),', '"type": "\\1",',
                                                 re.sub('([*a-zA-Z0-9_-]+)=', '"\\1": ', li)).replace("'", '"') + '}'
                            # print(json_)
                            self.registry_cache[set][name] = json.loads(json_)
                            if 'name' in self.registry_cache[set][name]:
                                del self.registry_cache[set][name]['name']
                except Exception:
                    if set == 'DFUNC':  # to keep unit testing quiet when there's no molpro
                        self.registry_cache[set] = {'B3LYP': {'set': 'DFUNC'}}
                    else:
                        return None
            return self.registry_cache[set]
        else:
            try:
                run = self.run_local_molpro(['--registry'])
                if not run.stdout:
                    return None
                return re.sub('.*: *', '', str(run.stdout)).replace('\\n', '').rstrip("'").split()
            except Exception:
                return None

    import builtins
    @builtins.property
    def local_molpro_root(self):
        r"""
        Get the directory of the Molpro installation in the local backend

        :return: directory
        :rtype: pathlib.Path
        """
        if self.local_molpro_root_ is None:
            logger.debug('setting local_molpro_root')
            try:
                run = self.run_local_molpro(['--registry'])
                logger.debug(f"run.stdout: {run.stdout}")
                logger.debug(f"run.stderr: {run.stderr}")
                logger.debug(f"run.returncode: {run.returncode}")
                if run.stdout:
                    logger.debug('Project.local_molpro_root makes stdout')
                    out1 = str(run.stdout).split('\n')[0]
                    path_to_lib = pathlib.Path(
                        re.sub(r'\\n.*', '', re.sub('.*registry at *', '', out1)).rstrip("'").replace('\\n',
                                                                                                                 ''))
                    self.local_molpro_root_ = path_to_lib.parent
                    logger.debug('Project.local_molpro_root sets local_molpro_root to '+str(self.local_molpro_root_))
            except Exception:
                return None
        return self.local_molpro_root_

    def procedures_registry(self):
        r"""
        Get the procedures registry from the Molpro pointed to in the local backend

        :return:
        :rtype: dict
        """
        entries = {}
        logger.debug('pymolpro.Project.procedures_registry()')
        try:
            all = ''
            with open(self.local_molpro_root / 'lib' / 'procedures.registry', 'r') as f:
                for line in f.readlines():
                    logger.debug('procedures.registry line: '+line)
                    if not re.match('^ *#', line) and line.strip(' ') != '':
                        all += line.replace('}', '').replace('\n', '')
            for record in all.split('{'):
                if not record:
                    continue
                in_quotations = False
                for i in range(len(record)):
                    if record[i] == "'":
                        in_quotations = not in_quotations
                    elif record[i] == ',' and in_quotations:
                        record = record[:i] + '§' + record[i + 1:]
                fields = record.split(',')
                entry = {}
                name = None
                for field in fields:
                    match = re.match(r"([^ =]+) *= *'?([^']+)'?", field)
                    if match:
                        entry[match.group(1)] = int(match.group(2)) if re.search('^-?[0-9]+$', match.group(2).strip(
                            ' ')) else match.group(2)
                        if match.group(1) == 'name':
                            name = match.group(2)
                if 'options' in entry:
                    entry['options'] = entry['options'].split('§')
                entries[name] = entry
        except Exception:
            pass
        return entries

    def basis_registry(self):
        r"""
        Get the basis registry from the Molpro pointed to in the local backend

        :return:
        :rtype: dict
        """
        entries = {}
        try:
            with open(self.local_molpro_root / 'lib' / 'basis.registry', 'r') as f:
                for line in f.readlines():
                    if not re.match('^ *#', line) and line.strip(' ') != '':
                        fields = line.split(',')
                        entry = {}
                        name = None
                        for field in fields:
                            match = re.match(r"([^ =]+) *= *'?([^']+)'?", field)
                            if match:
                                entry[match.group(1)] = int(match.group(2)) if re.search('^-?[0-9]+$',
                                                                                         match.group(2).strip(
                                                                                             ' ')) else match.group(2)
                                if match.group(1) == 'name':
                                    name = match.group(2)
                        if 'options' in entry:
                            entry['options'] = entry['options'].split(':')
                        if name is not None:
                            entries[name] = entry
        except Exception:
            pass
        return entries

    def gradient(self, instance=-1):
        try:
            grad = np.array([[float(x) for x in line.split()] for line in
                             self.xpath_search('//gradient')[instance].strip().split('\n')])
            return grad.reshape([grad.size])
        except:
            return None

    @builtins.property
    def vibrations(self, instance=-1):
        r"""
        Give information on normal vibrational modes

        :param instance:
        :type instance: int
        :return: Mass-weighted normal coordinates, vibrational wavenumbers, atomic masses, force constant matrix
        :rtype: dict
        """
        result = {}
        vibrations = self.xpath('//vibrations')
        if vibrations is None or len(vibrations) == 0:
            return None
        vibration_set = vibrations[instance]
        nodes = self.xpath('normalCoordinate', vibration_set)
        result['wavenumbers'] = np.array([float(k.get('wavenumber')) for k in nodes])
        if len(nodes) > 0 and nodes[0].get('IRintensity') is not None:
            result['IRintensity'] = np.array([float(k.get('IRintensity')) for k in nodes])
        normal_coordinates = np.array(
            [[float(x) for x in ''.join(node.text.strip().splitlines()).split()] for node in nodes])
        result['normal_coordinates'] = normal_coordinates
        result['atomic_masses'] = np.array(
            [float(x) for x in ''.join(self.xpath('masses', vibration_set)[0].text.strip().splitlines()).split()])
        sqrt_mass_matrix = np.diag(
            np.array([math.sqrt(mass * 1822.88848621731) for mass in result['atomic_masses'] for i in range(3)]))
        result['energies'] = result['wavenumbers'] / 219474.63
        mass_weighted_normal_coordinates = normal_coordinates @ sqrt_mass_matrix
        for i in range(len(mass_weighted_normal_coordinates)):
            mass_weighted_normal_coordinates[i, :] = mass_weighted_normal_coordinates[i, :] / np.linalg.norm(
                mass_weighted_normal_coordinates[i, :])
        result['mass_weighted_normal_coordinates'] = mass_weighted_normal_coordinates
        result['real_zero_imag'] = np.array(
            [0.0 if k.get('real_zero_imag') == 'Z' else -1.0 if k.get('real_zero_imag') == 'I' else 1.0 for k in nodes])
        result['force_constants'] = sqrt_mass_matrix @ mass_weighted_normal_coordinates.transpose() @ np.diag(
            result['energies']) @ np.diag(result['real_zero_imag']) @ np.diag(
            result['energies']) @ mass_weighted_normal_coordinates @ sqrt_mass_matrix
        return result

    def _anonymous_name(self, input: str | dict | None = None, specification: str | dict | None = None,
                        ansatz: str = None, **kwargs) -> str:
        import hashlib
        project_name = hashlib.sha256((str(input) + str(specification) + str(ansatz) +
                                       str(tuple(sorted(kwargs.items())))).encode(
            'utf-8')).hexdigest()[-8:]
        return project_name

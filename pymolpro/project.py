import math
import os
import pathlib

import pysjef
import subprocess
import re
import json
import numpy as np
import shutil


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

    :param str name: The base filename of the filesystem bundle carrying the project. If the bundle does not yet exist, it is created.
    :param str geometry: The geometry.
        If specified, the input for the job will be constructed, without the need for a subsequent call to :py:meth:`write_input()`.
        Any format recognised by Molpro can be used. This includes xyz, with or without the two header lines, or Z matrix, and lines can be separated either with newline or `;`. The geometry can be specified either as a string, or a filename or url reference, in which case the contents of the reference are resolved now.
    :param str method: The computational method for constructed input. Anything accepted as Molpro input, including parameters and directives, can be given.  If the method needs a preceding Hartree-Fock calculation, this is prepended automatically.
    :param str basis: The orbital basis set for constructed input. Anything that can appear after `basis=` in Molpro input is accepted.
    :param str ansatz: String of the form method/basis//geometry_method/geometry_basis or method/basis which is parsed to give the same effect as the
    method and basis parameters. If geometry_method/geometry_basis is specified, the calculation will be preceded by a geometry optimisation at that level of theory.
    :param str func: This should be one of

        * `energy` for a single geometry
        * `opt` for a geometry optimisation

    :param str extrapolate: If specified, carry out basis-set extrapolation. Anything that can appear after `extrapolate,basis=` in Molpro input is accepted.
    :param int spin: The spin multiplicity minus one
    :param int charge: Electrical charge of molecule

    """

    def __init__(self, name=None, geometry=None, method="hf", basis="cc-pVTZ", func="energy", extrapolate="", symm=True,
                 ansatz=None,
                 geometry_method=None,
                 geometry_basis=None,
                 preamble=None,
                 postamble=None,
                 initial=None,
                 charge=None,
                 spin=None,
                 files=[],
                 **kwargs):
        if not hasattr(self,'__initialized'):
            try:
                super().__init__(name=name, **kwargs)
                self.__initialized = True
            except Exception:
                raise FileNotFoundError("Cannot open project " + name)

        if ansatz is not None:
            self.__init__(name=name, geometry=geometry,
                                method=self.parse_ansatz(ansatz)['method'], basis=self.parse_ansatz(ansatz)['basis'], func=func,
                                extrapolate=extrapolate, symm=symm,
                                geometry_method=self.parse_ansatz(ansatz)['geometry_method'],
                                geometry_basis=self.parse_ansatz(ansatz)['geometry_basis'],
                                preamble=preamble, postamble=postamble, initial=initial, charge=charge, spin=spin,
                                **kwargs
                                )
            return
        self.local_molpro_root_ = None

        self.ansatz = None
        if geometry is not None:  # construct input
            # print('method',method,':',self.commandify_method(method),':',method_from_commands(method))
            # print('basis',basis)
            self.ansatz = (method_from_commands(method) + '/' + basis.strip()).strip()
            # print('self.ansatz', self.ansatz, 'ansatz', ansatz)
            if geometry_method is not None and geometry_basis is not None:
                self.ansatz += '//' + method_from_commands(geometry_method) + '/' + geometry_basis.strip()
            input = ''
            if initial is not None: input += initial + '\n'
            if not symm: input += 'symmetry, nosym\n'
            input += 'geometry={' + resolve_geometry(geometry) + '}\n'
            if preamble is not None: input +=  preamble + '\n'
            if charge is not None: input += 'charge=' + str(charge) + '\n'
            if spin is not None: input += 'spin=' + str(spin) + '\n'
            if geometry_method is not None:
                input += 'basis=' + geometry_basis + ';' + self.commandify_method(geometry_method) + '; optg\n'
            input += 'basis=' + basis + '\n'
            input += self.commandify_method(method) + '\n'
            if extrapolate != '': input += 'extrapolate,basis=' + extrapolate + '\n'
            if func[:3] == 'opt': input += 'optg\n'
            if postamble is not None: input += str(postamble) + '\n'
            if postamble is None or 'put,xml' not in postamble: input += '{put,xml;noorbitals,nobasis}\n'
            self.write_input(input)
        if files is not None and len(files) > 0:
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
                        shutil.copyfile(file,self.filename(path.suffix[1:]))
                        if path.suffix == '.xml':
                            input_from_output = self.input_from_output()
                            self.write_input('\n'.join(input_from_output))
                    elif path.suffix in ['.inp','.xyz'] and path.stem != 'optimised':
                        shutil.copyfile(file,pathlib.Path(self.filename('',run=-1))/path.name)
                    else:
                        if not rundir:
                            self.run_directory_new()
                            rundir = True
                        shutil.copyfile(file,pathlib.Path(self.filename(''))/path.name)

    # def set_method(self,method,basis="cc-pVTZ",geometry_method=None, geometry_basis=None):
    #     pass
    def commandify_method(self,method):
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

    def parse_ansatz(self,ansatz):
        parsed={}
        parts = ansatz.split('//')
        # if len(parts) == 1:
        #     parts.append(parts[0])
        result = []
        for part in parts:
            bits = part.split('/')
            # print('part',part,'bits',bits)
            if bits[0].strip().lower() == 'hf':
                bits[0] = 'df-hf'
            # elif bits[0].strip().lower() == 'mp2':
            #     bits[0] = 'df-hf;df-mp2'
            elif bits[0].strip().upper() in self.registry('DFUNC').keys():
                bits[0] = 'df-ks,' + bits[0]
            elif bits[0].strip().lower()[:3] == 'df-' and not re.match(bits[0].strip().lower()[3:], '[ur]?hf'):
                bits[0] = 'df-hf; ' + bits[0]
            else:
                bits[0] = 'hf; ' + bits[0]
            if len(bits) == 1:
                bits.append('cc-pVDZ')
            # print('bits',bits)
            result += bits
        parsed['ansatz'] = ansatz
        parsed['method'] = result[0]
        parsed['basis'] = result[1]
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

    def geometries(self, preamble='//'):
        """
         Obtain all geometries from the job output

         :param preamble: Prepend this to the xpath search expression to locate the geometry
         :return: list of dictionaries holding the geometry. Coordinates are in bohr
         """
        search = preamble + "*/cml:atomArray"
        nodes = self.xpath(search)
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

    def orbitals(self, instance=-1, minocc=1.0, ID=None):
        """
        Obtain some or all of the orbitals in the job output

        :param instance: Which set of orbitals
        :param minocc: Only orbitals with at least this occupation will be returned
        :param ID: Only the orbital whose ID attribute matches this will be selected
        :return: a list of Orbital objects
        """
        molecule = self.xpath('//*/molecule')[instance]
        orbitalSets = self.xpath('orbitals', molecule)
        if len(orbitalSets) != 1:
            raise Exception('something is wrong: there should be just one orbital set')
        result = []
        search = 'orbital'
        if ID:
            search += '[@ID="' + ID + '"]'
        from pymolpro import Orbital
        for orbital in self.xpath(search, orbitalSets[0]):
            if float(orbital.get('occupation')) >= minocc:
                result.append(Orbital(orbital))
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
        return pymolpro.grid.evaluateOrbitals(molecule, points, minocc=minocc, ID=ID, values=values)

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
            try:
                run = self.run_local_molpro(['--registry'])
                if run.stdout:
                    self.local_molpro_root_ = pathlib.Path(
                        re.sub(r'\\n.*', '', re.sub('.*registry at *', '', str(run.stdout))).rstrip("'").replace('\\n',
                                                                                                                 '')).parent
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
        try:
            all = ''
            with open(self.local_molpro_root / 'lib' / 'procedures.registry', 'r') as f:
                for line in f.readlines():
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
            np.array([math.sqrt(mass* 1822.88848621731) for mass in result['atomic_masses'] for i in range(3)]))
        result['energies'] = result['wavenumbers'] / 219474.63
        mass_weighted_normal_coordinates = normal_coordinates @ sqrt_mass_matrix
        for i in range(len(mass_weighted_normal_coordinates)):
            mass_weighted_normal_coordinates[i, :] = mass_weighted_normal_coordinates[i, :] / np.linalg.norm(
                mass_weighted_normal_coordinates[i, :])
        result['mass_weighted_normal_coordinates'] = mass_weighted_normal_coordinates
        result['real_zero_imag'] = np.array([0.0 if k.get('real_zero_imag')=='Z' else -1.0 if k.get('real_zero_imag')=='I' else 1.0 for k in nodes])
        result['force_constants'] = sqrt_mass_matrix @ mass_weighted_normal_coordinates.transpose() @ np.diag(
            result['energies']) @ np.diag(result['real_zero_imag']) @ np.diag(result['energies']) @ mass_weighted_normal_coordinates @ sqrt_mass_matrix
        return result

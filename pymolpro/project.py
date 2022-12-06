import pysjef
from pysjef.select import select


def no_errors(projects, ignore_warning=True):
    """
    Checks that none of the projects have any errors. Projects can by running.

    :param projects: list of projects to check for errors
    :param ignore_warning: Whether to count warnings as well as errors.
    :return: True/False whether any projects have errors
    """
    for p in projects:
        p.parse()
    errors = select(projects, '//error')
    if ignore_warning:
        errors = [x for x in errors if
                  not ('type' in x.attributes and x.attributes['type'].lower() == 'warning')]
    return len(errors) == 0


def element_to_dict(node, attributes=True):
    """
    Convert an lxml.etree node tree into a dict.
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


class Project(pysjef.project.Project):
    """
    Python binding to sjef, for managing molpro jobs.
    Project is a node with parsed molpro output as the only child.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def errors(self, ignore_warning=True):
        '''
        Return all error nodes

        :return: list of error nodes
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
        if principal: string += f', principal=true'
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
        return self.properties('[@name="Energy" or @name="total energy"]', *args, **kwargs)

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


    def evaluateOrbitals(self, points, instance=-1, minocc=1.0, ID=None, values=False):
        """
        Evaluate molecular orbitals on a grid of points

        :param points: List of geometries specified as 3-list of values in bohr
        :param instance: Which set of orbitals
        :param minocc: Only orbitals with at least this occupation will be returned
        :param ID: Only the orbital whose ID attribute matches this will be selected
        :param values:
        :return: array of dictionaries giving the occupation and values on the grid, or if ID is specified, a single dictionary, or if values==True, a numpy array
        """
        molecule = self.xpath('//*/molecule')[instance]
        import pymolpro.grid
        return pymolpro.grid.evaluateOrbitals(molecule, points, minocc=minocc, ID=ID, values=values)

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
        if len(matches[instance].xpath('@type')) > 0: type = matches[instance].xpath('@type')[0]
        length = 1
        if len(matches[instance].xpath('@length')) > 0: length = int(matches[instance].xpath('@length')[0])
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
        if len(matches) == 0: return []
        return [node.xpath('@name')[0] for node in (matches[instance].xpath('molpro-output:variable', namespaces={
            'molpro-output': 'http://www.molpro.net/schema/molpro-output'}))]

    def input_from_output(self, instance=-1):
        """
        Return a list of all defined input in the output stream

        :param instance: index of occurence of input node in output
        :return:
        """
        matches = self.xpath('//input')
        if len(matches) == 0: return []
        instance_ = matches[instance]
        return [node.text for node in (instance_.xpath('molpro-output:p', namespaces={
            'molpro-output': 'http://www.molpro.net/schema/molpro-output'}))]

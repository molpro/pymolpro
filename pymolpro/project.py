import pysjef
from pysjef.select import select


def no_errors(projects, ignore_warning=True):
    """
    Checks that none of the projects have any errors. Projects can by running.

    :param projects: list of projects to check for errors
    :return: True/False whether any projects have errors
    """
    for p in projects:
        p.parse()
    errors = select(projects, '//error')
    if ignore_warning:
        errors = [x for x in errors if
                  not ('type' in x.attributes and x.attributes['type'].lower() == 'warning')]
    return len(errors) == 0


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

    def properties(self, name='Energy', principal=True, *, value=False, **kwargs, ):
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

    def geometry(self, instance=-1):
        """
        Obtain the geometry from the job output

        :param instance: In the case of multiple geometries in the output stream, which occurence
        :return: dictionary holding the geometry. Coordinates are in bohr
        """
        return self.geometries()[instance]

    def geometries(self):
        """
         Obtain all geometries from the job output

         :return: list of dictionaries holding the geometry. Coordinates are in bohr
         """
        nodes = self.xpath("//*/cml:atomArray")
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

    def evaluateOrbitals(self, points, instance=-1, minocc=1.0):
        """
        Evaluate molecular orbitals on a grid of points

        :param points: List of geometries specified as 3-list of values in bohr
        :param instance: Which set of orbitals
        :param minocc: Only orbitals with at least this occupation will be returned
        :return:
        """
        molecule = self.xpath('//*/molecule')[instance]
        import pymolpro.grid
        return pymolpro.grid.evaluateOrbitals(molecule, points, minocc)

    def variable(self, name, instance=-1, list=False):
        """
        Return the value of a variable in the output xml stream

        :param name: The name of the variable
        :param instance: index of occurence in output
        :param list: Whether to force returning a list. If true, a list is always returned; otherwise
        if the result is a scalar, a scalar is returned, and if no match is found, None is returned
        :return:
        """
        matches = self.xpath('//variables/variable[@name="' + name + '"]')
        if len(matches) == 0:
            matches = self.xpath('//variables/variable[@name="_' + name + '"]')
        if len(matches) == 0:
            matches = self.xpath('//variables/variable[@name="!' + name + '"]')
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
        if list or length > 1:
            if type == 'xsd:double':
                return [float(k) for k in values]
            return values
        else:
            if type == 'xsd:double':
                return float(values[0])
            return values[0]

    def variable_units(self, name, instance=-1):
        """
        Return the value of the units of a variable in the output xml stream

        :param name:  The name of the variable
        :param instance: index of occurence in output
        :return:
        """
        matches = self.xpath('//variables/variable[@name="' + name + '"]')
        if len(matches) == 0:
            matches = self.xpath('//variables/variable[@name="_' + name + '"]')
        if len(matches) == 0:
            matches = self.xpath('//variables/variable[@name="!' + name + '"]')
        if len(matches) == 0 or len(matches) <= instance or len(matches) < abs(instance):
            return None
        units_nodes = matches[instance].xpath('@units')
        if len(units_nodes) == 0: return None
        return units_nodes[0]

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

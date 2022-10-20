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
                    'xyz': [Angstrom*float(atom.xpath("@x3")[0]), Angstrom*float(atom.xpath("@y3")[0]), Angstrom*float(atom.xpath("@z3")[0])]
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

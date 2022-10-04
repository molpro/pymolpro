from pysjef.project import Project
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


class MProject(Project):
    """
    Python binding to molpro-project, for managing molpro jobs.
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
        Shorthand for `p.select('//properties[name= {}, principal ={}, ... ]')`

        :param name: name of property
        :param principal: principal property
        :param value: return by value
        :param kwargs: any other attribute selectors for the select function
        :return: list of properties
        """
        string = f'name={name}, principal={principal}'
        for key, val in kwargs.items():
            string += f',{key}={val}'
        selector = f'//property[{string}]'
        if value:
            selector += '.value'
        return self.select(selector)

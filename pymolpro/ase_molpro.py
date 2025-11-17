import numpy as np
import ase
from ase.calculators.calculator import BaseCalculator, CalculationFailed
from debugpy.common.util import force_str
from networkx.drawing import forceatlas2_layout

from pymolpro import Project
import hashlib
from ase.units import Bohr, Ha
import tempfile
from pathlib import Path


class ASEMolpro(BaseCalculator):
    """ASE calculator interface to Molpro. Implemented using pysjef and pymolpro, which use project bundles to contain input, output and other files, as described at https://github.com/molpro/sjef and https://github.com/molpro/pymolpro.

    """
    implemented_properties = ['energy', 'forces']

    def __init__(self, name='ASEMolpro', method="hf", basis="cc-pVDZ", project_location=None, **kwargs):
        """

        :param name: The file name for the Molpro project
        :type name: str
        :param method: Molpro input fragment that specifies the method
        :type method: str
        :param basis: Basis set specification in the form of arguments to Molpro basis directive
        :type basis: str
        :param project_location: Directory where the project bundle will be placed. Defaults to a unique temporary directory
        :type project_location: str
        :param kwargs: Any additional arguments to pass through to pymolpro.Project()
        :type kwargs: dict
        """
        super().__init__()
        self.molpro_project_name = name
        self.method = method
        self.basis = basis
        self.kwargs = kwargs
        self.projects = {}
        self.project_location = project_location if project_location else tempfile.mkdtemp()

    def calculate(self, atoms, properties, system_changes):
        if True or self.molpro_project_name is None:
            elements = atoms.get_chemical_symbols()
            positions = atoms.get_positions()
            geom_string = str(len(elements)) + '\n'
            for index, c in enumerate(positions):
                geom_string += '\n ' + elements[index]
                for xyz in c:
                    geom_string += ' ' + str(xyz)
            self.molpro_project_name = hashlib.sha256((geom_string +
                                                       str(tuple(sorted(self.kwargs.items())))).encode(
                'utf-8')).hexdigest()[-8:]
        Path(self.project_location).mkdir(parents=True, exist_ok=True)
        forces = 'forces' in properties
        project = Project(self.molpro_project_name, location=self.project_location, geometry=geom_string,
                          method=self.method, basis=self.basis, **self.kwargs,
                          epilogue='{force;varsav}' if forces else None)
        self.projects[self.molpro_project_name] = project
        project.run(wait=True)
        if project.status == 'failed': raise CalculationFailed
        self.results['energy'] = float(project.variable('ENERGY')) * Ha
        if forces:
            gradx = project.variable('GRADX')
            grady = project.variable('GRADY')
            gradz = project.variable('GRADZ')
            self.results['forces'] = -np.array(
                [[float(gradx[i]), float(grady[i]), float(gradz[i])] for i in
                 range(len(gradx))]) * Ha / Bohr

    def clean(self):
        """
        Remove all Molpro projects from file system
        """
        for project in self.projects.values():
            project.erase()

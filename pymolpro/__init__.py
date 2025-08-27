import pysjef
from pysjef import *
from pymolpro.project import no_errors, element_to_dict, resolve_geometry
import pymolpro.node_xml as node_xml
from pymolpro.orbital import Orbital
from pymolpro.tuple import Tuple, Single, Pair
import pymolpro.database
from ._version import __version__, __version_tuple__
from pymolpro.geometry import xyz_to_zmat
from pymolpro.ase_molpro import ASEMolpro
from pymolpro.registry import local_molpro_root, allowed_methods, basis_registry, registry, procedures_registry

__all__ = ['Project', 'Orbital', 'Pair', 'Single', 'Tuple', 'no_errors', 'element_to_dict', 'Database', 'xyz_to_zmat', 'ASEMolpro']

Settings.set('project_default_suffix', 'molpro')
pysjef.project_factory.PROJECT_FACTORY["molpro"] = project.Project

RootXml.TAG_TO_NAME['molpro'] = node_xml.tag_to_name
RootXml.SPECIAL_NODES['molpro'] = {'property': node_xml.PropertyXml,
                                   'plot': node_xml.PlotXml}

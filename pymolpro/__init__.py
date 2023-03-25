import pysjef
from pysjef import *
from pymolpro.project import no_errors, element_to_dict, resolve_geometry
import pymolpro.node_xml as node_xml
from pymolpro.orbital import Orbital
from pymolpro.database import Database, library_database, run_database, compare_database_runs
from pymolpro.tuple import Tuple, Single, Pair
from . import _version

__all__ = ['Project', 'Orbital', 'Pair', 'Single', 'Tuple', 'no_errors', 'element_to_dict', 'Database']

Settings.set('project_default_suffix', 'molpro')
pysjef.project_factory.PROJECT_FACTORY["molpro"] = project.Project

RootXml.TAG_TO_NAME['molpro'] = node_xml.tag_to_name
RootXml.SPECIAL_NODES['molpro'] = {'property': node_xml.PropertyXml,
                                   'plot': node_xml.PlotXml}

__version__ = _version.get_versions()['version']

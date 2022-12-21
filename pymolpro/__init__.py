import pysjef
from pysjef import *
from pymolpro.project import no_errors, element_to_dict
import pymolpro.node_xml as node_xml
from pymolpro.orbital import Orbital
from . import _version

Settings.set('project_default_suffix', 'molpro')
pysjef.project_factory.PROJECT_FACTORY["molpro"] = project.Project

RootXml.TAG_TO_NAME['molpro'] = node_xml.tag_to_name
RootXml.SPECIAL_NODES['molpro'] = {'property': node_xml.PropertyXml,
                                   'plot': node_xml.PlotXml}

__version__ = _version.get_versions()['version']

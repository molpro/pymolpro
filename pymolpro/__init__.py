import pysjef
from pysjef import *
from pymolpro.project import MProject, no_errors
import pymolpro.node_xml as node_xml

Settings.set('project_default_suffix', 'molpro')
pysjef.project_factory.PROJECT_FACTORY["molpro"] = MProject

RootXml.TAG_TO_NAME['molpro'] = node_xml.tag_to_name
RootXml.SPECIAL_NODES['molpro'] = {'property': node_xml.PropertyXml,
                                   'plot': node_xml.PlotXml}

from . import _version
__version__ = _version.get_versions()['version']

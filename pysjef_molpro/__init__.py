import pysjef
from pysjef import *
from pysjef_molpro.project import MProject, no_errors
import pysjef_molpro.node_xml as node_xml

pysjef.PROJECT_FACTORY["molpro"] = MProject
pysjef.PROJECT_FACTORY_DEFAULT_SUFFIX = 'molpro'

pysjef.RootXml.TAG_TO_NAME['molpro'] = node_xml.tag_to_name
RootXml.SPECIAL_NODES['molpro'] = {'property': node_xml.PropertyXml,
                                   'plot': node_xml.PlotXml}
pysjef.RootXml.DEFAULT_SUFFIX = 'molpro'

import pysjef
from pysjef import *
from pymolpro.project import no_errors, element_to_dict
import pymolpro.node_xml as node_xml
from pymolpro.orbital import Orbital
from pymolpro.tuple import Tuple, Single, Pair
from . import _version

Settings.set('project_default_suffix', 'molpro')
pysjef.project_factory.PROJECT_FACTORY["molpro"] = project.Project

RootXml.TAG_TO_NAME['molpro'] = node_xml.tag_to_name
RootXml.SPECIAL_NODES['molpro'] = {'property': node_xml.PropertyXml,
                                   'plot': node_xml.PlotXml}

__version__ = _version.get_versions()['version']


def xpath(element, query):
    """
    Run xpath search on an element in the job xml, with support for default namespace

    :param element: The root element for the search
    :param query: Any xpath search expression supported by lxml.etree.Element
    :return: list of etree.Element objects or of strings, depending on whether the xpath expression results in an attribute
    """
    from lxml import etree
    from io import StringIO
    default_ns_name = '__default__'
    ns = {k if k is not None else default_ns_name: v for k, v in element.getroottree().getroot().nsmap.items()}
    # print("ns:",ns)
    import re
    queryns = re.sub(r'(::|/|^)([_a-zA-Z][-._a-zA-Z0-9]*)(?=/|$|\[)', r'\1' + default_ns_name + r':\2', query)
    # print("queryns",queryns)
    try:
        return element.xpath(queryns, namespaces=ns)
    except Exception as e:
        print("xpath query failed:", e)
        print("query =", query, ", queryns =", queryns, ", namespaces =", ns)

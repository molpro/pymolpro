from pysjef import RootXml
import regex as re


def tag_to_name(tag):
    """
    Takes XML tag with namespace, outputs raw node name.
    """
    # TODO could we use etree.QName for this?
    _NAMESPACE_RE = re.compile(r"{(.+)}(.*)")
    match = _NAMESPACE_RE.match(tag)
    tag = tag if match is None else match.group(2)
    namespace = None if match is None else match.group(1)
    if namespace == "http://www.xml-cml.org/schema":
        tag = "cml" + tag.capitalize()
    return tag


class PropertyXml(RootXml):
    """
    molpro-output:property - Output about a single property of a single molecule
    Attributes:
        str name, any value
    """


class CmlMoleculeXml(RootXml):
    """cml:molecule - Output about a single molecule
    Children:
        cml:atomarry, cml:symmetry
    """

    # cmlMolecule has 2 jobs:
    # 1. steal information from nodes within it (e.g. cmlsymmetry)
    # 2. allow its own information to be stolen by <molecule>
    # def finalise(self, **options):
    #     # find my child named CmlSymmetry and take point group
    #     child = [ch for ch in self.children if ch.nodename == 'cmlSymmetry'][0]
    #     self.attributes.update(child.attributes)
    #     self.children.remove(child)
    #     # take coordinates from my atoms and make a geometry array
    #     # regenerate my own list of properties
    #     # TODO eventually take attributes from parent <molecule> as well


class CmlAtomXml(RootXml):
    """cml:atom - Output about a single atom
    """

    # def finalise(self, **options):
    #     # take properties from myself and give them more intuitive names
    #     # I have: id, elementType, x3, y3, z3
    #     # x3, y3 and z3 are all bad
    #     # rename to x, y, z AND provide a coords tuple
    #     x = self._rename_attribute('x3', 'x')
    #     y = self._rename_attribute('y3', 'y')
    #     z = self._rename_attribute('z3', 'z')
    #     # append coords
    #     self.attributes['coords'] = (x, y, z)


class InputXml(RootXml):
    """molpro-output:input - Molpro input dataset
    Children:
        p, include
    """

    # def finalise(self, **options):
    #     # <input> needs to take its <p> children, combine them into a single
    #     # list, and then put that list onto <molpro> as an attribute
    #     # should be <molpro><job><input>
    #     assert self.parent.parent.nodename == "molpro"
    #     # at finalise, my <p> children are initialised
    #     inp = [child.attributes['.'] for child in self.children]
    #     # say farewell to the kids
    #     self.children = []
    #     # add inp as one of <molpro>'s attributes
    #     # self.parent.parent.attributes['input'] = inp
    #     # <molpro> will regenerate its own properties after this
    #     # now delete myself
    #     # self.parent.children = [child for child in self.parent.children
    #     # if child is not self]
    #     self.attributes['value'] = inp
    #     # my parent does not know who I am and I will soon be garbage collected
    #     # XXX change this if XML regeneration is a desired feature


class JobstepXml(RootXml):
    """molpro-output:jobstep - Output from a single job step
    Children:
        error, gradient, opt, property, cube, time, storage, vibrations,
        jobstep, cml:molecule, variables, instanton
    Attributes:
        str command, str commandset, bool displaced
    """


class MoleculeXml(RootXml):
    """molpro-output:molecule - Summary output about the job
    Children:
        cml:molecule, basisSet, vibrations, orbitals, stm:metadataList,
        platform
    Attributes:
        str id, str index, str InChI, str InChIKey, str SMILES, str title, str
        method, str basis, str geometryMethod, str geometryBasis, double energy
    """
    # <molecule> is very poorly named - it contains output about the final
    # result of the job, but the actual molecular geometry is still contained
    # within <cml:molecule>
    # therefore the filter for <molecule> will be named summary, and the
    # <cml:molecule> that is contained within <molecule> will not be accessible
    # by the molecules filter from outside the summary
    # however, <molecule> shouldn't actually have to do anything


class VariableXml(RootXml):
    """
    TODO Doc
    """
    # def finalise(self, **options):
    #     # steal value as attribute
    #     values = []
    #     for child in self.children:
    #         try:
    #             values.append(child.attributes['.'])
    #         except KeyError:
    #             values.append(None)
    #     if len(values) == 1:
    #         value = to_numerical(values[0])
    #     else:
    #         # join the numbers so that _to_float will parse them
    #         value = np.array([to_numerical(v) for v in values])
    #         # returns np.array
    #     self.attributes['value'] = value
    #     self.children = []


class VariablesXml(RootXml):
    """molpro-output:variables - Container for Molpro internal variables
    Children:
        variable
    """


class OrbitalsXml(RootXml):
    """molpro-output:orbitals - Container for Molpro internal orbitals
    Children:
        orbital
    Attributes:
        str basis, str angular, str spin, str method
    """

    # The <orbitals> node contains a bunch of <orbital>s.
    # Rather than accessing each orbital individually, it would be preferred to
    # compile a list of energies/etc that are placed as an attribute onto the
    # <orbitals> node.
    # def finalise(self, **options):
    #     # <orbital> content is a list of floats - these should be copied and
    #     # taken as an attribute on this node
    #     # convert orbital data to 2d arrays
    #     data = []
    #     for child in self.children:
    #         try:
    #             data.append(to_numerical(child.attributes['.']))
    #         except KeyError:
    #             data.append(None)
    #     symmetry_ids = [child.attributes['symmetryID'] for child in
    #                     self.children]
    #     occupations = [child.attributes['occupation'] for child in
    #                    self.children]
    #     energies = [child.attributes['energy'] for child in self.children]
    #     self.attributes['data'] = np.array(data)
    #     self.attributes['energies'] = np.array(energies)
    #     self.attributes['symmetryIDs'] = np.array(symmetry_ids)
    #     self.attributes['occupations'] = np.array(occupations)


class OrbitalXml(RootXml):
    """molpro-output:orbital - A Molpro molecular orbital
    Children:
        [string]
    Attributes:
        str symmetryID, double energy, double occupation
    """


# class PXml(MolproXml):
#     """molpro-output:p - Contains text.
#     Attributes:
#         int depth
#     """

class PlatformXml(RootXml):
    """molpro-output:platform - Container for metadata about the platform
    Children:
        version, licence, parallel, machine, dimensions
    """


class VersionXml(RootXml):
    """molpro-output:version
    Children:
        [string], date
    Attributes:
        int major, int minor, str SHA, int integer_bits, str parallelism
    """


class DateXml(RootXml):
    """molpro-output:date - The time of the calculation.
    Children:
        [string datetime]
    Attributes:
        int year, int month, int day, int hour, int minute, int second
    """


class LicenceXml(RootXml):
    """molpro-output:licence
    Attributes:
        str id
    """


class ParallelXml(RootXml):
    """molpro-output:parallel
    Attributes:
        int processes, int nodes, int all_processes, int openmp
    """


class MachineXml(RootXml):
    """molpro-output:machine
    Attributes:
        str hostname
    """


# This class is already defined as part of Cube
# class DimensionsXml(MolproXml):
#     """molpro-output:dimensions
#     Attributes:
#         [undocumented]
#     """

class DiagnosticsXml(RootXml):
    """molpro-output:diagnostics - Diagnostic info from the job
    Attributes:
        int warnings
    """


class BasisSetXml(RootXml):
    """molpro-output:basisSet - Container for a basis set
    Children:
        basisGroup, association
    Attributes:
        str id, str type, str angular, int groups, str primitives, int length,
        int cartesianLength
    """


class VibrationsXml(RootXml):
    """molpro-output:vibrations - Container for info on molecular vibrations
    Children:
        normalCoordinate
    Attributes:
        str name, str type, str units, int length
    """


class NormalCoordinateXml(RootXml):
    """molpro-output:normalCoordinate - A normal vibrational coordinate
    Children:
        [list of doubles]
    Attributes:
        double wavenumber, str units, double IRintensity, str IRintensityunits,
        str symmetry, str real_zero_imag
    """


class CubeXml(RootXml):
    """molpro-output:cube - Metadata for a cube dataset
    Children:
        field, dimensions, origin, axes, step
    Attributes:
        str method
    """


class FieldXml(RootXml):
    """molpro-output:field
    Attributes:
        str quantity, str type, str number, str symmetry, double occupancy,
        double energy, str file
    """


class DimensionsXml(RootXml):
    """molpro-output:dimensions
    Children:
        [list of doubles]
    """


class OriginXml(RootXml):
    """molpro-output:origin
    Children:
        [list of doubles]
    """


class AxesXml(RootXml):
    """molpro-output:axes
    Children:
        [list of doubles]
    """


class StepXml(RootXml):
    """molpro-output:step
    Children:
        [list of doubles[
    """


class PlotXml(RootXml):
    """molpro-output:plot - Metadata for plot from a Molpro table
    Attributes:
        str table, str plot, str type
        Also attributes from http://www.w3.org/1999/xlink
    """



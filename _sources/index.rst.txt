.. pymolpro documentation master file, created by
   sphinx-quickstart on Thu Oct 13 06:35:46 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pymolpro (version |release|)
============================

**pymolpro** is a Python library that provides support
for working with the `Molpro quantum chemistry package <https://www.molpro.net/>`_.

The principal feature is
the :py:meth:`Project` class that provides access to a complete Molpro job, including input
and output files together with metadata such as job status information.
The project is stored as a bundle implemented as a directory in the file system.
The class is a Python binding of the
`sjef <https://molpro.github.io/sjef/>`_ library with additional Molpro-specific customisation,
and the project bundle can be accessed
also through a command-line interface and the `gmolpro <https://www.molpro.net/manual/doku.php?id=gmolpro_graphical_user_interface>`_
graphical user interface.

For technical reasons, :py:meth:`Project` invokes the submodular class :py:meth:`pymolpro.project.Project`.
Normally all that is required to instantiate a class, on either a new or existing bundle,
is to pass a single argument which is the path of the bundle, with or without the compulsory
suffix `.molpro`.
Simply running Molpro and inspecting the output can be achieved as::

  from pymolpro import Project
  p=Project("Neon")
  p.write_input("""
  geometry={Ne}
  rhf
  ccsd(t)
  """)
  p.run()
  print(p.out)

Normally, :py:meth:`pymolpro.project.Project.run()` launches the job in the background, potentially on a remote machine via a suitably configured `backend <https://molpro.github.io/sjef/md____w_sjef_sjef_src_sjef_backends.html>`_,
and its status is available through the :py:attr:`pymolpro.project.Project.status` attribute.
Additionally, the job launching request is ignored if the project has already run the job successfully using the same input; this introduces the convenience of being able to re-run an entire workflow, for example in a Jupyter notebook, without needless recomputation.

Molpro produces an xml file that contains all the essential results, marked up using the `molpro-output <https://www.molpro.net/schema/molpro-output>`_ schema.
The pymolpro Project class contains functions that interpret that output, and at the lowest level this can be
achieved through a suitable `XPath <https://www.w3.org/TR/1999/REC-xpath-19991116/>`_ search.
For extracting calculated properties, including energies, convenience functions
:py:meth:`pymolpro.project.Project.properties()`,
:py:meth:`pymolpro.project.Project.energies()`,
that wrap
suitable XPath searches, are provided.
These functions can return either a Python dictionary containing all information about the property, or just its value.
The following is a simple example where a job is created, run and
analysed.::

   from pymolpro import Project
   p = Project("Neon")
   p.write_input("geometry={Ne}; rhf; ccsd(t)")
   p.run(wait=True)
   assert p.status == 'completed' and not p.errors()
   energies = {}
   for node in p.xpath("//property[@name='Energy' or @name='total energy']"):
     energies[node.xpath("@method")[0]] = float(node.xpath("@value")[0])
   energy_values = p.energies()
   all_principal_properties_dict = p.properties(principal=True, value=False)
   final_energy = p.energy(method='CCSD(T)')

Installation
------------
pymolpro is on `conda forge <https://conda-forge.org>`_ and can be installed on most systems using ``conda install -c conda-forge pymolpro``.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   modules
   examples

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

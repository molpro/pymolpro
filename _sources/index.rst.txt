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
suffix `.molpro`. The following is a simple example where a job is created, run and
analysed.::

   from pymolpro import Project
   p = Project("Neon")
   p.write_input("geometry={Ne}; rhf; rks,b3lyp")
   p.run(wait=True)
   print(p.xpath_search("//property[@name='Energy']", "value"))
   print(p.properties('Energy', value=True))



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

Molpro python support
=====================

**pymolpro** is a Python library that provides support
for working with the `Molpro quantum chemistry package <https://www.molpro.net/>`_.

The principal feature is
the `Project` class that provides access to a complete Molpro job, including input
and output files together with metadata such as job status information.
The project is stored as a bundle implemented as a directory in the file system.
The class is a Python binding of the
`sjef <https://molpro.github.io/sjef/>`_ library with additional Molpro-specific customisation,
and the project bundle can be accessed
also through a command-line interface and the `gmolpro <https://www.molpro.net/manual/doku.php?id=gmolpro_graphical_user_interface>`_
graphical user interface.

pymolpro is on `conda forge <https://conda-forge.org>`_ and can be installed on most systems using ``conda install -c conda-forge pymolpro``.

Documentation at https://molpro.github.io/pymolpro.

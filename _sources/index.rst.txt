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
   all_principal_properties_dict = p.properties(principal=True, dict=True)
   final_energy = p.energy(method='CCSD(T)')

Installation
------------
pymolpro is on `conda forge <https://conda-forge.org>`_ and can be installed on most systems using ``conda install -c conda-forge pymolpro``.
On Microsoft Windows, you need to also install `msys2 <https://www.msys2.org>`_,
and then, in an msys command window, ``pacman -S rsync openssh``.

For a complete set-up of pymolpro within Jupyter notebooks driven from the comand line on linux or macOS,

* Install `Molpro <https://www.molpro.net>`_
* Install `Miniconda <https://www.anaconda.com>`_. On macOS with homebrew, this can be done with ``brew install miniconda``. You may need to then type ``conda init`` and restart the shell.
* Create a conda environment to contain at least these packages; of course you can add any others you need, e.g. ``matplotlib``: ``conda create -n pymolpro-jupyter -c conda-forge pymolpro jupyter nb_conda_kernels``.
* If you intend to run Molpro jobs on a remote machine, set up password-free ssh access to it, using the following or otherwise:  ``ssh-keygen; ssh-copy-id`` user@host. Then edit the file ``~/.sjef/molpro/backends.xml``  to `specify the backend <https://molpro.github.io/sjef/md____w_sjef_sjef_src_sjef_backends.html>`_.
*  Activate the environment: ``conda activate pymolpro-jupyter``. To create a new notebook, ``jupyter notebook``; to open an existing notebook ``jupyter notebook existing_file.ipynb``

For a complete set-up of pymolpro within Jupyter notebooks on Windows,

* Install `Molpro <https://www.molpro.net>`_, and msys2 as described above.
* Install `Anaconda3 <https://www.anaconda.com>`_
* Open Anaconda Navigator
* Select the ``Environments`` tab and press the ``Create`` button. Give your new conda environment a name, e.g. ``pymolpro-jupyter``.
* With the ``Channels`` button, add the ``conda-forge`` channel, and then select ``Not installed`` and search for, and install ``pymolpro``, ``jupyter``, ``nb_conda_kernels``, as well as any other packages you might want, such as ``matplotlib``.
* If you intend to run Molpro jobs on a remote machine, from the command menu run ``MSYS2``, and (a) create a password-free ssh key by running ``ssh-keygen``; (b) authorize that key for login on the desired remote machine: ``ssh-copy-id`` user@host; (c) edit the file ``.sjef\molpro\backends.xml`` (relative to your Windows home directory) to `specify the backend <https://molpro.github.io/sjef/md____w_sjef_sjef_src_sjef_backends.html>`_.
* Select the ``Home`` tab and run the Jupyter Notebook application, ensuring that the selected environment is ``pymolpro-jupyter``.

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

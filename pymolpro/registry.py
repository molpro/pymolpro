import atexit
import shutil
import tempfile
import pathlib

from pymolpro import Project

_registry_project = None
_registry_project_directory = None


def ensure_registry_project():
    global _registry_project
    global _registry_project_directory
    if _registry_project is None:
        _registry_project_directory = tempfile.mkdtemp()
        _registry_project = Project((pathlib.Path(_registry_project_directory) / 'registry_project').as_posix())
        atexit.register(_registry_project_delete)

def _registry_project_delete():
    _registry_project.erase()
    shutil.rmtree(_registry_project_directory, ignore_errors=True)

def local_molpro_root():
    r"""
    Get the directory of the Molpro installation in the local backend

    :return: directory
    :rtype: pathlib.Path
    """
    ensure_registry_project()
    return _registry_project.local_molpro_root


def procedures_registry():
    r"""
    Get the procedures registry from the Molpro pointed to in the local backend

    :return:
    :rtype: dict
    """
    ensure_registry_project()
    return _registry_project.procedures_registry()


def basis_registry():
    r"""
    Get the basis registry from the Molpro pointed to in the local backend

    :return:
    :rtype: dict
    """
    ensure_registry_project()
    return _registry_project.basis_registry()


def registry(set=None):
    """
    Access the registry belonging to Molpro in the `local` backend.

    :param set: If present, obtain as a dictionary the specified registry set. If absent, obtain a list of all sets
    :return:
    """
    ensure_registry_project()
    return _registry_project.registry(set)


def allowed_methods():
    result = []
    _procedures_registry = procedures_registry()
    for keyfound in _procedures_registry.keys():
        if _procedures_registry[keyfound]['class'] == 'PROG':
            result.append(_procedures_registry[keyfound]['name'])
    return result

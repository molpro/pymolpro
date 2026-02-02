import pathlib
import logging
logger = logging.getLogger(__name__)

_registry_project = None


def _ensure_registry_project():
    from .project import Project
    global _registry_project
    if _registry_project is None:
        _registry_project = Project()
        logger.debug('pymolpro._ensure_registry_project() creates')


def local_molpro_root() -> pathlib.Path:
    r"""
    Get the directory of the Molpro installation in the local backend

    :return: directory
    """
    _ensure_registry_project()
    logger.debug('pymolpro.local_molpro_root() -> ' + str(_registry_project.local_molpro_root))
    return _registry_project.local_molpro_root


def procedures_registry():
    r"""
    Get the procedures registry from the Molpro pointed to in the local backend

    :return:
    :rtype: dict
    """
    _ensure_registry_project()
    logger.debug('pymolpro.procedures_registry() -> ' + str(_registry_project.procedures_registry()))
    return _registry_project.procedures_registry()


def basis_registry():
    r"""
    Get the basis registry from the Molpro pointed to in the local backend

    :return:
    :rtype: dict
    """
    _ensure_registry_project()
    return _registry_project.basis_registry()


def registry(set=None):
    """
    Access the registry belonging to Molpro in the `local` backend.

    :param set: If present, obtain as a dictionary the specified registry set. If absent, obtain a list of all sets
    :return:
    """
    _ensure_registry_project()
    return _registry_project.registry(set)


def allowed_methods():
    result = []
    _procedures_registry = procedures_registry()
    for keyfound in _procedures_registry.keys():
        if _procedures_registry[keyfound]['class'] == 'PROG':
            result.append(_procedures_registry[keyfound]['name'])
    return result

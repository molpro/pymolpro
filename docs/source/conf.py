# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'pymolpro'
copyright = '2022, Marat Sibaev, Peter Knowles'
author = 'Marat Sibaev, Peter Knowles'

import sys
from pathlib import Path
# sys.path.insert(0, Path(__file__).parent.parent.absolute())
import subprocess

subprocess.check_call(
    [sys.executable, "-m", "pip", "install", "--no-deps", "--force-reinstall", Path(__file__).parent.parent.parent])
import pymolpro

release = pymolpro.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'nbsphinx',
]

templates_path = ['_templates']
exclude_patterns = []

variables_to_export = ["release"]
frozen_locals = dict(locals())
rst_epilog = '\n'.join(map(lambda x: f".. |{x}| replace:: {frozen_locals[x]}", variables_to_export))
del frozen_locals

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']

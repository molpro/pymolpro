import sys
import subprocess
import pymolpro
import os
from pathlib import Path

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'pymolpro'
copyright = '2022, Marat Sibaev, Peter Knowles'
author = 'Marat Sibaev, Peter Knowles'

# sys.path.insert(0, Path(__file__).parent.parent.absolute())

subprocess.check_call(
    [sys.executable, "-m", "pip", "install", "--no-deps", "--force-reinstall", Path(__file__).parent.parent.parent])

release = pymolpro.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'nbsphinx',
    'sphinx-jsonschema',
]

templates_path = ['_templates']
exclude_patterns = []

variables_to_export = ["release"]
__frozen_locals = dict(locals())
rst_epilog = '\n'.join(map(lambda x: f".. |{x}| replace:: {__frozen_locals[x]}", variables_to_export))

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
# html_static_path = ['_static']


for dbname in pymolpro.database.library():
    dbname_pretty = dbname.replace('_', ' ')
    db = pymolpro.database.load(dbname)
    if not os.path.exists('database'):
        os.makedirs('database')
    with open('database/' + dbname + '.rst', 'w') as f:
        f.write(db.__str__(rst=True, geometry=False, title=dbname_pretty))

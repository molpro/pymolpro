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
sys.path.insert(0, Path(__file__).parent.parent.absolute())
# import pymolpro
import pip
# pip.main(['install']+[Path(__file__).parent.parent.as_posix(),"--no-deps","--force-reinstall"])
# import os.path
# BASE_DIR = os.path.abspath(os.path.join(__file__, '../../../'))
# BASE_DIR=Path(__file__).parent.parent.parent
# BASE_DIR='pymolpro'
# pip.main(['install',"--no-deps","--force-reinstall",BASE_DIR])
# pip.main(['install'])
import subprocess
subprocess.check_call([sys.executable,"-m","pip","install","--no-deps","--force-reinstall",Path(__file__).parent.parent.parent])
import pymolpro
release = pymolpro.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

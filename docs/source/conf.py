# Configuration file for the Sphinx documentation builder.

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))  # <-- Esto es esencial, añade la raíz del proyecto (python/)

# -- Project information -----------------------------------------------------
project = 'DFToy'
copyright = '2025, Juan Esteban Neira Díaz, Julián Cogua Arévalo'
author = 'Juan Esteban Neira Díaz, Julián Cogua Arévalo'
release = '1.0'

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary"
]
autosummary_generate = True

templates_path = ['_templates']
exclude_patterns = []

language = 'es'

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']


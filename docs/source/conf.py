# Configuration file for the Sphinx documentation builder.

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))  # raíz del proyecto, donde están scf/ y gui_app/

# -- Project information -----------------------------------------------------
project = 'dftoy'
author = 'Julian Cogua, Juan Esteban Neira'
release = '1.2'

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",      # Extrae docstrings automáticamente
    "sphinx.ext.napoleon",     # Soporte para docstrings estilo NumPy/Google
    "sphinx.ext.viewcode"      # Enlaces al código fuente
]

templates_path = ['_templates']
exclude_patterns = []

language = 'es'

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'   # Tema moderno tipo ReadTheDocs
html_static_path = ['_static']

# Configuration file for the Sphinx documentation builder.
# See https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Agrega la raíz del proyecto al path (donde están scf/ y gui_app/)
sys.path.insert(0, os.path.abspath('../..'))  # sube dos niveles desde docs/source

# -- Project information -----------------------------------------------------
project = 'DFToy'
copyright = '2025, Julian Cogua, Juan Esteban Neira'
author = 'Julian Cogua, Juan Esteban Neira'
release = '1.2'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',     # Extrae docstrings automáticamente
    'sphinx.ext.napoleon',    # Soporte para docstrings estilo NumPy
    'sphinx.ext.viewcode',    # Enlaces al código fuente
    'sphinx.ext.todo',        # Permite directivas TODO
]

templates_path = ['_templates']
exclude_patterns = []

language = 'es'

# -- HTML output -------------------------------------------------------------
html_theme = 'sphinx_rtd_theme'  # Tema moderno similar a ReadTheDocs
html_static_path = ['_static']

# -- Autodoc options ---------------------------------------------------------
autodoc_member_order = 'bysource'
autodoc_typehints = 'description'
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True
}

# -- Napoleon settings -------------------------------------------------------
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_param = True
napoleon_use_rtype = True

# -- Todo extension ----------------------------------------------------------
todo_include_todos = True


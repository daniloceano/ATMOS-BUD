# Configuration file for the Sphinx documentation builder.

import os
import sys

# Add the project's source directory to sys.path
sys.path.insert(0, os.path.abspath('../src'))
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------
project = 'ATMOS-BUD'
copyright = '2025, Danilo Couto de Souza'
author = 'Danilo Couto de Souza'
release = '0.1.6'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
]

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False

# Autodoc settings
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'undoc-members': True,
    'exclude-members': '__weakref__'
}

# Mock imports for ReadTheDocs
autodoc_mock_imports = [
    'numpy', 'pandas', 'xarray', 'matplotlib', 'cartopy', 
    'cmocean', 'cdsapi', 'scipy', 'netCDF4', 'shapely', 
    'metpy', 'dask', 'sklearn'
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# Theme options
html_theme_options = {
    'canonical_url': '',
    'analytics_id': '',
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}

html_title = f"ATMOS-BUD v{release} Documentation"
html_short_title = "ATMOS-BUD Docs"

# Intersphinx mapping
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
    'xarray': ('https://docs.xarray.dev/en/stable/', None),
}

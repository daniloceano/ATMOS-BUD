# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Mock imports to prevent actual code execution during documentation build
class MockModule:
    """Mock module to prevent imports from executing code"""
    def __getattr__(self, name):
        return MockModule()
    
    def __call__(self, *args, **kwargs):
        return MockModule()
    
    def __getitem__(self, key):
        return MockModule()
    
    def __setitem__(self, key, value):
        pass

# Mock scientific computing packages that might execute code on import
MOCK_MODULES = [
    'numpy', 'pandas', 'xarray', 'matplotlib', 'matplotlib.pyplot',
    'cartopy', 'cartopy.crs', 'cartopy.feature', 'cmocean', 'cmocean.cm',
    'cdsapi', 'sklearn', 'sklearn.preprocessing', 'scipy', 'netCDF4',
    'shapely', 'shapely.geometry', 'metpy', 'metpy.calc', 'metpy.units',
    'dask', 'dask.array'
]

for mod_name in MOCK_MODULES:
    sys.modules[mod_name] = MockModule()

# Add the src directory to the Python path
sys.path.insert(0, os.path.abspath('../src'))
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'ATMOS-BUD'
copyright = '2024, Danilo Couto de Souza'
author = 'Danilo Couto de Souza'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinx.ext.autosummary'
]

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True

# Autodoc settings
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
    'exclude-members': '__weakref__'
}

# Prevent autodoc from executing code
autodoc_mock_imports = MOCK_MODULES
autodoc_preserve_defaults = True

# Additional safety settings
autoclass_content = 'class'  # Only include class docstring, not __init__
autodoc_typehints = 'description'

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# Theme options for ReadTheDocs
html_theme_options = {
    'canonical_url': '',
    'analytics_id': '',
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'vcs_pageview_mode': '',
    'style_nav_header_background': '#2980B9',
    # Toc options
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}

# Additional HTML options
html_title = f"ATMOS-BUD v{release} Documentation"
html_short_title = "ATMOS-BUD Docs"
html_logo = None
html_favicon = None

# Links and metadata
html_context = {
    "display_github": True,
    "github_user": "daniloceano",
    "github_repo": "ATMOS-BUD",
    "github_version": "main",
    "conf_py_path": "/docs/",
}

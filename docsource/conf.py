import os
import sys
import sphinx_rtd_theme

# Add the path to your package to the system path
sys.path.insert(0, os.path.abspath('../'))

# Project information
project = 'foldtree'
author = 'Dave Moi'

# Extensions to use
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.coverage',
    'sphinx.ext.napoleon',
    'sphinx_rtd_theme'
]

# Theme settings
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 3,
    'style_external_links': True
}

# Add any additional options for autodoc
autodoc_default_options = {
    'member-order': 'bysource'
}

# Add any modules to be excluded from the documentation
exclude_patterns = []

# The master toctree document.
master_doc = 'index'

# Mock import modules that may not be available in the documentation build environment
autodoc_mock_imports = ['requests' , 'Bio' , 'numpy' , 'pandas' , 'matplotlib' , 'seaborn' , 'scipy' , 'wget' , 'statsmodels' , 'toytree' , 'pandas' , '' ]
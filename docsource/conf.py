import os
import sys

sys.path.insert(0, os.path.abspath('../../'))
sys.path.insert(0, os.path.abspath('../src/'))

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'sphinx.ext.viewcode']
html_theme = 'sphinx_rtd_theme'

import AFDB_tools, foldseek2tree, treescore , corecut 


# List of modules to mock import when building documentation
#autodoc_mock_imports = ['treescore' , 'foldseek2tree' , 'AFDB_tools']


# Default order for members
autodoc_member_order = 'bysource'
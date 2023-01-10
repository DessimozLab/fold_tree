import os
import sys

sys.path.insert(0, os.path.abspath('../../'))
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon']
html_theme = "classic"


# List of modules to mock import when building documentation
autodoc_mock_imports = ['treescore' , 'foldseek2tree' , 'AFDB_tools']

# Default order for members
autodoc_member_order = 'bysource'
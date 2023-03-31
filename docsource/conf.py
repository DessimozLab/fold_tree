import os
import sys

sys.path.append(os.path.abspath('../src/'))

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'sphinx.ext.viewcode']
html_theme = 'sphinx_rtd_theme'



autodoc_mock_imports = ['scipy' , 'numpy' , 'wget' , 'statsmodels' , 'toytree', 'pandas' , 're' , 'os' , 'scipy.stats' , 'argparse' , 'subprocess' , 'shlex' , 
'scipy.spatial.distance' , 'scipy.stats' , 'statsmodels', 'Bio' , 'requests'
]

import AFDB_tools, foldseek2tree, treescore , corecut 


# List of modules to mock import when building documentation


# Default order for members
autodoc_member_order = 'bysource'

.. toctree::
   :maxdepth: 2
   
   treeinspector
   foldseek2tree
   AFDBtools
   treescore


Welcome to My Project
=====================

This is the documentation for foldtree, it's a combination of some utility functions and a snakemake workflow to make trees from alphafold structures.

.. image:: https://img.shields.io/pypi/v/foldtree.svg
   :target: https://pypi.python.org/pypi/foldtree





Installation
------------

To use foldtree, clone the git rep and install it with pip:

test project
.. code-block:: bash

   $ git clone 
   $ cd foldtree
   $ pip install -e .


Usage
-----

You can use the functions in the different source module to design your own workflows or use the snakemake workflow.

To run the snakemake workflow on the test dataset try using 

.. code-block:: bash

   $ snakemake --cores 4 --use-conda --configfile config.yaml --config input_dir=test_data/ output_dir=test_output/ --use-conda --conda-prefix /tmp/conda


Or if you are using a slurm cluster you can use the slurm profile:

.. code-block:: bash

   $ snakemake --profile slurmsimple --use-conda --configfile config.yaml --config input_dir=test_data/ output_dir=test_output/ --use-conda --conda-prefix /tmp/conda

This workflow will try to create a tree for each of the uniprot identifiers in the identifier.txt file in the test_data folder. The output will be in the test_output folder.
You can use this workflow on your own data by creating a folder with an identifier file. The workflow will attempt to to download all the necesary data to produce a sequence and structure based tree.


To use the foldtree utility functions, import the different modules in your code:

.. code-block:: python

   from foldtree.src import foldseek2tree
   from foldtree.src import AFDBtools
   from foldtree.src import treescore

   
   myproject.do_x()






Troubleshooting
---------------


If you encounter any issues while using My Project, please file a bug report on our GitHub repository: https://github.com/user/repo/issues


Credits
-------

My Project was created by John Doe.
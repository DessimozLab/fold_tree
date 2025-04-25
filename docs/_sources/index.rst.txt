
.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   treeinspector
   foldseek2tree
   AFDBtool
   treescore



Foldtree
=====================

This is the documentation for foldtree, it's a combination of some utility functions and a snakemake workflow to make trees from alphafold structures.

Installation
------------

To use foldtree, first install snakemake:
https://snakemake.readthedocs.io/en/stable/getting_started/installation.html


Now we can clone the repo and create a conda environment with the software need to run fold tree


.. code-block:: bash
   $git clone git@github.com:DessimozLab/fold_tree.git
   $cd cd fold_tree
   $mamba create -n foldtree --file= ./workflow/config/fold_tree.yaml
   $mamba activate foldtree


Now, we're ready to run the pipeline. For most users, using the fold_tree pipeline should be sufficient for their needs. 
You can setup a fold_tree run by creating a folder for your output. 


.. code-block:: bash
   $mkdir myfam

Now we can either add an identifiers.txt file containing the uniprot identifiers of all of the proteins we would like to make a tree with.

.. code-block:: bash
   $mkdir myfam


.. code-block:: bash
   └── myfam
      └── identifiers.txt




Or we can run the pipeline on our own set of structures. Please note that discontinuities or other defects in the PDBs may adversly affect results.
Let's make our structure directory and add some PDB files to it. In this case the identifiers file is blank.

.. code-block:: bash
   └── myfam
      ├── identifiers.txt
      └── structs
         ├── struct1.pdb
         ├── struct2.pdb
         └── struct3.pdb

Now we're ready to build our trees. Let's run the pipeline.



Usage
-----

To run the snakemake workflow on the test dataset try using. You can change the folder variable to the location of your data. 

.. code-block:: bash

   $ snakemake --cores 4 --use-conda -s ./workflow/fold_tree --config folder=./testdata filter=False customstructs=False  --use-conda 



Or if you are using a slurm cluster you can use the slurm profile:

.. code-block:: bash

   $ snakemake --cores 4 --use-conda -s ./workflow/fold_tree --config folder=./testdata filter=False customstructs=False  --profile slurmsimple --use-conda 

The fold_tree workflow will create a tree for each of the uniprot identifiers in the identifier.txt file in the input folder. 

To use custom structures leave a blank identifier file and set the customstructs variable to True. 

.. code-block:: bash

   $ snakemake --cores 4 --use-conda -s ./workflow/fold_tree --config folder=./myfam filter=False customstructs=True  --profile slurmsimple --use-conda 


To use the foldtree utility functions in your own work first install the repo as a python library.

.. code-block:: bash

   $ git clone 
   $ cd foldtree
   $ pip install -e .

Then import the libraries in your script or notebooks

.. code-block:: python

   from foldtree.src import foldseek2tree
   from foldtree.src import AFDBtools
   from foldtree.src import treescore

   
   # use the functions somehow.
   # comments/help are provided in the code

There are also examples of how to use the different functions in the notebooks in the notebooks folder.

Troubleshooting
---------------

If you encounter any issues while using My Project, please file a bug report on our GitHub repository: https://github.com/DessimozLab/fold_tree/issues


Credits
-------

This project was created by Dave Moi, Yannis Nevers and Charles Bernard at DessimozLab (DBC at the university of Lausanne).
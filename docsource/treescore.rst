docstrings
---------------

The tree string module has a few functions that are useful for
scoring trees based on the taxonomic information in the tree. 
Trees that are more taxonomically plausible will have higher scores.

The module uses the taxonomic information from uniprot to score trees
using a recursive algorithm.  The algorithm is described in the
following paper:

    Simple chained guide trees give poorer multiple sequence alignments than inferred trees in simulation and phylogenetic benchmarks
    Tan G, Gil M, Löytynoja AP, Goldman N, Dessimoz C
    Proc. Natl. Acad. Sci. U. S. A., 2015

    https://pubmed.ncbi.nlm.nih.gov/25564672/

.. automodule:: src.treescore
    :members:    
    :undoc-members:
    :show-inheritance:


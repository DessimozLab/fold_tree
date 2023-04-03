
AFDB tools
=====================

The AFDB tools module contains utility functions to download target AlphaFold structural models, perform actions on them and interact with the UniProt server.

Briefly, the user can:

- Download a set of target AlphaFold structural models using the `grab_struct` function.
- Get descriptive statistics for an AlphaFold model's pLDDT metrics (minimum and maximum value, mean, variance, etc) using the `descr` function.
- Filter AlphaFold structural models based on their residue's pLDDT metric using the `filter_plddt` function, discarding models presenting mean and/or minimum pLDDT values below a specified score.
- Request protein information from the UniProt server for a given protein name or set of proteins using the `unirequest` and `grab_entries` functions. 
- Convert data frames obtained from the UniProt server into FASTA format using the `res2fasta` function.

.. automodule:: src.AFDB_tools
    :members:
    :undoc-members:
    :show-inheritance:

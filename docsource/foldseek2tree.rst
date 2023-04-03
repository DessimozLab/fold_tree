
foldseek2tree
=====================

The `foldseek2tree` module has the main functions needed to infer a structural phylogenetic tree given a target set of structural models in PDB format. 

Briefly, a structural phylogenetic tree can be infered under this framework performing the following steps: 
- Create a Foldseek database for the target set of structural models using the `runFoldseekdb` function.
- Run pairwise all-vs-all structural comparisons with Foldseek using the `runFoldseek_allvall` function (use `runFoldseek_allvall_EZsearch` to run Foldseek under the easy-search scheme of work).
- Create a distance matrix suitable for phylogenetic inference using the `distmat_to_txt`, using Foldseek comparison results as input data.
- Employ the above mentiohed distance matrix to infer a phylogenetic tree employing distance-based methods of inference such as FastME (`runFastme`).

This task is automated using `foldtree`'s main pipeline of work, but the user can use the `structblob2tree` function as a wrapper function to perform them without running the whole pipeline. Please note that under `foldtree`'s pipeline distance matrices are computed as (1-x), where x is the corresponding similarity metric (_e.g._ Foldseek's 3Di/aa score) and infered trees are post-processed for subsequent analyses, making all branches positive with the `postprocess` function.

The `foldseek2tree` module also comprises a set of tools to perform minor accompanying analyses, such as:
- Testing if a distance matrix behaves in a clock-like matter using a $\chi$^2 distribution, based on Tajima's test (`clock_test` function; Tajima, 1993).
- Get a consensus tree for a list of input trees (`consensustree` function).

.. automodule:: src.foldseek2tree
    :members:
    :undoc-members:
    :show-inheritance:


Retrieving a candidate set of structural homologs via CLI
==========================================================

A candidate set of structural homologs from the [UniRef90](https://www.uniprot.org/help/uniref) clustered database can be retrieved for a set of _seed_ structures using the `retrieve_uniref90_homologs.py` script.

Briefly, the candidate set of homologs is retrieved performing the following steps:
1. Remote Foldseek structural alignment is performed against the UniProt/AlphaFold50 database using the _seed_ structures provided by the user (as PDB files). 
2. The resulting list of candidate hits is filtered based on criteria given by the user (see Parameters section).
3. The list of UniRef90 clusters integrated by candidate hits is retrieved from the UniProtKB server.
4. Sequence IDs for proteins integrating this UniRef90 clusters are retrieved. Sequences are filtered according to their length to avoid including artefactual sequences (see Parameters section).
5. When available, structural models are retrieved from the alphaFold Structure DataBase.  

Importantly, the candidate set of retrieved sequences should be consider a preliminary set of structures to work with. It is highly recommendably for the user to check the quality of the retrieved structures (_i.e._ pLDDT scores) or the taxonomic representation of the set.
Filtering using functions from the `AFDB tools` module might be useful, as well as iterative searches relaxing some of the parameters. When possible, the user might also consider enriching the set with structural models retrieved from other databases (using members of the set as _query_ structures), or with locally infered AlphaFold structural models. 

Usage
======

Minimal usage of the script only requires specifying a folder with _seed_ structures for which a candidate set of homologs is desired.
This is achieved by running the following command-line:

``
./retrieve_uniref90_homologs.py -s seeds_folder
``

where, `seeds_folder` is the path to the folder containing the _seed_ structures. This must be in PDB format and ".pdb" extension.

The user can specify additional parameters in order to constrain the output sets of candidate homolog structures.

As an example, one might want to consider only hits showing at least 70% of both query and target structures during the remote Foldseek structural alignment phase, in order to avoid including spurious hits due to usage of similar domains in non-related structures. Also, in order to further ensure avoiding false positives, a maximum _e-value_ of 1e^{-15} might be specified.
This set of conditions is met by running:

``
./retrieve_uniref90_homologs.py -s seeds_folder -qcov 70 -scov 70 -eval 1e-15 -o ouput_folder
``

Note that the previous line specifies results should be saved in the specified `output_folder`.


Parameters
===========
-s, --seed-structures : str
	Input folder allocating structural models employed in homologues search. Files must be in PDB format, with a .pdb extension.
-o, --output-folder : str
	Output folder where results are saved. [Default: <seed_structures>/output].
-prob, --prob-threshold : float
	Minimum probabily score for a FoldSeek hit to be considered as a candidate homologous structure. [Default: 0.9].
-qcov, --query-cov-threshold : float
	Minimum percentage of the query structure (i.e. seed structure) that must be covered by a Foldseek alignment to consider the hit as a candidate homologous structure. [Default: 70].
-scov, --subject-cov-threshold : float
	Minimum percentage of the subject structure (i.e. candidate homolog) that must be covered by a Foldseek alignment to consider the hit as a candidate homologous structure. [Default: 0].
-eval, --evalue-threshold: float
	Maximum e-value for a FoldSeek hit to be considered as a candidate homologous structure. [Default: 1e-05].
-len_std, --length-std-filtering : float
	Number of length standard deviations to consider around the mean when filtering UniRef90 clusters by sequence length. [Default: 1].

Output
========

For each _seed_ structure provided by the user, an output folder with the candidate set of structural homologs (in PDB format) is provided in an output folder named following the pattern:

``
<seed_name>_prob_<prob>_qcov_<qcov>_scov_<scov>_eval_<eval>_uniref90_homologs_structs
``

Where `seed_name` refers to the prefix in the name of each seed's PDB file, and `prob`, `qcov`, `scov` and `eval` refer to the parameters employed when running the script.
This output folders are by default generated inside a directory named `output` inside the `seed_structures` path. The user can specify alternative output folders by using the `-o` flag when running (see Usage above)

.. automodule:: src.retrieve_uniref90_homologs
    :members:
    :undoc-members:
    :show-inheritance:

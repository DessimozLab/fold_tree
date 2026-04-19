
<p align="center" width="100%">
  <img src="foldtree_logo.png" style="width:30%; height:auto;">
</p>


[![Snakemake](https://img.shields.io/badge/snakemake-≥7.8.0-brightgreen.svg)](https://snakemake.github.io)


Foldtree creates phylogenetic trees from protein structures by using Foldseek to align protein structures and generate distance matrices for tree construction. 

This repository contains the scripts, the main Snakemake workflow, and additional code and data reported in the paper.

# Cite

If you use foldtree in your work, please cite: 

https://www.nature.com/articles/s41594-025-01649-8

DOI: https://doi.org/10.1038/s41594-025-01649-8


# Try on Colab

We created a Colab notebook that allows you to try Foldtree without installing anything. You can find it here:

https://colab.research.google.com/github/DessimozLab/fold_tree/blob/main/notebooks/Foldtree.ipynb

# Install from conda

```
conda install -c bioconda foldtree
```

# Install from sources

Alternatively, you can install Foldtree from sources.

## conda

Foldtree is a Snakemake pipeline that runs multiple tools inside conda environments. First, you need to install conda. We recommend [Miniforge3](https://github.com/conda-forge/miniforge), as it ships both conda and mamba (a faster native replacement for conda).


> [!WARNING]
> 1) Micromamba is not supported, as Snakemake currently does not work reliably with it. If you are using micromamba, we recommend switching to mamba. You *can try* to make it work; good luck. See [one](https://github.com/snakemake/snakemake/issues/3490), [two](https://github.com/snakemake/snakemake/issues/2322).
> 2) We discourage using older version of conda, i.e. versions that do not use libmamba as the default dependency solver. Foldtree creates conda environments under the hood, and solving them with outdated conda may take forever. Check that 
> `conda info | grep solver`
> shows `libmamba`, and update conda if it does not.


## snakemake

Create a conda environment for Foldtree and install Snakemake:

```bash
mamba create -n foldtree python=3.10
mamba activate foldtree
mamba install click snakemake-minimal
```

## foldtree

Clone the repository:

```
git clone git@github.com:DessimozLab/fold_tree.git
cd fold_tree
```

Install Foldtree:

```
python -m pip install . --no-deps --no-build-isolation --no-cache-dir
```

Check that the installation worked:

```
foldtree --help
```


# Usage

## Prepare input
Now we're ready to run the workflow. First, create a folder that will contain the results. Here we call it `myfam`.

```
mkdir myfam
```

Next, choose one of two input options:

### Option 1: UniProt identifiers

Create a text file containing the UniProt identifiers of the structures you would like to include in the tree. Each line should contain one identifier.

Example `identifiers.txt`:
```
A0A060X1A3
A0A087X979
A0A091CQV1
A0A091F6W6
...
```

If identifiers are provided, Foldtree will automatically fetch available protein structures from the AlphaFold DB.

### Option 2: Precompiled structures

If you already have structure files (for example, structures not available in AlphaFold DB), create a `structs` folder inside `myfam` and place all structure files (PDB format) there. In this case, leave the identifier file empty.

```
mkdir myfam/structs
```


## Run Foldtree

All we need to do now is run the workflow. In this example, we use 4 cores. We also have the option of filtering out structures with an average pLDDT < 50 (disabled by default). This was shown to improve the quality of trees relative to sequence-based trees (see the manuscript). Here we disable filtering to include all structures in the tree.

```
foldtree --folder myfam --no-filter -p 
```

We can also use more cores for the Foldseek all-vs-all comparison step using the `--foldseek-cores` flag.

```
foldtree --cores 4 --folder myfam -p --foldseek-cores 4
```

If we have a custom set of structures, use the following command and leave the identifier file empty:

```
foldtree --cores 4 --folder myfam -p --custom-structs
```


Foldtree produces 3 trees: the Foldtree metric, LDDT, and TM tree. Also, it outputs some information retreived from the Uniprot API on the proteins, foldseek comparisons, distance matrices and descriptors of the pLDDT of the structures.

## Run on a cluster

Foldtree passes any additional arguments via `--extra-snakemake-arguments` to snakemake. For example, to run the pipeline on a Slurm cluster, prepare a slurm config 

```
git clone https://github.com/jdblischak/smk-simple-slurm slurmsimple
```

Adapt `slurmsimple/simple/config.yaml` correspondingly. Then pass it to snakemake:

```
foldtree <other args> -esa "--profile slurmsimple/simple"
```

## Managing Foldtree conda environments

To run Foldtree, Snakemake automatically creates additional conda environments for the required tools. By default, these environments are stored in `~/.foldtree`.

The environments are created only once (about 2 GB total) and reused in subsequent runs.

On some systems this may be inconvenient (for example, on HPC clusters the home directory might have limited storage). In such cases, you can specify an alternative location using `--conda-prefix`:

```
foldtree <other-args> --conda-prefix <path-to-env-storage>
```



# Known issues

## "libmamba Non-conda folder exists at prefix"

You are likely using an older version of mamba (<2.4). Update mamba or force snakemake to use conda (should not take much longer):

```
foldtree <other args> -esa "--conda-frontend conda"
```

The [legends say](https://github.com/snakemake/snakemake/issues/3249) you can also use `mamba=1`, but we have not tested it.


# Other
## Repo contents

In the examples folder you can find the data we used to make the RRNPPA phylogeny which is presented in detail in the manuscript as a motivating case study alongside other case studies.

In the notebooks folder there are jupyter notebooks used to collect and plot the results of the benchmarking experiments we performed as well as generate the RRNPPA phylogeny.

The src folder contains python code used to interact with the Uniprot and Alphafold databases as well as building structural.


## Benchmarking experiments

### Prepare

To reproduce the benchmarking experiments, use the benchmarking workflow. This workflow is not directly wrapped by the `foldtree` command. To run it, you need to run it with snakemake from sources.

The experiment should be done on a cluster environment, since thousands of trees will be generated. We performed experiments on a slurm cluster and the `slurmsimple` module is included as a module in this repo to allow snakemake to schedule jobs. Other cluster approaches to using snakemake on a cluster should also work (snakemake provides extensive documentation on this). 

As the benchmarking pipeline has more dependencies than the main one, it was isolated into a different conda environment. To install benchmarking dependencies, run:

```
mamba env create -n foldtree-bench --file=src/foldtree/workflow/config/foldtree.yaml
mamba activate
```

### Run

The benchmarking pipeline is located in `src/foldtree/workflow`. It will output 4 types of trees for sequence-based analysis for each family using 2 possible aligners (`muscle5` or `clustalo`) and 2 possible tree building approaches (`iqtree` of `fasttree`). The pipeline also output 12 structural trees for each family. These either use only 3di or 3di and amino acid alignments, 3 different strutural distances ( Fident or Foldtree metric, LDDT and TM score) and with or without statistical correction. The trees are then rooted with MAD and scored using ultrametricity and taxonomic congruence metrics.

To rerun the OMA HOG experiments download the identifier data from zenodo. (https://doi.org/10.5281/zenodo.8346286)

Then unzip the archives to create the folder structure

Then use snakemake on each directory of the data. Here we only derive trees for families starting at the LUCA level.

```
snakemake  --profile slurmsimple/simple --use-conda -s ./workflow/Benchmarking_workflow --directory  ./path/to/your/data/OMA_data/LUCA
```

To rerun the CATH/CAT experiments you will need to install openMM to deal with discontinuous structures and the Pebble multiprocessing librarx. These packages can be install in the foldree environment

Let's activate the environment and install openmm, pdbfixer and pebble

```
mamba activate foldtree
mamba install -c conda-forge openmm
mamba install -c conda-forge pdbfixer
mamba install pebble
```

Now we can run the data preparation script.

```
python ./src/dataprep/prepare_protsets_CATH.py
```

This will set up folders for each family at the CAT and CATH levels, download PDB files for each of the families and correct them using PDBfixer from OpenMM.
Now you can run the workflow on the CATH or CAT dataset in the same way we did with the OMA dataset.

```
snakemake  --profile slurmsimple/simple --use-conda -s ./workflow/Benchmarking_workflow --directory  ./path/to/your/data/CATH_data/
```

The code to generate the figures from the manuscript is available in the notebooks folder of the repo.

## Funding

Funded by NIH* through the Pathogen Data Network. 
<p align="center" width="100%">
  <img src="pdn_logo.png" style="width:30%; height:auto;">
</p>

*This resource is supported as a whole or in part by the National Institute Of Allergy And Infectious Diseases of the National Institutes of Health under grant n°U24AI183840, awarded to the SIB Swiss Institute of Bioinformatics. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.

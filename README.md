# Snakemake workflow: `snake_tree`
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)

A Snakemake workflow for building stupid structure trees woohooo... finish this

## Usage

install snakemake
Clone repo
cd to repo 

download and install up to date foldseek binaries

## Installation

    # static Linux AVX2 build (check using: cat /proc/cpuinfo | grep avx2)
    wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz; tar xvzf foldseek-linux-avx2.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH
    # static Linux SSE4.1 build (check using: cat /proc/cpuinfo | grep sse4_1)
    wget https://mmseqs.com/foldseek/foldseek-linux-sse41.tar.gz; tar xvzf foldseek-linux-sse41.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH
    # static macOS build (universal binary with SSE4.1/AVX2/M1 NEON)
    wget https://mmseqs.com/foldseek/foldseek-osx-universal.tar.gz; tar xvzf foldseek-osx-universal.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH

If you dont have root privileges then make sure the binary is located at foldtree/foldseek/bin/foldseek

test with

snakemake .testdata/ --cores 1 --use-conda snakemake  ./testdata/RFdistances.json --use-conda --cores 3

This will generate struct trees based on different alignment modes and scores as well as the consensus tree and compare their pairwise RF distances. The structures will be downloaded from alpha fold db and put into a structs folder. The sequence_dataset.csv file will be generated from the uniprot entries of each protein. It contains lineage information, the protein description and its amino acid sequence. The allvall_*_.csv files contain the foldseek all vs all comparison results for alignment modes 0 and 2. 


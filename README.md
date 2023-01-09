# Snakemake workflow: `snake_tree`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for building stupid structure trees woohooo... finish this


## Usage

install snakemake
test with
snakemake --cores 1 --use-conda

this should use the testdata folder to generate a structure and sequence tree comparison.
otherwise generate just the struct tree using 

snakemake .testdata/struct_tree.nwk.PP.nwk.rooted --cores 1 --use-conda

This should automatically root it using MAD

snakemake -s ./workflow/benchmark_all --profile slurmsimple/simple/ --config custom_structs=True cath=True filter=False iqtree_cores=1 iqtree_redo=True -k --touch --rerun-incomplete --use-conda --directory ./CAT_data

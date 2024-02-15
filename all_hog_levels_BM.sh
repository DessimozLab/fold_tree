#!/bin/bash

program="snakemake -s ./workflow/benchmark_all --profile slurmsimple/simple/ -T 5 --config filter=False -k --rerun-incomplete --use-conda --directory   "
folders=($(ls -d /work/FAC/FBM/DBC/cdessim2/default/dmoi/projects/snake_tree/OMA_data_unfiltered/OMA_data/* ) )
for folder in "${folders[@]}"; do
  echo $program "$folder"
  $program "$folder" --unlock
  $program "$folder" 
done



#!/bin/bash

program="snakemake --profile slurmsimple/simple/ --config filter=True -k --rerun-incomplete -T 10 --use-conda --directory  "
folders=($(ls -d /work/FAC/FBM/DBC/cdessim2/default/dmoi/projects/snake_tree/OMA_data/* ) )
for folder in "${folders[@]}"; do
  echo $program "$folder"
  $program "$folder"
done

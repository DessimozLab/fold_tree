#!/bin/bash

program="snakemake --profile slurmsimple/simple/ --config filter=True -k --rerun-incomplete --use-conda --touch --directory  "
unlock="snakemake --profile slurmsimple/simple/ --config filter=True -k --rerun-incomplete  --use-conda --unlock --directory  "

folders=($(ls -d /work/FAC/FBM/DBC/cdessim2/default/dmoi/projects/snake_tree/OMA_data/* ) )
for folder in "${folders[@]}"; do
  echo $unlock "$folder"
  $unlock "$folder"

  echo $program "$folder"
  $program "$folder"
done

#!/bin/bash

program="snakemake --profile slurmsimple/simple/  --use-conda --directory  "
folders=($(ls -d /work/FAC/FBM/DBC/cdessim2/default/dmoi/datasets/Structure_Trees_mk2/* ) )
for folder in "${folders[@]}"; do
  echo $program "$folder"
  $program "$folder"
done
#!/bin/bash

# The name of the program to run
program="./snakemake --profile slurmsimple/simple/  --use-conda --directory  "

# Get the list of folders using ls and store it in an array
folders=($(ls -d /work/FAC/FBM/DBC/cdessim2/def│ault/dmoi/datasets/Structure_Trees/))

for folder in "${folders[@]}"; do
  # Remove the trailing slash from each folder name
  # Use the folder variable in the program
  $program "$folder"

  # Go back to the parent folder
done
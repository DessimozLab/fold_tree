#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=10G
#SBATCH --time=24:00:00
#SBATCH --job-name=python
source /work/FAC/FBM/DBC/cdessim2/default/dmoi/miniconda3/etc/profile.d/conda.sh 
#conda init bash
conda activate ML2
python prepare_protsets_CATH.py
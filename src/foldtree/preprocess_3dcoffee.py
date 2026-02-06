# preprocess 3dcoffee alignment for nextflow pipeline input
# create samplesheet.csv and toolsheet.csv in folder 3dcoffee
# pdb folder too! No! they are in {wildcards.folder}/structs/
import os
from Bio import SeqIO
import glob

folder = snakemake.wildcards.folder
structfolder = folder + '/structs/' 
print('structfolder: ', structfolder)

pdbs = glob.glob(structfolder + '*.pdb')

#find absolute path of pdbs
pdbs = [os.path.abspath(pdb) for pdb in pdbs]

with open( snakemake.output[0], 'w') as outfile:
    for pdb in pdbs:
        pdbid = os.path.basename(pdb).split('.')[0].split('/')[-1]
        outfile.write('>'+pdbid + ' _P_ ' + pdb + '\n')
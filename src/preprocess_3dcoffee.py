# preprocess 3dcoffee alignment for nextflow pipeline input
# create samplesheet.csv and toolsheet.csv in folder 3dcoffee
# pdb folder too! No! they are in {wildcards.folder}/structs/
import os

toolsheet_str = """tree,args_tree,aligner,args_aligner
,,3DCOFFEE,"""
samplesheet_str = """id,fasta,reference,structures
test,{0},{0},{1}""".format(snakemake.input[0], snakemake.wildcards.folder)  # check how to use wildcards in script

# write files
with open(snakemake.output[0]) as samplefile:
    samplefile.write(samplesheet_str)
with open(snakemake.output[1]) as toolfile:
    toolfile.write(toolsheet_str)

# convert clustal format to fasta and save in output
from Bio import SeqIO
in_alg = SeqIO.read(snakemake.input[0],format='clustal')
SeqIO.write(in_alg,snakemake.output[0],format='fasta')

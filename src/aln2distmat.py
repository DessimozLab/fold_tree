
import toytree
import structalns
import os
import pandas as pd
import glob
from Bio import SeqIO
import structalns
from foldseek2tree import distmat_to_txt
import numpy as np

#take a fasta as input and transform the msa into a distance matrix
#meant as a control against ml methods

print(snakemake.input)
Fasta = snakemake.input[0]
#parse fasta
fasta = SeqIO.parse(Fasta, 'fasta')
sequences = {i.id:i.seq for i in fasta}
#separate dictionary into two lists
seqnames = list(sequences.keys())
seqs = list(sequences.values())
#calculate all pairwise identities
distmat = np.zeros((len(sequences), len(sequences)))
for i in range(len(seqs)):
    for j in range(len(seqs)):
        if i < j:
            distmat[i,j] = 1-structalns.Fident(seqs[i], seqs[j])
distmat = distmat + distmat.T
#output the distmat
distmat_to_txt( seqnames, distmat, snakemake.output[0] )

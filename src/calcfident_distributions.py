import structalns
import itertools
from Bio import SeqIO
from scipy.stats import describe

alnstats = {}

for alnfile in snakemake.input[1:]:
    aln = SeqIO.parse(alnfile, 'fasta')
    sequences = {i.id:i.seq for i in aln}
    #calculate all pairwise identities    
    fident_vec = { (i,j): structalns.Fident(sequences[i], sequences[j]) for i,j in itertools.combinations(sequences.keys(),2)}
    
    description = describe(list(fident_vec.values()))
    #add the frac ident for each alignment type to the dataframe

    
    alnstats[alnfile] = description
    print(alnfile, description)



#output the stats as json
import json

with open(snakemake.output[0] , 'w') as out:
    json.dump(alnstats, out)

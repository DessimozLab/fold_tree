
import toytree
import structalns
import os
import pandas as pd
import glob


print(snakemake.input)
print(snakemake.output)

alndf = pd.read_table(snakemake.input[0], header = None)

infolder = snakemake.input[0].split('/')[:-1]
infolder = ''.join( [i + '/' for i in infolder])

mapper3di, mapperAA = structalns.read_dbfiles3di( snakemake.input[1] , snakemake.input[2])

#add the 3di alignment to the dataframe
columns = 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,qaln,taln,cigar,alntmscore'.split(',')
alndf.columns = columns

print('submat', snakemake.params.submat)


alndf['query'] = alndf['query'].map(lambda x :x.replace('.pdb', ''))
alndf['target'] = alndf['target'].map(lambda x :x.replace('.pdb', ''))

alndf['3diq']= alndf['query'].map(mapper3di)
alndf['3dit']= alndf['target'].map(mapper3di)
alndf['AAq']= alndf['query'].map(mapperAA)
alndf['AAt']= alndf['target'].map(mapperAA)

#output a fasta with the 3di sequences
res = alndf.apply(structalns.calc_fident_crossaln , axis = 1)
alndf = pd.concat([alndf,res] , axis = 1)

with open(snakemake.output[0] , 'w') as out:
    for seq in alndf['query'].unique():
        out.write('>'+seq.replace('.pdb', '' )+'\n')
        out.write(mapperAA[seq.replace('.pdb', '')]+'\n')

with open(snakemake.output[1] , 'w') as out:
    for seq in alndf['query'].unique():
        out.write('>'+seq.replace('.pdb', '' )+'\n')
        out.write(mapper3di[seq.replace('.pdb', '')]+'\n')

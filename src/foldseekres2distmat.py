import foldseek2tree
import numpy as np
import pandas as pd

res = pd.read_table(snakemake.input[0], header = None)
print(res.head())

res[0] = res[0].map(lambda x :x.replace('.pdb', ''))
res[1] = res[1].map(lambda x :x.replace('.pdb', ''))
res.columns = 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,lddtfull,alntmscore'.split(',')
ids = list( set(list(res['query'].unique()) + list(res['target'].unique())))
pos = { protid : i for i,protid in enumerate(ids)}
kernels = ['fident', 'alntmscore', 'lddt']
matrices = { k:np.zeros((len(pos), len(pos))) for k in kernels }
#calc kernel for tm, aln score, lddt
for idx,row in res.iterrows():
    for k in matrices:
        matrices[k][pos[row['query']] , pos[row['target']]]= row[k]

for i,kernel in enumerate(matrices):
    matrices[k] += matrices[k].T
    matrices[k] /= 2
    np.save(k + '_distmat.npy' , matrices[k])
    distmat_txt = foldseek2tree.distmat_to_txt( ids , matrices[k] , snakemake.output[i] )
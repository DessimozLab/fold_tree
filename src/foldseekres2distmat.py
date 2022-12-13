import foldseek2tree
import numpy as np
import pandas as pd

res = pd.read_table(snakemake.input[0], header = None ,delim_whitespace=True)
res[0] = res[0].map(lambda x :x.replace('.pdb', ''))
res[1] = res[1].map(lambda x :x.replace('.pdb', ''))
ids = list( set(list(res[0].unique()) + list(res[1].unique())))
pos = { protid : i for i,protid in enumerate(ids)}
self_dists = res[res[0] == res[1]]
self_distmap = dict(zip(self_dists[0] , self_dists[2] ) )    
kernel_distmat = np.zeros((len(pos), len(pos)))
for idx,row in res.iterrows():
    kernel_distmat[pos[row[0]] , pos[row[1]]] = foldseek2tree.kernelfun(self_distmap[row[0]] , self_distmap[row[1]] , row[2])
kernel_distmat += kernel_distmat.T
kernel_distmat /= 2

distmat_txt = foldseek2tree.distmat_to_txt( ids , kernel_distmat , snakemake.output[0] )

import foldseek2tree
import numpy as np
import pandas as pd




res = pd.read_table(snakemake.input[0], header = None ,delim_whitespace=True)
res[0] = res[0].map(lambda x :x.replace('.pdb', ''))
res[1] = res[1].map(lambda x :x.replace('.pdb', ''))
res.columns = 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,lddtfull,alntmscore'.split(',')


ids = list( set(list(res[0].unique()) + list(res[1].unique())))
pos = { protid : i for i,protid in enumerate(ids)}
self_dists = res[res[0] == res[1]]
self_distmap = dict(zip(self_dists[0] , self_dists[2] ) )    
kernel_distmat = np.zeros((len(pos), len(pos)))
distmat = np.zeros((len(pos), len(pos)))

for distance in ['lddt' , 'bits' , 'fident' , 'alntmscore' , 'evalue' ]:
    for idx,row in res.iterrows():
    #calc kernel for tm, aln score, lddt
        kernel_distmat[pos[row[0]] , pos[row[1]]] = foldseek2tree.kernelfun(self_distmap[row[0]] , self_distmap[row[1]] , row[2])
        distmat[pos[row[0]] , pos[row[1]]]= row[3]

        kernel_distmat += kernel_distmat.T
        kernel_distmat /= 2

        distmat += distmat.T
        distmat /= 2

        #0 to 1, lower is better
        if distance in [  'evalue' ]:
            
        #0 to 1, higher is better
        if distance in [ 'alntmscore' , 'lddt' ]:
            distmat = 1 - distmat
            kernel_distmat = 1 - kernel_distmat
        
        #0 to N , higher is better
        if distance in [ 'bits' ]:
            distmat /= np.max(distmat)
            kernel_distmat /= np.max(kernel_distmat)
            distmat = 1 - distmat
            kernel_distmat = 1 - kernel_distmat
        
        np.save(distance + '_kernel_distmat.npy' , kernel_distmat)
        np.save(distance + '_distmat.npy' , distmat)
        
        distmat_txt = foldseek2tree.distmat_to_txt( ids , distmat , snakemake.output[0]  )
        distmat_txt = foldseek2tree.distmat_to_txt( ids , distmat , snakemake.output[0]  )

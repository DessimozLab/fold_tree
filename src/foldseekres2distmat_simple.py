import foldseek2tree
import numpy as np
import pandas as pd


try:
    res = pd.read_table(snakemake.input[0], header = None)

    print(res.head())
    #get the folder of the input file
    infolder = snakemake.input[0].split('/')[:-1]
    infolder = ''.join( [i + '/' for i in infolder])+'/'
    res[0] = res[0].map(lambda x :x.replace('.pdb', ''))
    res[1] = res[1].map(lambda x :x.replace('.pdb', ''))
    res.columns = 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,lddtfull,alntmscore'.split(',')
    ids = list( set(list(res['query'].unique()) + list(res['target'].unique())))
    pos = { protid : i for i,protid in enumerate(ids)}
    kernels = ['fident', 'alntmscore', 'lddt' ]

    #set kernel columns to float
    for k in kernels:
        res[k] = res[k].astype(float)
    #change nan to 0
    res = res.fillna(0)
    matrices = { k:np.zeros((len(pos), len(pos))) for k in kernels }
    print(res)
    #calc self score

    #calc kernel for tm, aln score, lddt
    for idx,row in res.iterrows():
        for k in matrices:
            matrices[k][pos[row['query']] , pos[row['target']]] += row[k]
            matrices[k][pos[row['target']] , pos[row['query']]] += row[k]

    for i,k in enumerate(matrices):
        matrices[k] /= 2
        matrices[k] = 1-matrices[k]
        np.fill_diagonal(matrices[k], 0)    
        print(k , matrices[k], 'variance' ,  np.var(matrices[k]) , 'maxdist' , np.amax(matrices[k]))
        np.save( infolder + k + '_distmat.npy' , matrices[k])
        print(snakemake.output[i]  , k )

        if snakemake.params.varfilter == True:
            if np.amax(matrices[k]) > 0.35 and np.var(matrices[k]) > 0.0005 :
                distmat_txt = foldseek2tree.distmat_to_txt( ids , matrices[k] , snakemake.output[i] )
            else:
                print('dist too low, not writing distmat')
                with open(snakemake.output[i] , 'w') as handle:
                    handle.write('')
                    handle.write('\n')
        else:
            distmat_txt = foldseek2tree.distmat_to_txt( ids , matrices[k] , snakemake.output[i] )
        
except pd.errors.EmptyDataError:
    print(' err empty file found for all vs all')
    for o in snakemake.output:
        with open(o, 'w') as handle:
            handle.write('')
    raise SystemExit(0)
import pandas as pd
import Bio.PDB
import numpy as np
import os
import tqdm

def extract_core(resdf , hitthresh = .8 , corefolder = 'core_structs/' , structfolder = 'structs/'):
    #read all results
    folder =''.join([ sub + '/' for sub in resdf.split('/')[:-1] ])
    print(folder)
    resdf = pd.read_table(resdf , header = None)
    #map hits to each struc
    hits = {}
    #proportion of structures that need to map to a residue fo
    with tqdm.tqdm(total=len(resdf[0].unique())) as pbar:
        for i,q in enumerate(resdf[0].unique()):
            sub = resdf[resdf[0] == q]
            hitvec = np.zeros((1,sub.iloc[0][7]))
            for idx,r in sub.iterrows():
                hitvec[0,r[5]:r[6]] = hitvec[0,r[5]:r[6]]+1
            hitvec /= len(resdf[0].unique())
            core = np.where(hitvec>hitthresh)[1]
            hits[q]= { 'min': np.amin(core), 'max': np.amax(core)}
            pbar.set_description('processed: %d' % (1 + i))
            pbar.update(1)
        
    #make core struct folder
    try:
        os.mkdir(folder+corefolder)
    except:
        print('already')

    #parse each struct and output core to folder 
    parser = Bio.PDB.PDBParser()
    with tqdm.tqdm(total=len(hits)) as pbar:
        for i,q in enumerate(hits):
            print(folder+structfolder+q)
            struct = parser.get_structure(q.split('.')[0], folder+structfolder+q )
            #zero based indexing...
            struct_core = Bio.PDB.Dice.extract( struct ,'A' , hits[q]['min']+1 , hits[q]['max']+1 ,folder+corefolder+q  )
            pbar.set_description('processed: %d' % (1 + i))
            pbar.update(1)
    hitsdf = pd.DataFrame.from_dict( hits , orient='index'  )
    hitsdf.to_csv(folder +'struct_cores.csv')
    return folder +'struct_cores.csv'



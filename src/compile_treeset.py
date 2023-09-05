
import glob
import pandas as pd
import json
import os
import tqdm
from itertools import combinations
import plotly.express as px
import plotly.figure_factory as ff
from scipy.stats import wilcoxon
import numpy as np


def compile_folder_resdict(rootfolder , scorefunc = 'score_x_frac' , verbose = False):

    """
    this function compiles the treescores for all the trees in a folder
    it checks that the number of sequences and structures are the same
    it also checks that the treescores are present

    params:
        rootfolder: folder with the treescores
        scorefunc: treescore to use
    returns:
        resdf: dataframe with the treescores
        refclols: list of the columns in the dataframe

    """
    
    print(rootfolder)
    res = {}
    folders = set(glob.glob(rootfolder + '*/' ))-set([rootfolder+'logs/'])
    print(len(folders))
    if len(folders)>0:
        with tqdm.tqdm(total=len(folders)) as pbar:
            for i,folder in enumerate(folders):
                nstructs = len(glob.glob(folder+'structs/*.pdb'))
                if os.path.isfile(folder+'treescores_sequences.json'):
                    treescores = glob.glob(folder + '*_treescores_struct_tree.json' ) + list(glob.glob(folder+'treescores_sequences*.json')) + list(glob.glob( folder+'treescores_sequences_iq*.json'))
                    if len(treescores)>0 and os.path.isfile(folder + 'sequences.fst'):
                        with open(folder + 'sequences.fst') as fstin:
                            nseqs = fstin.read().count('>')
                        pbar.set_description('processed: %d' % (1 + i))
                        pbar.update(1)
                        if nseqs == nstructs :
                            for score in treescores:
                                #check if score exists
                                if os.path.isfile(score):

                                    with open(score) as taxin:
                                        tax_res = json.load(taxin)
                                    tax_res= {s.split('/')[-1]:tax_res[s] for s in tax_res}
                                    if folder not in res:
                                        res[folder] = { s:tax_res[s][scorefunc] for s in tax_res if  scorefunc  in tax_res[s]}
                                    else:
                                        res[folder].update({ s:tax_res[s][scorefunc] for s in tax_res if scorefunc in tax_res[s]})

                            res[folder].update({ 'nseqs':   nseqs})
                        else:
                            if verbose == True:
                                print('nseqs != nstructs', folder)
                                print(nseqs, nstructs)
            return res
                            
def compile_folder(rootfolder , scorefunc = 'score_x_frac' , verbose = False):
    '''
    This function compiles the treescores for all the trees in a folder
    it checks that the number of sequences and structures are the same
    it also checks that the treescores are present
    
    '''

    res = compile_folder_resdict(rootfolder , scorefunc = scorefunc , verbose = verbose)
    if len(res)>0:
        resdf = pd.DataFrame.from_dict(res, orient = 'index')
        resdf.columns = [ c.replace('.PP.nwk.rooted', '').replace('.aln.fst.nwk.rooted' , '' ) for c in  resdf.columns]
        if verbose == True:
            print(resdf.head(), resdf.shape)
        refcols = list(resdf.columns)
        refcols.remove('nseqs')

        #divide the scores by the number of sequences
        for c in refcols:
            resdf[c+'_norm'] = resdf[c] / resdf['nseqs']

        for c1,c2 in combinations(refcols,2):
            resdf[c1+'_'+c2+'_delta'] = resdf[c1] - resdf[c2]
            resdf[c1+'_'+c2+'_max'] = resdf[[c1,c2]].apply( max , axis = 1)

            resdf[c1+'_'+c2+'_delta_norm'] = resdf[c1+'_'+c2+'_delta'] / resdf[c1+'_'+c2+'_max']
        resdf['clade'] = rootfolder.split('/')[-2]
        resdf['family'] = resdf.index.map( lambda x :  x.split('/')[-2])

        return resdf, refcols

def compile_folder_treestats(rootfolder , scorefunc = 'ultrametricity_norm' , verbose = False):
    res = compile_folder_resdict(rootfolder , scorefunc = scorefunc , verbose = verbose)
    if len(res)>0:
        resdf = pd.DataFrame.from_dict(res, orient = 'index')
        #resdf.columns = [ c.replace('.PP.nwk.rooted', '').replace('.aln.fst.nwk.rooted' , '' ) for c in  resdf.columns]
        if verbose == True:
            print(resdf.head(), resdf.shape)
        return resdf
    
def compare_treesets(tree_resdf , colfilter= 'sequence' , display_lineplot = False , display_distplot = True , verbose = False):

    '''
    this function compares the treescores for all the trees in a folder
    it uses ploty to plot the results and performs a wilcoxon test

     params:
        tree_resdf: dataframe with the treescores
        colfilter: string to filter the columns to compare
    returns:
        None
    '''
    rescols = [ 'lddt_1_raw_struct_tree' , 'fident_1_raw_struct_tree', 'alntmscore_1_raw_struct_tree',   'sequences' ]
    refcols = [ 'lddt_1_raw_struct_tree' , 'fident_1_raw_struct_tree', 'alntmscore_1_raw_struct_tree','lddt_0_raw_struct_tree' , 'fident_0_raw_struct_tree', 'alntmscore_0_raw_struct_tree', 'sequences' ]


    for c1,c2 in combinations(refcols,2):
        if colfilter in c1 or colfilter in c2:
            try:
                if verbose == True:
                    print(c1,c2)
                    print('delta:', tree_resdf[c1+'_'+c2+'_delta'].dropna().sum(),
                        'delta norm:',   tree_resdf[c1+'_'+c2+'_delta_norm'].dropna().sum(),wilcoxon(tree_resdf[c1+'_'+c2+'_delta'].dropna()))
            except:
                print('error', c1,c2)

            try:
                sub = tree_resdf.dropna(subset = [c1+'_'+c2+'_delta'])
                maxval = sub[[c1, c2]].max().max()
                if display_lineplot == True:
                    fig = px.scatter(sub, x=c1, y=c2 , hover_data=[c1+'_'+c2+'_delta_norm' , c1+'_'+c2+'_delta'  , 'family'])
                    fig.add_shape(type="line",
                        x0=0, 
                        y0=0, 
                        x1=maxval, 
                        y1=maxval)
                    fig.show()
            except:
                print('error', c1,c2)
    
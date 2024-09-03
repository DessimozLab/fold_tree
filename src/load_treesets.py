#!/usr/bin/env python
# coding: utf-8

import glob
import os
import numpy as np
import sys
sys.path.append( '../src/')
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
from compile_treeset import compile_folder , compare_treesets
import plotly.express as px
import plotly.figure_factory as ff
import plotly.express as px
from scipy.stats import wilcoxon, ks_2samp , describe

#construct structure and sequence feature dataset
import tqdm
import glob
import os
import numpy as np
import tqdm
import json

from compile_treeset import compile_folder_treestats
import pandas as pd


def compile_results(cladefolders , plotall = False, scorefunc = 'root_score'):
    rescols = ['lddt_1_raw_struct_tree' , 'fident_1_raw_struct_tree', 'alntmscore_1_raw_struct_tree', 'sequences' , 'sequences.aln.fst.treefile.rooted']
    rescols_norm = ['lddt_1_raw_struct_tree'+'_norm' , 'fident_1_raw_struct_tree'+'_norm', 'alntmscore_1_raw_struct_tree'+'_norm' , 'sequences'+'_norm' , 'sequences.aln.fst.treefile.rooted']
    dfs = []
    plotall = False
    for folder in cladefolders:
        #try:
        if 'logs' not in folder:
            print(folder)
            res = compile_folder(folder, scorefunc = scorefunc, verbose = True)
            if res :
                tree_resdf , refcols = res
                compare_treesets(tree_resdf  , colfilter= 'sequence' , display_lineplot = False , verbose = True)
                tree_resdf['folder'] = folder

                #add filtered bool column
                tree_resdf['filtered'] = tree_resdf['folder'].apply(lambda x : 'unfiltered' not in x)

                if 'OMA' in folder:
                    dfs.append(tree_resdf)
                    if plotall == True:
                        graph_treedf(tree_resdf , rescols, rescols_norm)
                else:
                    dfs.append(tree_resdf)
    total_df = pd.concat(dfs)
    try:
        if 'OMA' in folder:
            graph_treedf(total_df[total_df.filtered == False] , rescols, rescols_norm , prefix = 'OMA' )
            graph_treedf(total_df[total_df.filtered == True] , rescols, rescols_norm , prefix = 'OMA' )
    except:
        print('graphing err' )
    return total_df


cladefolders = [ '../CAT_mk3/']# , '../SCOP_data/' ]
catdf = compile_results(cladefolders, scorefunc = 'root_score_nr' )
catdf_ntcs = compile_results(cladefolders, scorefunc = 'score' )
catdf_lineage = compile_results(cladefolders, scorefunc = 'taxdegree_score' )


cladefolders = [ '../CATH_mk3/']# , '../SCOP_data/' ]
cathdf_root = compile_results(cladefolders, scorefunc = 'root_score_nr' )
cathdf_ntcs = compile_results(cladefolders, scorefunc = 'score' )
cathdf_lineage = compile_results(cladefolders, scorefunc = 'taxdegree_score' )

cladefolders = set(glob.glob( '../OMA_data_unfiltered/OMA_data/*/' )) - set([ '../OMA_data_unfiltered/OMA_data/logs/' ]) 
OMADF = compile_results(cladefolders, scorefunc = 'root_score_nr' )
OMADF_ntcs = compile_results(cladefolders, scorefunc = 'score' )
OMADF_lineage = compile_results(cladefolders, scorefunc = 'taxdegree_score' )

#set columns to the intersection of all the dataframes
#find intersection of all OMADF colums

omacols = set(OMADF.columns).intersection(OMADF_ntcs.columns).intersection(OMADF_lineage.columns)
catcols = set(catdf.columns).intersection(catdf_ntcs.columns).intersection(catdf_lineage.columns)
cathcols = set(cathdf_root.columns).intersection(cathdf_ntcs.columns).intersection(cathdf_lineage.columns)

OMADF = OMADF[omacols]
OMADF_ntcs = OMADF_ntcs[omacols]
OMADF_lineage = OMADF_lineage[omacols]

catdf = catdf[catcols]
catdf_ntcs = catdf_ntcs[catcols]
catdf_lineage = catdf_lineage[catcols]

cathdf_root = cathdf_root[cathcols]
cathdf_ntcs = cathdf_ntcs[cathcols]
cathdf_lineage = cathdf_lineage[cathcols]

#save cat , cath and omadf
catdf.to_pickle('catdf.pkl')
catdf_ntcs.to_pickle('catdf_ntcs.pkl')
catdf_lineage.to_pickle('cathf_lineage.pkl')

cathdf_root.to_pickle('cathdf.pkl')
cathdf_ntcs.to_pickle('cathdf_ntcs.pkl')
cathdf_lineage.to_pickle('cathdf_lineage.pkl')

OMADF.to_pickle('omadfmk2.pkl')
OMADF_ntcs.to_pickle('omadfmk2_ntcs.pkl')
OMADF_lineage.to_pickle('omadfmk2_lineage.pkl')


OMADF = OMADF[OMADF.filtered == False]
unfiltereddf = OMADF
unfiltereddf_ntcs = OMADF_ntcs
unfiltereddf_lineage = OMADF_lineage

cladefolders = set(glob.glob( '../OMA_data_unfiltered/OMA_data/*/' ))-set([ '../OMA_data/logs/' ])

dfs = []
for folder in cladefolders:
    print(folder)
    OMA_treestat_DF = compile_folder_treestats( folder, verbose = False  )
    dfs.append(OMA_treestat_DF)
OMA_treestats = pd.concat(dfs)
OMA_treestats.to_pickle('OMA_treestats.pkl')

cladefolders = [ '../CAT_mk3/' , '../CATH_data/' ]
dfs = []
for folder in cladefolders:
    print(folder)
    CATH_treestat_DF = compile_folder_treestats( folder, verbose = False  )
    dfs.append(CATH_treestat_DF)
CATH_treestats = pd.concat(dfs)
CATH_treestats.to_pickle('CATH_treestats.pkl')

allfolders = []
cladefolders = set(glob.glob( '../OMA_data_unfiltered/OMA_data/*/' )) - set([ '../OMA_data_unfiltered/OMA_data/logs/' ])
for folder in cladefolders:
    hogs = glob.glob(folder+'*/')
    cladestats = []
    print(folder)
    allstats = {}
    for h in tqdm.tqdm(hogs):
        alns = glob.glob(h+'*dist_*_*.json')
        for a in alns:
            allstats[h+a] = {}
            with open(a) as f:
                alnstats = json.load(f)
            allstats[h+a]['fam'] = h
            allstats[h+a]['aln'] = a.split('/')[-1]
            allstats[h+a]['stats'] = list( alnstats.values() )[0]
    cladestatsdf = pd.DataFrame.from_dict(allstats , orient = 'index')
    allfolders.append(cladestatsdf)
allfoldersdf = pd.concat(allfolders)
#save allfoldersdf
allfoldersdf.to_pickle('alnstats_allfoldersdf.pkl')


def ret_description(vec , label = ''):
    return { label+'_mean': np.mean(vec) , label+'_max': np.amax(vec) , label+'_min': np.amin(vec) ,label+'_var': np.var(vec) }

#cladefolders = set(glob.glob( '../OMA_data/*/' ))-set([ '../OMA_data/logs/' ])
cladefolders =  set(glob.glob( '../OMA_data_unfiltered/OMA_data/*/' )) - set([ '../OMA_data_unfiltered/OMA_data/logs/' ])
allfolders = cladefolders

dfs = []
print(allfolders)
for clade in allfolders:
    print(clade)
    res = {}
    folders = glob.glob(clade + '*/' )
    #for folder in tqdm.tqdm_notebook(folders):
    for folder in folders:
        if 'logs' not in folder:
            
            nstructs = len(glob.glob(folder+'structs/*.pdb'))
            try:
                with open(folder + 'sequences.fst') as fstin:
                    nseqs = fstin.read().count('>')
            except:
                nseqs = 0
            if nstructs == nseqs:
                if os.path.isfile(folder + 'plddt.json' ) and os.path.isfile(folder + 'sequence_dataset.csv' ):
                    plddt_df = pd.read_json(folder + 'plddt.json').T
                    if len(plddt_df)>0:
                        plddt_df.columns = 'nobs,minmax,mean,variance,skewness,kurtosis'.split(',')
                        plddt_df['min'] = plddt_df.minmax.map( lambda r: r[0] )
                        plddt_df['max'] = plddt_df.minmax.map( lambda r: r[1] )
                        res[folder] = {}
                        for col in ['nobs', 'min' , 'max' , 'mean' , 'variance' , 'skewness' , 'kurtosis' ]:
                            descriptors = ret_description(plddt_df[col] , label = col)
                            for l in descriptors:
                                res[folder][l] = descriptors[l]
                        #add in some descriptors of the taxonomic spread and sequence set
                        
                        seqdf = pd.read_csv(folder+'sequence_dataset.csv' )
                        res[folder]['nprots'] = len(seqdf)
                        cladesets = [ set(l.split(',')) for l in  seqdf['Taxonomic lineage (Ids)'] ]
                        union_all = cladesets[0]
                        intersection_all = cladesets[0]
                        for c in cladesets:
                            union_all = union_all.union(c)
                            intersection_all.intersection(c)
                        res[folder]['n_clades'] = len(union_all)
                        res[folder]['inter_clades'] = len(intersection_all)
                        res[folder]['inter/nc'] =  res[folder]['inter_clades'] /  res[folder]['n_clades'] 
                        res[folder]['nc/np'] = res[folder]['nprots'] / res[folder]['n_clades']
    resdf = pd.DataFrame.from_dict(res, orient = 'index')
    dfs.append(resdf)
seqset_resdf = pd.concat(dfs)
seqset_resdf = pd.to_pickle('structure_and_sequence_qcmetrics.pkl')

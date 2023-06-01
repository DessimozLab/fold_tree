#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')


# In[2]:


import sys
sys.path.append('../')
from src import AFDB_tools, treescore , foldseek2tree
import pandas as pd
import torch


# In[3]:


#load cath dataset from CATH data folder with comment lines as # and white space seperator
domain_df = pd.read_csv('../CATH_data/cath-domain-list.txt' , comment='#' , header = None , delim_whitespace=True)
print(domain_df.head())


# In[4]:


#change the type to string
domain_df = domain_df.astype(str)
domain_df['superfam'] = domain_df[[1,2,3,4]].apply(lambda x : '.'.join(x), axis = 1) 
domain_df['CAT'] = domain_df[[1,2,3]].apply(lambda x : '.'.join(x), axis = 1) 


# In[5]:


siftsdf = pd.read_csv('../CATH_data/pdb_chain_cath_uniprot.csv', header=1 )


# In[6]:


#merge with domain df on cath id
siftsdf = siftsdf.merge(domain_df, left_on = 'CATH_ID' , right_on = 0 , how = 'left')
siftsdf = siftsdf.dropna()
print(siftsdf.head(), len(siftsdf) , print(len(siftsdf.superfam.unique())))


# In[7]:


from matplotlib import pyplot as plt
print(siftsdf.superfam.unique(),len(siftsdf.superfam.unique()))
#make a histogram of the superfam counts
famsizes = siftsdf.CAT.value_counts().to_numpy()
plt.hist( famsizes[famsizes>10] , bins = 100)
plt.show()


# In[8]:


import urllib.request
from Bio import PDB
import warnings
import os


def dlchain(pdb_id , chainID , savepath):
    
    try:
        if  os.path.isfile(savepath)==False:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pdbl = PDB.PDBList()
                pdb_file = pdbl.retrieve_pdb_file(pdb_id , file_format = 'pdb' )
                structure = PDB.PDBParser().get_structure(pdb_id, pdb_file )
                chain_A = None
                for model in structure:
                    for chain in model:
                        if chain.id == chainID:
                            chain_A = chain
                            break
                    if chain_A:
                        break

                if chain_A:
                    io = PDB.PDBIO()
                    io.set_structure(chain_A)
                    io.save(savepath)
                else:
                    print("err", savepath)
    except:
        print('err' , savepath)


# In[ ]:


from src import AFDB_tools
import os
import shutil
import random
import glob



#datapath = '../CATH_data/'

datapath = '../CAT_data/'

nprots = 250
fams = 2000
clear = False
total = 0
lengths = []
#iterate over all superfamilies and create a tree for each
if not os.path.exists(datapath):
    os.mkdir(datapath)

if clear == True:
    print('clearing folders')
    for folder in [datapath]:
        subfolders = glob.glob(folder+'*/')
        for sub in subfolders:
            shutil.rmtree(sub)

for i,superfam in enumerate(siftsdf.CAT.unique()):
    if not os.path.exists(datapath+superfam):
        os.mkdir(datapath+superfam)
    #create a tree for each superfam
    #sample 1000 proteins from the superfam
    sub = siftsdf[siftsdf.CAT == superfam]
    prots = list( set(sub['SP_PRIMARY'].unique()) )
    #output the sequence dataset to a file
    #create a folder for the superfam if it does not exist
    #output the uniport ids to a file
    if len(prots)> 100:
        random.shuffle(prots)
        prots = prots[:nprots]
        print(len(prots))
        
        
        lengths.append(len(prots))
        total+=1
        with open(datapath+superfam+'/identifiers.txt', 'w') as f:
            f.write('\n'.join( prots ))
            
        if not os.path.exists(datapath+superfam+'/structs/'):
            os.mkdir(datapath+superfam+'/structs/')
        protdict = dict(zip( sub['SP_PRIMARY'] , sub.PDB ))
        protdict = {p:protdict[p] for p in prots }
        chaindict = dict(zip( sub['SP_PRIMARY'] , sub.CHAIN ))
        [ dlchain(pdb, chaindict[unid] , datapath+superfam+'/structs/'+unid.upper()+'.pdb') for unid,pdb in protdict.items() ] 
print('done', total)
plt.hist(lengths)
plt.show()


# In[10]:


from src import AFDB_tools
import os
import shutil
import random
import glob

datapath = '../CATH_data/'
unfiltered = '../CATH_data_unfiltered/'

nprots = 500
fams = 2000
clear = True
total = 0
lengths = []
#iterate over all superfamilies and create a tree for each
if not os.path.exists(datapath):
    os.mkdir(datapath)
if not os.path.exists(unfiltered):
    os.mkdir(unfiltered)
    
if clear == True:
    print('clearing folders')
    for folder in [datapath,unfiltered]:
        subfolders = glob.glob(folder+'*/')
        for sub in subfolders:
            shutil.rmtree(sub)

for i,superfam in enumerate(siftsdf.superfam.unique()):
    if not os.path.exists(datapath+superfam):
        os.mkdir(datapath+superfam)
    if not os.path.exists(unfiltered+superfam):
        os.mkdir(unfiltered+superfam)
    #create a tree for each superfam
    #sample 1000 proteins from the superfam
    sub = siftsdf[siftsdf.superfam == superfam]
    prots = list( set(sub['SP_PRIMARY'].unique()) )
    #output the sequence dataset to a file
    #create a folder for the superfam if it does not exist
    #output the uniport ids to a file
    if len(prots)> 100:
        random.shuffle(prots)
        prots = prots[:nprots]
        print(len(prots))
        lengths.append(len(prots))
        total+=1
        with open(datapath+superfam+'/identifiers.txt', 'w') as f:
            f.write('\n'.join( prots ))
        with open(unfiltered+superfam+'/identifiers.txt', 'w') as f:
            f.write('\n'.join( prots ))
print('done', total)
plt.hist(lengths)
plt.show()


# In[ ]:


scopdf = pd.read_csv('../SCOP_data/pdb_chain_scop_uniprot.csv')
print(scopdf)


# In[ ]:


scopdf = pd.read_csv('../SCOP_data/scop-cla-latest.txt', header = None, comment= '#' , sep = ' ')
scopdf.columns = 'FA-DOMID FA-PDBID FA-PDBREG FA-UNIID FA-UNIREG SF-DOMID SF-PDBID SF-PDBREG SF-UNIID SF-UNIREG SCOPCLA'.split()
scopdf['SF'] = scopdf['SCOPCLA'].map( lambda x : ''.join(x.split(',')[0-2]) )
print(scopdf)

vc = scopdf.SF.value_counts()
plt.hist( vc[vc>10])


# In[ ]:


from src import AFDB_tools
import os
datapath = '../SCOP_data/'
nprots = 500
fams = 2000
total = 0
lengths = []
#iterate over all superfamilies and create a tree for each

if not os.path.exists(datapath):
    os.mkdir(datapath)
    
for i,superfam in enumerate(scopdf.SF.unique()):
    if not os.path.exists(datapath+superfam):
        os.mkdir(datapath+superfam)

    #create a tree for each superfam
    #sample 1000 proteins from the superfam
    sub = scopdf[scopdf.SF == superfam]
    prots = list(sub['SF-UNIID'].unique())
    prots = prots[:nprots]
    #output the sequence dataset to a file
    #create a folder for the superfam if it does not exist
    #output the uniport ids to a file
    if len(prots)> 100:
        lengths.append(len(prots))
        total+=1
        with open(datapath+superfam+'/identifiers.txt', 'w') as f:
            f.write('\n'.join( prots ))
print('done', total)
plt.hist(lengths)
plt.show()


# In[38]:


#subsample the oma folders

import os
import shutil
import random
import glob
from tqdm import tqdm
from matplotlib import pyplot as plt



#path to the oma folder
path = '../../../datasets/Structure_Trees_mk2/*/'
#path to the output folder
outpath = '../OMA_data/'

#make the output folder if it does not exist
if not os.path.exists(outpath):
    os.mkdir(outpath)

#number of proteins to sample
min_prots = 10
#number of superfamilies to sample
fams = 2000

superfams={}
#iterate over all superfamilies and compile an identifiers.txt file in the output folder
for clade in glob.glob(path):
    print(clade)
    if clade not in superfams:
        superfams[clade] = {}
    #stores length of each superfam
    superfam = set(glob.glob(clade+'*/'))-set([clade + 'logs/'])
    #use tqdm progress bar
    for s in tqdm(superfam):
        
        with open(clade+s+'identifiers.txt', 'r') as f:
                nprots = len(f.readlines())
        if nprots > min_prots:
            superfams[clade][s] = nprots
lengths = { c:len(superfams[c]) for c in superfams}
print(lengths)
plt.hist(list(lengths.values()))
plt.show()


# In[58]:


minlen = min(list(lengths.values()) )
print(minlen)
#clear output directory


#sample the superfamilies
for clade in superfams:
    cladestr = clade.split('/')[-2]+'/'
    if not os.path.isdir(outpath+cladestr):
        os.mkdir(outpath+cladestr)
    print(clade)
    superfams[clade] = random.sample( list(superfams[clade]) , minlen)
    for s in tqdm(superfams[clade]):
        #copy the files to the output folder
        superfam = s.split('/')[-2]+'/'
        if not os.path.exists(outpath+cladestr+superfam):
            os.mkdir(outpath+cladestr+superfam)
            shutil.copy(s+'identifiers.txt', outpath+cladestr+superfam+'identifiers.txt')


# In[ ]:





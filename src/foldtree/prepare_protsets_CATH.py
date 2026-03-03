sys.path.append('../../src')
from src import AFDB_tools, treescore , foldseek2tree
import urllib.request
from Bio import PDB
import warnings
import os
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import shutil
import random
import glob
import sys
import pandas as pd

#load cath dataset from CATH data folder with comment lines as # and white space seperator
domain_df = pd.read_csv('../../CATH_data/cath-domain-list.txt' , comment='#' , header = None , delim_whitespace=True)
print(domain_df.head())
#change the type to string
domain_df = domain_df.astype(str)
domain_df['superfam'] = domain_df[[1,2,3,4]].apply(lambda x : '.'.join(x), axis = 1) 
domain_df['CAT'] = domain_df[[1,2,3]].apply(lambda x : '.'.join(x), axis = 1)
siftsdf = pd.read_csv('../../CATH_data/pdb_chain_cath_uniprot.csv', header=1 )

#merge with domain df on cath id
siftsdf = siftsdf.merge(domain_df, left_on = 'CATH_ID' , right_on = 0 , how = 'left')
siftsdf = siftsdf.dropna()
print(siftsdf.head(), len(siftsdf) , print(len(siftsdf.superfam.unique())))
print(siftsdf.superfam.unique(),len(siftsdf.superfam.unique()))


nprots = 250
fams = 2000
clear = False
total = 0
lengths = []

def dlchain(pdb_id , chainID , savepath):
    try:
        if  os.path.isfile(savepath)==False:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pdbl = PDB.PDBList()
                pdb_file = pdbl.retrieve_pdb_file(pdb_id , file_format = 'pdb' )
                structure = PDB.PDBParser().get_structure(pdb_id, savepath)
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
                    fixer = PDBFixer(filename=savepath)
                    fixer.findNonstandardResidues()
                    fixer.replaceNonstandardResidues()
                    fixer.removeHeterogens(True)
                    fixer.findMissingResidues()
                    fixer.findMissingAtoms()
                    fixer.addMissingAtoms()
                    fixer.addMissingHydrogens(7.0)
                    PDBFile.writeFile(fixer.topology, fixer.positions, open(savepath, 'w'))

                else:
                    print("err", savepath)
    except:
        print('err' , savepath)


datapaths = {'../../CAT_data/':'CAT' , '../../CATH_data/' : 'superfam' }

#iterate over all superfamilies and create a tree for each
for datapath,category in datapaths.items():
    if not os.path.exists(datapath):
        os.mkdir(datapath)
    if clear == True:
        print('clearing folders')
        for folder in [datapath]:
            subfolders = glob.glob(folder+'*/')
            for sub in subfolders:
                shutil.rmtree(sub)

    for i,fam in enumerate(siftsdf[category].unique()):
        if not os.path.exists(datapath+fam):
            os.mkdir(datapath+fam)
        #create a tree for each superfam
        #sample 1000 proteins from the superfam
        sub = siftsdf[siftsdf[category] == fam]
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
            with open(datapath+fam+'/identifiers.txt', 'w') as f:
                f.write('\n'.join( prots ))
            if not os.path.exists(datapath+fam+'/structs/'):
                os.mkdir(datapath+fam+'/structs/')
            protdict = dict(zip( sub['SP_PRIMARY'] , sub.PDB ))
            protdict = {p:protdict[p] for p in prots }
            chaindict = dict(zip( sub['SP_PRIMARY'] , sub.CHAIN ))
            [ dlchain(pdb, chaindict[unid] , datapath+fam+'/structs/'+unid.upper()+'.pdb') for unid,pdb in protdict.items() ] 
print('done', total)

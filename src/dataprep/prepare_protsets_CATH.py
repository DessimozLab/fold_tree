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
import time
import tqdm
import pandas as pd
import multiprocessing as mp
sys.path.append('../')
import AFDB_tools, treescore , foldseek2tree
#argument parser for the script
import argparse
import errno
import warnings
from concurrent.futures import TimeoutError
from pebble import ProcessPool, ProcessExpired


import os


if __name__ == '__main__':
    nprots = 50
    fams = 2000
    clear = False
    total = 0
    parser = argparse.ArgumentParser(description='Prepare CAT and CATH datasets.')
    parser.add_argument('--nprots', type=int, default=250, help='number of proteins to download')
    parser.add_argument('--fams', type=int, default=2000, help='number of families to download')
    parser.add_argument('--clear', type=bool, default=False, help='clear the data folder')
    args = parser.parse_args()

    #if the arguments are present then set the vars
    #otherwise use the default values
    if args.nprots:
        nprots = args.nprots
    if  args.fams:
        fams = args.fams
    if args.clear:
        clear = args.clear

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
    print(siftsdf.head() , len(siftsdf))
    lengths = []

    
    def dlchain(pdb_id , chainID , savepath , unid , verbose = False):
        if  os.path.isfile(savepath) == False:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pdbl = PDB.PDBList(verbose = verbose )
                pdb_file = pdbl.retrieve_pdb_file(pdb_id , file_format = 'pdb')
                structure = PDB.PDBParser().get_structure(pdb_id, pdb_file)
                if len(structure)>1:
                    for model in structure:
                            for chain in model:
                                if chain.id == chainID:
                                    chain_A = chain
                                    break
                else:
                    chain_A = structure[0][chainID]
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
                    #fixer.addMissingHydrogens(7.0)
                    PDBFile.writeFile(fixer.topology, fixer.positions, open(savepath, 'w'))
                    return None
                else:
                    print(model)
                    print("err", savepath)
                    return None


    
    datapaths = { '../../CATH_data_mk2/' : 'superfam' ,  '../../CAT_data_mk2/':'CAT'  }

    #iterate over all superfamilies and create a tree for each
    for datapath,category in datapaths.items():
        if not os.path.exists(datapath):
            os.mkdir(datapath)    
        if clear == True:
            print('clearing folders')
            subfolders = glob.glob(datapath+'*/')
            for sub in subfolders:
                shutil.rmtree(sub)
        
        for fam in tqdm.tqdm(siftsdf[category].unique()):
            if not os.path.exists(datapath+fam):
                os.mkdir(datapath+fam)
            


            if not os.path.exists(datapath+fam+ datapath+fam+'/sequences.fst'):
                #create a tree for each superfam
                sub = siftsdf[siftsdf[category] == fam]
                
                
                prots = list( set(sub['SP_PRIMARY'].unique()) )
                #output the uniport ids to a file
                if len(prots)> 10:
                    print(fam)
                    if not os.path.exists(datapath+fam + 'identifiers.txt'):
                        random.shuffle(prots)
                        prots = prots[:nprots]
                        with open(datapath+fam + 'identifiers.txt', 'w') as f:
                            f.write('\n'.join( prots ))
                    else:
                        with open(datapath+fam + 'identifiers.txt') as f:
                            prots = [ i.strip() for i in f if len(i.strip())>0 ]
                    
                    total+=1
                    protdict = dict(zip( sub['SP_PRIMARY'] , sub.PDB ))
                    protdict = {p:protdict[p] for p in prots }
                    chaindict = dict(zip( sub['SP_PRIMARY'] , sub.CHAIN ))

                    if not os.path.exists(datapath+fam+'/structs/'):
                        os.mkdir(datapath+fam+'/structs/')
                    with ProcessPool() as pool:
                        futures = [ pool.schedule( dlchain,  ( protdict[unid], chaindict[unid] , datapath+fam+'/structs/'+unid.upper()+'.pdb' , unid.upper() , False ) , timeout = 60) for unid in prots  ] 
                        for future in tqdm.tqdm(futures, total=len(prots)):
                            try:
                                results = future.result()
                            except TimeoutError as error:
                                print("unstable_function took longer than %d seconds" % error.args[1])
                            except ProcessExpired as error:
                                print("%s. Exit code: %d" % (error, error.exitcode))
                            except Exception as error:
                                print("unstable_function raised %s" % error)
                                print(error.traceback)  # Python's traceback of remote process

                        pool.close()
                        pool.join()
                        
                    with open(datapath+fam+'/identifiers.txt', 'w') as f:
                        f.write('\n'.join( prots ))
                 

    print('done', total)

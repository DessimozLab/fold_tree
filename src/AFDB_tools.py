
import os 
import wget 
from io import StringIO
import pandas as pd
import numpy as np
import requests

def grab_struct(uniID, structfolder, overwrite=False):
    try:
        os.mkdir(structfolder)
    except:
        pass
    print(uniID)
    try:
        prefix = 'https://alphafold.ebi.ac.uk/files/AF-'
        post = '-F1-model_v3.pdb'
        url = prefix+uniID.upper()+post
        if not os.path.isfile(structfolder + uniID +'.pdb'):
            wget.download(url, structfolder + uniID +'.pdb')
    except:
        print('structure not found', uniID)
        return uniID
    return None

#lets pull in the tor pathway members with some gene names
def unirequest_tab(name, verbose = False):
    #we query first by protein name and then gene name
    url = 'http://rest.uniprot.org/uniprotkb/stream?'
    params = [
    #
    'query={}'.format(name),
    'fields=id,gene_names,protein_name,reviewed,protein_name,organism_name,lineage_ids,sequence',
    'format=tsv',
    ]
    params = ''.join([ p+'&' for p in params ])[:-1]
    data = requests.get(url+params).text
    #only return the first hit for each query    
    try:
        data =  pd.read_table(StringIO(data)).iloc[0]
        data['query'] = name
        if verbose == True:
            print(data)

        return data    
    except:
        print('err uniprot API' , name)

def grab_entries(ids, verbose = False):
    name_results = pd.concat([unirequest_tab(name) for name in ids] , axis = 1 , ignore_index= True).transpose()
    if verbose == True:
        print(name_results)
    return name_results

def res2fasta(unires_df):
    unires_df['fasta'] = unires_df[ ['query' , 'Sequence']].apply( lambda r : '> '+ r.query + '\n'+ r.Sequence+ '\n' , axis = 1)
    fasta = ''.join(unires_df.fasta)
    return fasta

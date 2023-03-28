#!/usr/bin/env python
# coding: utf-8

# # Overview
# Given a PDB file, performe a remote search against the AFDB50 using FoldSeek

# # Load libraries
from requests import get, post
import sys, time, glob, pandas as pd, tqdm, os


# # Define auxiliary functions
def foldseek_search(pdb_path, output_file):
    try:
        # load PDB
        time.sleep(1)
        if not os.path.exists(output_file):
            with open(pdb_path, 'r') as file:
                data = file.read()
            # starting search
            ticket = post('https://search.foldseek.com/api/ticket', {
                            'q' : data,
                            'database[]' : ['afdb50', 'afdb-swissprot', 'pdb100', 'afdb-proteome', 'mgnify_esm30', 'gmgcl_id'],
                            'mode' : '3diaa',
                        }).json()
            # collecting results
            repeat = True
            while repeat:
                status = get('https://search.foldseek.com/api/ticket/' + ticket['id']).json()
                if status['status'] == "ERROR":
                    # handle error
                    sys.exit(0)
            
                # wait a short time between poll requests
                time.sleep(1)
                repeat = status['status'] != "COMPLETE"
            # get all hits for the first query (0)
            result = get('https://search.foldseek.com/api/result/' + ticket['id'] + '/0').json()
            # iterating over results
            foldseek_results_rows = []
            # voy por base de datos
            for results_dbs in result['results']:
                # recorro los alineamientos
                for alignment in results_dbs['alignments']:
                    # agrego la fila
                    foldseek_results_rows.append(pd.DataFrame.from_dict({'pdb_file': [pdb_path],
                                                                         'db': [results_dbs['db']], 
                                                                         'query': [alignment.get('query', '-')],
                                                                         'target': [alignment.get('target', '-')],
                                                                         'seqId': [alignment.get('seqId', '-')],
                                                                         'alnLength': [alignment.get('alnLength', '-')],
                                                                         'missmatches': [alignment.get('missmatches', '-')],
                                                                         'gapsopened': [alignment.get('gapsopened', '-')],
                                                                         'qStartPos': [alignment.get('qStartPos', '-')],
                                                                         'qEndPos': [alignment.get('qEndPos', '-')],
                                                                         'dbStartPos': [alignment.get('dbStartPos', '-')],
                                                                         'dbEndPos': [alignment.get('dbEndPos', '-')],
                                                                         'prob': [alignment.get('prob', '-')],
                                                                         'eval': [alignment.get('eval', '-')],
                                                                         'score': [alignment.get('score', '-')],
                                                                         'qLen': [alignment.get('qLen', '-')],
                                                                         'dbLen': [alignment.get('dbLen', '-')],
                                                                         'qAln': [alignment.get('qAln', '-')],
                                                                         'dbAln': [alignment.get('dbAln', '-')],
                                                                         'tCa': [alignment.get('tCa', '-')],
                                                                         'tSeq': [alignment.get('tSeq', '-')],
                                                                         'taxId': [alignment.get('taxId', '-')],
                                                                         'taxName': [alignment.get('taxName', '-')]}))
            # return result
            return pd.concat(foldseek_results_rows)
        else:
            pass
    except Exception as e:
        print(e)

def create_dir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)


# # Perform remote search on seed structures
for input_pdb in tqdm.tqdm(snakemake.input):
    # get the PDB ID
    pdb_id = input_pdb.rpartition('/')[2].rpartition('.')[0]
    # create directories to allocate results
    create_dir('{0}/foldseek_searches'.format(snakemake.wildcards.seed_folder)) # modif para agregar el tema del {folder}
    # get the output name
    output_file = '{0}/foldseek_searches/foldseek_search_{1}.tsv'.format(snakemake.wildcards.seed_folder, pdb_id) # modif para agregar el tema del {folder}
    # perform search
    result_table = foldseek_search(pdb_path = input_pdb, output_file = output_file)
    # save TSV file
    result_table.to_csv(output_file, sep = '\t', index = False)


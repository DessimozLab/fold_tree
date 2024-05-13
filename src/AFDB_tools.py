
import os 
import wget 
from io import StringIO
import pandas as pd
import requests
from scipy.stats import describe
from Bio.PDB import *
from Bio.SeqUtils import seq1

import time
import numpy as np



def get_amino_acid_sequence(pdb_filename):
	""" This function extracts the amino acid sequence from a PDB file."""
	# Initialize the parser
	parser = PDBParser()
	# Parse the structure from file
	structure = parser.get_structure('PDB_structure', pdb_filename)
	# Extract amino acid sequences from the structure

	for model in structure:
		for chain in model:
			sequence = ''
			for residue in chain.get_residues():
				if residue.get_resname() in ['HOH', 'WAT']:
					continue  # Ignore water molecules
				sequence += seq1(residue.get_resname(), undef_code='X')
			return sequence



def descr(pdb_path):
	'''
	Extracts the plddt (in the beta factor column) of the first atom of each residue in a PDB file and returns a descriptive statistics object.
	Parameters:
		pdb_path (str): The path to the PDB file.'''

	lppd=[]
	parser = PDBParser()
	struc = parser.get_structure("a", pdb_path)
	for res in struc.get_residues():
		for at in res.get_atoms():
		   lppd.append(at.get_bfactor())
		   break
	return describe(lppd)


def filter_plddt(pdb_path , thresh= .6 , minthresh = .5 ):
	'''
	Extracts the plddt (in the beta factor column) of the first atom of each residue in a PDB file and returns bool if the pdb is accepted or not.

	Parameters:
		pdb_path (str): The path to the PDB file.'''
	
	lddt=[]
	parser = PDBParser()
	struc = parser.get_structure("a", pdb_path)
	for res in struc.get_residues():
		for at in res.get_atoms():
		   lddt.append(at.get_bfactor())
		   break
	if np.mean(lddt) < thresh or np.amin(lddt) < minthresh:
		return False
	else:
		return True


def grab_struct(uniID, structfolder, rejected = None, overwrite=False):

	"""
	Downloads a protein structure file from the AlphaFold website and saves it to the specified folder.
	
	Parameters:
	uniID (str): The UniProt ID of the protein for which the structure is being downloaded.
	structfolder (str): The path to the folder where the structure file should be saved.
	overwrite (bool, optional): A flag indicating whether to overwrite an existing file with the same name in the specified folder. Defaults to False.
	
	Returns:
	None: If the file is successfully downloaded or if overwrite is set to True and a file with the same name is found in the specified folder.
	str: If an error occurs during the download or if a file with the same name is found in the specified folder and overwrite is set to False.
	
	Examples:
	>>> grab_struct('P00533', '/path/to/structures/')
	None
	>>> grab_struct('P00533', '/path/to/structures/', overwrite=True)
	None
	"""

	try:
		os.mkdir(structfolder)
	except:
		pass
	try:
		prefix = 'https://alphafold.ebi.ac.uk/files/AF-'
		post = '-F1-model_v4.pdb'
		url = prefix+uniID.upper()+post
		if not os.path.isfile(structfolder + uniID +'.pdb'):
			if rejected is None or (rejected and not os.path.isfile(structfolder + uniID +'.pdb')):
				wget.download(url, structfolder + uniID +'.pdb')
	except:
		print('structure not found', uniID)
		return uniID
	return None



#lets pull in the tor pathway members with some gene names

def chunk(data,csize):
	return [data[x:x+csize] for x in range(0, len(data), csize)]

def unirequest_tab(name, verbose = False):

	"""
	Makes a request to the UniProt API and returns information about a protein in tab-separated format.
	
	Parameters:
	name (str): The name of the protein for which information is being requested.
	verbose (bool, optional): A flag indicating whether to print the returned data to the console. Defaults to False.
	
	Returns:
	pd.DataFrame: A DataFrame containing information about the protein, with one row for each hit in the search.
	
	Examples:
	>>> unirequest_tab('P00533')
															 id  ...                                            sequence
	0  sp|P00533|1A2K_HUMAN RecName: Full=Alpha-2-...  ...  MPTSVLLLALLLAPAALVHVCRSRFPKCVVLVNVTGLFGN...
	"""
	#we query first by protein name and then gene name
	url = 'http://rest.uniprot.org/uniprotkb/stream?'
	params = [
	'query=accession:{}'.format(name),
	'fields=id,accession,gene_names,protein_name,reviewed,protein_name,organism_name,lineage_ids,sequence',
	'format=tsv',
	]
	params = ''.join([ p+'&' for p in params ])[:-1]
	data = requests.get(url+params).text
	#only return the first hit for each query    
	try:
		data =  pd.read_table(StringIO(data))
		data['query'] = data['Entry']
		data = data[ data['Entry'].isin(name.split('+OR+'))]
		if verbose is True:
			print(data)
		return data    
	except:
		print('error', data )
		time.sleep(10)
		unirequest_tab(name, verbose = True)

def grab_entries(ids, verbose = False):
	"""
	Makes requests to the UniProt API for information about proteins with the given IDs.
	
	Parameters:
	ids (list): A list of UniProt IDs for the proteins for which information is being requested.
	verbose (bool, optional): A flag indicating whether to print the returned data to the console. Defaults to False.
	
	Returns:
	pd.DataFrame: A DataFrame containing information about the proteins, with one row for each hit in the search.
	
	Examples:
	>>> grab_entries(['P00533', 'P15056'])
															 id  ...                                            sequence
	0  sp|P00533|1A2K_HUMAN RecName: Full=Alpha-2-...  ...  MPTSVLLLALLLAPAALVHVCRSRFPKCVVLVNVTGLFGN...
	1  sp|P15056|1A01_HUMAN RecName: Full=Alpha-1-...  ...  MAAARLLPLLPLLLALALALTETSCPPASQGQRASVGDRV...
	
	Notes:
	This function makes requests to the UniProt API for information about proteins with the given IDs. If a request is successful, the returned data is processed and added to a DataFrame. If a request is unsuccessful, an error message is printed to the console.
	"""
	try:
		name_results = pd.concat([unirequest_tab( '+OR+'.join(c) , verbose = verbose) for c in chunk(ids, 50 )] , ignore_index= True)
	except:
		print('error', ids )
		time.sleep(10)
		name_results = pd.concat([unirequest_tab( '+OR+'.join(c) , verbose = verbose) for c in chunk(ids, 50 )] , ignore_index= True)
		
	if verbose == True:
		print(name_results)
	return name_results

def res2fasta(unires_df):
	"""
	Converts a DataFrame containing protein information into a FASTA format string.
	
	Parameters:
	unires_df (pd.DataFrame): A DataFrame containing information about proteins, with columns 'query' and 'Sequence' representing the name and sequence of each protein, respectively.
	
	Returns:
	str: A string in FASTA format representing the proteins in the input DataFrame.
	
	Examples:
	>>> unires_df = pd.DataFrame([{'query': 'P00533', 'Sequence': 'MPTSVLLLALLLAPAALVHVCRSRFPKCVVLVNVTGLFGN'}])
	>>> res2fasta(unires_df)
	'> P00533\nMPTSVLLLALLLAPAALVHVCRSRFPKCVVLVNVTGLFGN\n'
	"""
	
	unires_df = unires_df.drop_duplicates(subset=['query'])
	unires_df['query'] = unires_df['query'].map( lambda x : '>' + x + '\n')
	unires_df['fasta'] = unires_df['query'] + unires_df.Sequence
	unires_df['fasta'] = unires_df['fasta'].map( lambda x : x + '\n')
	fasta = ''.join(unires_df.fasta)
	return fasta

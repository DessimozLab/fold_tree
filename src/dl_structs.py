import AFDB_tools
import os
import shutil
import glob
import pandas as pd
from Bio import PDB as pdb

'''
This script is used to download the structures from the afdb
it also filters the structures based on the plddt score

'''



infolder = snakemake.input[0].split('/')[:-1]
infolder = ''.join( [i + '/' for i in infolder])
structfolder = infolder+'structs/'
rejectedfolder = infolder+'rejected/'

#remove tmp folder

try:
	shutil.rmtree(infolder+'tmp/')
except:
	pass

custom_structs = snakemake.params.custom_structs

if custom_structs == True:
	print('custom structures, skipping download of structures')
	found = glob.glob(structfolder+'*.pdb')
	finalset = { f.replace('.pdb', '' ).split('/')[-1] : AFDB_tools.get_amino_acid_sequence(f) for f in found }
	with open(snakemake.output[0] , 'w') as outfile:
		outfile.write(''.join(['>'+i+'\n'+finalset[i]+'\n' for i in finalset]))
	with open(snakemake.output[1] , 'w') as outfile:
		outfile.write(''.join(['>'+i+'\n'+finalset[i]+'\n' for i in finalset]))

else:
	try:
		os.mkdir(structfolder)
	except:
		print(structfolder , 'already exists ')

	try:
		os.mkdir(rejectedfolder)
	except:
		print(rejectedfolder , 'already exists ')



	#with open(snakemake.input[0]) as infile:
	#	ids = [ i.strip() for i in infile if len(i.strip())>0 ]
	seqdf = pd.read_csv(snakemake.input[0])
	ids = list(seqdf['query'].unique())

	missing = [	AFDB_tools.grab_struct(i, structfolder, rejectedfolder) for i in ids ]
	found = glob.glob(structfolder+'*.pdb') + glob.glob(rejectedfolder+'*.pdb')
	found = { i.split('/')[-1].replace('.pdb',''):i for i in found}
	missing_structs = set(ids)-set(found.keys())
	filtervar = snakemake.params.filtervar

	filtervar_min = snakemake.params.filtervar_min
	filtervar_avg = snakemake.params.filtervar_avg


	#get plddt from afdb structures and remove those with avg plddt < 0.4
	if filtervar == True:
		plddt = { i:AFDB_tools.filter_plddt( found[i] , thresh= filtervar_avg , minthresh = filtervar_min ) for i in found}
	else:
		plddt = { i:True for i in found}

	for i in list(found.keys()):
		if i not in ids or plddt[i] is False:
			#move to rejected folder
			if not os.path.isfile(rejectedfolder + i + '.pdb'):
				shutil.move(found[i], rejectedfolder)
			else:
				os.remove(found[i])
			missing_structs.add(i)
			del found[i]
	
	
	missing_sequences = set(ids)-set(seqdf['query'].unique())
	print('missing in afdb:',missing_structs)
	print( 'missing in sequences:',missing_sequences)
	finalset = set(ids)-set(missing_sequences)
	finalset = set(finalset)-set(missing_structs)

	resdf = seqdf[seqdf['query'].isin(finalset)]
	#fasta = AFDB_tools.res2fasta(resdf)

	found = glob.glob(structfolder+'*.pdb') + glob.glob(rejectedfolder+'*.pdb')
	finalset = { f.replace('.pdb', '' ).split('/')[-1] : AFDB_tools.get_amino_acid_sequence(f) for f in found if f.replace('.pdb', '' ).split('/')[-1] in finalset }
	
	assert len(finalset) == len(resdf['query'].unique()) , 'finalset and resdf do not have the same length'

	#with open(snakemake.output[0] , 'w') as outfile:
	#	outfile.write(''.join(['>'+i+'\n'+finalset[i]+'\n' for i in finalset]))
	resdf.to_csv(snakemake.output[0])
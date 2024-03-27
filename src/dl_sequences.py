import AFDB_tools
import os
import shutil	

custom_structs = snakemake.params.custom_structs

print('custom_structs: ', custom_structs)

if custom_structs == True and snakemake.params.cath == False:
	print('custom structures, skipping download of sequences')
	#writes a dummy file with the structs included
	with open(snakemake.output[0] , 'w') as outfile:
		outfile.write('')

else:
	with open(snakemake.input[0]) as infile:
		ids = [ i.strip() for i in infile if len(i.strip())>0 ]

	resdf = AFDB_tools.grab_entries(ids, verbose = True)
	resdf.to_csv(snakemake.output[0])

if snakemake.params.clean == True:
	basedir = snakemake.params.basedir
	#move all structs from rejected to structs
	for i in ids:
		shutil.move(basedir + 'rejected/*.pdb', basedir + 'structs/*.pdb')
	#delete the tmp dir
	shutil.rmtree(basedir + 'tmp')
	#delete all files except identifiers.txt in basedir
	for file in os.listdir(basedir):
		if file != 'identifiers.txt':
			os.remove(basedir + file)
	
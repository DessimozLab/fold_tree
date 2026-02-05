import AFDB_tools
import os
import shutil	

custom_structs = snakemake.params.custom_structs
print('custom_structs: ', custom_structs)

clean = snakemake.params.clean_folder
print( 'clean: ', clean	)

basedir = snakemake.input[0].split('/')[:-1]
basedir = ''.join( [i + '/' for i in basedir])

if  clean == True:
	#move all structs from rejected to structs
	pass
	"""
	if os.path.exists(basedir + 'structs/') == False:
		os.mkdir(basedir + 'structs/')
	if os.path.exists(basedir + 'rejected/') == False:
		os.mkdir(basedir + 'rejected/')
	else:
		ids = os.listdir(basedir + 'rejected/')
		for i in ids:
			shutil.move(basedir + 'rejected/'+i, basedir + 'structs/' +i)
	#delete all files except identifiers.txt in basedir
	for file in os.listdir(basedir):
		if file != 'identifiers.txt' and os.path.isfile(basedir + file):
			os.remove(basedir + file)
	"""

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

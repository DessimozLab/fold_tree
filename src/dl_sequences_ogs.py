import AFDB_tools
import os



custom_structs = snakemake.params.custom_structs

print('custom_structs: ', custom_structs)

if custom_structs == True:
	print('custom structures, skipping download of sequences')
	#writes a dummy file with the structs included
	with open(snakemake.output[0] , 'w') as outfile:
		outfile.write('')




else:
	for i,file in enumerate(snakemake.input):
		with open(file) as infile:
			ids = [ i.strip() for i in infile if len(i.strip())>0 ]
		#make output directory for structures if it doesn't exist
		outdir = os.path.dirname( '/'.join(file.split('/')[:-1]) + file.split('_')[0] )
		if not os.path.exists(outdir):	
			os.makedirs(outdir)
		resdf = AFDB_tools.grab_entries(ids, verbose = True)
		resdf.to_csv(snakemake.output[i])
		
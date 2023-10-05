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
		with open(snakemake.input[0]) as infile:
			ids = [ i.strip() for i in infile if len(i.strip())>0 ]
		#make output directory if it doesn't exist
		outdir = os.path.dirname(snakemake.output[i].split('_')[0])
		if not os.path.isdir(outdir):
			os.mkdir(outdir)
		resdf = AFDB_tools.grab_entries(ids, verbose = True)
		resdf.to_csv(snakemake.output[0])

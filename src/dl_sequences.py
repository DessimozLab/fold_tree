import AFDB_tools

with open(snakemake.input[0]) as infile:
	ids = [ i.strip() for i in infile if len(i.strip())>0 ]
print(ids)

resdf = AFDB_tools.grab_entries(ids, verbose = True)


fasta = AFDB_tools.res2fasta(resdf)

resdf.to_csv(snakemake.output[1])

with open(snakemake.output[0] , 'w') as outfile:
	outfile.write(fasta)

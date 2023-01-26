import AFDB_tools

with open(snakemake.input[0]) as infile:
	ids = [ i.strip() for i in infile if len(i.strip())>0 ]
print(ids)

resdf = AFDB_tools.grab_entries(ids, verbose = True)
resdf.to_csv(snakemake.output[0])

import AFDB_tools
import os
import json
import pandas as pd

infolder = snakemake.input[0].split('/')[:-1]
infolder = ''.join( [i + '/' for i in infolder])
structfolder = infolder+'structs/'
try:
	os.mkdir(structfolder)
except:
	print(structfolder , 'already exists ')
#with open(snakemake.input[0]) as infile:
#	ids = [ i.strip() for i in infile if len(i.strip())>0 ]
seqdf = pd.read_csv(snakemake.input[0])
ids = list(seqdf.query)
print(ids)
missing = [	AFDB_tools.grab_struct(i, structfolder) for i in ids ]
missing = [ i for i in missing if i]
print('missing in afdb:',missing)
finalset = set(ids)-set(missing)
resdf = seqdf[seqdf.query.isin(finalset)]
fasta = AFDB_tools.res2fasta(resdf)
with open(snakemake.output[0] , 'w') as outfile:
	outfile.write(fasta)

with open(snakemake.output[1] , 'w') as outfile:
	outfile.write(json.dumps(finalset))
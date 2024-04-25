import pandas as pd
from scipy.stats import describe
import json

alnstats = {}
res = pd.read_table(snakemake.input[0], header = None)
print(res.head())

res[0] = res[0].map(lambda x :x.replace('.pdb', ''))
res[1] = res[1].map(lambda x :x.replace('.pdb', ''))

if snakemake.params.fmt is None:
    res.columns = 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,lddtfull,alntmscore'.split(',')
else:
    res.columns = snakemake.params.fmt.split(',')

ids = list( set(list(res['query'].unique()) + list(res['target'].unique())))
pos = { protid : i for i,protid in enumerate(ids)}
k = ['fident']
res[k] = res[k].astype(float)
#change nan to 0
res = res.fillna(0)
#exclude self alignments
res = res[res['query'] != res['target']]

#output a the alignment stats in json
alnstats[snakemake.input[0]] = describe(list(res['fident'].values))
print(snakemake.input[0], describe(list(res['fident'].values)))

with open(snakemake.output[0] , 'w') as out:
    json.dump(alnstats, out)


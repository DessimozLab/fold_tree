import ujson as json
import AFDB_tools
import glob

path = ''.join( [ p + '/' for p in snakemake.input[0].split('/')[:-1]]) +'structs/*.pdb'
all_structs ={}
for struct in glob.glob(path):
    all_structs[struct] = AFDB_tools.descr(struct)
    #unpack scipy stats describe object into a dictionary
    all_structs[struct] = { k:v for k,v in zip( ['nobs', 'minmax', 'mean', 'variance', 'skewness', 'kurtosis'], [ float(x) if type(x) is not tuple else (float(x[0]),float(x[1])) for x in list(all_structs[struct])]  ) }
#write to log file
with open(snakemake.log[0], 'w' ) as lfout:
    lfout.write(str(all_structs))

with open(snakemake.output[0], 'w' ) as rfout:
    rfout.write(json.dumps(all_structs))
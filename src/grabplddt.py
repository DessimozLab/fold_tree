import ujson as json
import AFDB_tools
import glob

path = ''.join( [ p + '/' for p in snakemake.input[0].split('/')[:-1]]) +'structs/*.pdb'
all_structs ={}
for struct in glob.glob(path):
    all_structs[struct] = AFDB_tools.descr(struct)
print(all_structs)
with open(snakemake.output[0], 'w' ) as rfout:
    rfout.write(json.dumps(all_structs))
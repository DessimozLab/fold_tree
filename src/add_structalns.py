
import toytree
import structalns
import os
import pandas as pd
import glob


print(snakemake.input)
print(snakemake.output)

alndf = pd.read_table(snakemake.input[0], header = None)
tre = toytree.tree(snakemake.input[1] ) 

infolder = snakemake.input[0].split('/')[:-1]
infolder = ''.join( [i + '/' for i in infolder])

mapper3di, mapperAA = structalns.read_dbfiles3di( snakemake.input[2] , snakemake.input[3])

#add the 3di alignment to the dataframe
columns = 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,qaln,taln,cigar,alntmscore'.split(',')
alndf.columns = columns

print('submat', snakemake.params.submat)


alndf['query'] = alndf['query'].map(lambda x :x.replace('.pdb', ''))
alndf['target'] = alndf['target'].map(lambda x :x.replace('.pdb', ''))

alndf['3diq']= alndf['query'].map(mapper3di)
alndf['3dit']= alndf['target'].map(mapper3di)
alndf['AAq']= alndf['query'].map(mapperAA)
alndf['AAt']= alndf['target'].map(mapperAA)

#output a fasta with the 3di sequences
res = alndf.apply(structalns.calc_fident_crossaln , axis = 1)
alndf = pd.concat([alndf,res] , axis = 1)

with open(snakemake.output[0] , 'w') as out:
    for seq in alndf['query'].unique():
        out.write('>'+seq.replace('.pdb', '' )+'\n')
        out.write(mapperAA[seq.replace('.pdb', '')]+'\n')

with open(snakemake.output[1] , 'w') as out:
    for seq in alndf['query'].unique():
        out.write('>'+seq.replace('.pdb', '' )+'\n')
        out.write(mapper3di[seq.replace('.pdb', '')]+'\n')

#prepare tree attributes
for i,n in enumerate(tre.treenode.traverse()):
    n.aln = None
    n.aln3di = None
    n.leafset = None
    if len(n.name) == 0:
        n.name = 'internal_'+str(i)


"""
alnfolder = infolder+'alnscratch/'
if not os.path.exists(alnfolder):
    os.mkdir(infolder+'alnscratch/')




finalaln, finalaln3di = structalns.traverse_tree_merge_mafft( tre.treenode , structalns.get_leafset(tre.treenode) , alndf , infolder+'alnscratch/' , submat = snakemake.params.submat , verbose = True ) 

print('finalaln',finalaln)
#print the final alignments
print('nsequences' , len(tre.get_tip_labels()))


finalaln = structalns.remove_seeds(finalaln)
finalaln3di = structalns.remove_seeds(finalaln3di)

finalaln = structalns.remove_redundant(finalaln)
finalaln3di = structalns.remove_redundant(finalaln3di)

with open(finalaln) as f:
    aacount = f.read().count('>')

with open(finalaln3di) as f:
    count3di = f.read().count('>')

assert len(tre.get_tip_labels()) == aacount
assert len(tre.get_tip_labels()) == count3di


for fasta,data in {snakemake.output[2]:finalaln, snakemake.output[3]:finalaln3di}.items():
    with open(fasta , 'w') as fastout:
        with open( data ) as fastin:
            fastout.write( fastin.read())

#cleanup the aln files
for f in glob.glob(infolder+'alnscratch/*inter*'):
    os.remove(f)

"""
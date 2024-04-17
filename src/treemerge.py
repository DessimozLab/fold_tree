import os
import glob
import toytree
import structalns
import pandas as pd

tre = toytree.tree(snakemake.input[0] ) 
alndf = pd.read_table(snakemake.input[1], header = None)

#prepare tree attributes
for i,n in enumerate(tre.treenode.traverse()):
    n.aln = None
    n.aln3di = None
    n.leafset = None
    if len(n.name) == 0:
        n.name = 'internal_'+str(i)
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

for fasta,data in {snakemake.output[0]:finalaln, snakemake.output[1]:finalaln3di}.items():
    with open(fasta , 'w') as fastout:
        with open( data ) as fastin:
            fastout.write( fastin.read())
#cleanup the aln files
for f in glob.glob(infolder+'alnscratch/*inter*'):
    os.remove(f)

import toytree
import pandas as pd

#open the tree
t = toytree.tree(snakemake.input[0])
#open the labels
labels = pd.read_csv(snakemake.input[1])
print('swapping' , labels.head())
#read the species names
labels['Organism'] = labels['Organism'].map( lambda x : x.replace(' ', '_').replace('(', '').replace(')', '').replace('.', '_') )
mapper = dict(zip(labels['Entry'], labels['Organism']))
print(mapper)
print(t)

#change leaf names to organism names
for leaf in t.treenode.traverse():
    if leaf.is_leaf():
        leaf.name = mapper[leaf.name]
#write the tree
print(t)
with open( snakemake.output[0] , 'w') as outfile:
    outfile.write(t.write())
    #add newlines to the end of the file
    outfile.write('\n')
    
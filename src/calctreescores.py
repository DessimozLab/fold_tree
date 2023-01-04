import json
import treescore
import pandas as pd
import toytree


#calc the taxscore 
uniprot_df = pd.readcsv(snakemake.input[0])
tree = toytree.tree(snakemake.input[1])
lineages = treescore.make_lineages(uniprot_df)
tree = treescore.label_leaves( tree , leaf_lineages)
overlap = treescore.getTaxOverlap(tree.treenode)
taxscore = tree.treenode.score
with open(snakemake.output[0], 'w') as snakeout:
    snakeout.write( json.dumps( {snakemake.input[0] : taxscore} ) )
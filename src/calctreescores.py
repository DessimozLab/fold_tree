import json
import treescore
import pandas as pd
import toytree


#calc the taxscore 
uniprot_df = pd.read_csv(snakemake.input[0])
scores = {}
for t in snakemake.input[1:]:
    tree = toytree.tree(t)
    lineages = treescore.make_lineages(uniprot_df)
    tree = treescore.label_leaves( tree , lineages)
    overlap = treescore.getTaxOverlap(tree.treenode)
    taxscore = tree.treenode.score
    scores[t] = taxscore
with open(snakemake.output[0], 'w') as snakeout:
    snakeout.write( json.dumps( scores ) )
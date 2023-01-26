import json
import treescore
import pandas as pd
import toytree
import numpy as np
from scipy.stats import describe

#calc the taxscore 
uniprot_df = pd.read_csv(snakemake.input[0])
scores = {}
stats = {}
for t in snakemake.input[1:]:
    tree = toytree.tree(t)
    lineages = treescore.make_lineages(uniprot_df)
    tree = treescore.label_leaves( tree , lineages)
    overlap = treescore.getTaxOverlap(tree.treenode)
    taxscore = tree.treenode.score
    #calc descriptive stats on normalized branch lens    
    lengths = np.array([node.dist for node in t.treenode.traverse()])
    lengths /= np.sum(lengths)
    scores[t] = {'score': taxscore, 'stats': describe(lengths)}

with open(snakemake.output[0], 'w') as snakeout:
    snakeout.write( json.dumps( scores ) )
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
    tree = treescore.labelwRED(tree.treenode)
    overlap = treescore.getTaxOverlap(tree)

    taxscore = tree.score
    redscore = tree.score_x_red

    #calc descriptive stats on normalized branch lens to see if trees are balanced  
    lengths = np.array([node.dist for node in tree.traverse()])
    lengths /= np.sum(lengths)

    scores[t] = {'score': taxscore, 'stats': describe(lengths) , 'red_x_score': redscore }
with open(snakemake.output[0], 'w') as snakeout:
    snakeout.write( json.dumps( scores ) )
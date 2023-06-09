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
    tree = toytree.tree(t , format = 0)
    lineages = treescore.make_lineages(uniprot_df)
    tree = treescore.label_leaves( tree , lineages)
    #tree = treescore.labelwRED(tree.treenode)
    overlap = treescore.getTaxOverlap(tree.treenode)
    taxscore = tree.treenode.score
    redscore = tree.treenode.score_x_frac
    #calc descriptive stats on normalized branch lens to see if trees are balanced  
    lengths = np.array([node.dist for node in tree.treenode.traverse()])
    lengths /= np.sum(lengths)

    #calc the root first taxscore
    treescore.getTaxOverlap_root(tree.treenode)
    root_score = treescore.sum_rootscore(tree.treenode)

    scores[t] = {'score': taxscore, 'stats': describe(lengths) , 'score_x_frac': redscore , 'root_score': root_score }

with open(snakemake.output[0], 'w') as snakeout:
    snakeout.write( json.dumps( scores ) )
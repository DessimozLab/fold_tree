import toytree
from itertools import combinations
import json

trees = {}
for t in snakemake.input:
    trees[t] = toytree.tree(t)

distances={}
for t1,t2 in combinations(trees, 2):
    distances[t1+'_'+t2] = trees[t1].treenode.robinson_foulds(trees[t2].treenode)[0:2]

with open(snakemake.output[0], 'w' ) as rfout:
    rfout.write(json.dumps(distances))

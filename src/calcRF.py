import toytree
from itertools import combinations
import json

trees = {}
for t in snakemake.input:
    trees[t] = toytree.tree(t)
distances={}
for t1,t2 in combinations(trees, 2):
    distances[(t1,t2)] = trees[t1].treenode.robinson_foulds(trees[2].treenode)
with open(snakemake.output[0]) as rfout:
    rfout.write(json.dumps(distances))

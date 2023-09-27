import foldseek2tree
import toytree


#load the tree
tre = toytree.tree(snakemake.input[0])
#output the 1st tree

tre.write(snakemake.output[0], tree_format=0)

print( 'postprocessing MAD-root done')
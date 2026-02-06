import foldseek2tree
import toytree

treefile = foldseek2tree.postprocess(snakemake.input[0] , snakemake.output[0] ) 
print(treefile , 'done')
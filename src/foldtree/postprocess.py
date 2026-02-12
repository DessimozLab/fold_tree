import foldseek2tree

treefile = foldseek2tree.postprocess(snakemake.input[0] , snakemake.output[0] ) 
print(treefile , 'done')
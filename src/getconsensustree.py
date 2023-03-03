import foldseek2tree
ct = foldseek2tree.consensustree(snakemake.input)
ct.write(snakemake.output[0])
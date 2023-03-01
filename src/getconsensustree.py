import foldseek2tree

ct = foldseek2tree.consensustree(snakemake.input)
ct.write(outfile = snakemake.output[0], format_root_node=True)


with open(snakemake.output[0], 'w') as outfile:
    for fname in snakemake.input:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)
            outfile.write('\n')
            
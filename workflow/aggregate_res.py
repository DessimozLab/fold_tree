with open('output.txt') as snakeout:
	for l in snakemake.input:
		snakeout.write(l)

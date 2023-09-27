import foldseek2tree
import toytree

with open(snakemake.input[0]) as madrootinfile:
    with open(  snakemake.output[0]  , 'w') as madrootoutfile:
        infile = madrootinfile.read()
        trees = infile.split(';\n')
        if len( trees) > 1:
            print('selecting first MAD-root tree')
            madrootoutfile.write( trees[0] + ';\n')
        else:
            madrootoutfile.write( infile  )
print( 'postprocessing MAD-root done')
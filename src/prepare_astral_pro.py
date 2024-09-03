

import sys
sys.path.append( '../src/')
import treescore


speciesmap, speciesmap_file , ncbitree = treescore.prepare_astral_input(snakemake.input.finalset , speciestreeout = snakemake.output.species_tree , mapperout = snakemake.output.map )

print( 'speciesmap' , speciesmap)
print( 'speciesmap_file' , speciesmap_file)
print( 'ncbitree' , ncbitree)
print( 'done' )
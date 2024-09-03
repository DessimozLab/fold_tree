

import json
import sys
sys.path.append( '../src/')
import treescore

df , score = treescore.return_astral_score(snakemake.input.st , snakemake.input.logout )

print( 'score' , score)
print( 'done' )

#cast all values to to float
for k in score:
    score[k] = float(score[k])
    
json.dump( score , open( snakemake.output.json , 'w') )

df.to_csv( snakemake.output.csv , index = False )
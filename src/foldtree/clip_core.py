import json
import corecut
import glob
import pandas as pd

print('cutting core')
corecut.extract_core(snakemake.input[0] , snakemake.output[0] )
print('done')
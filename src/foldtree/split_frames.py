import glob
import os
import sys
from Bio.PDB import PDBParser, PDBIO, Select
import tqdm 
import random

#find all pdb files in the input directory

inputdir   = snakemake.params['input_dir']
outputdir  = snakemake.params['output_dir']
nframes = snakemake.params['nframes']


for i in range(nframes):
    if not os.path.exists(outputdir + 'frame_'+str(i)):
        os.makedirs(outputdir + 'frame_'+str(i))
        os.makedirs(outputdir + 'frame_'+str(i) + '/structs')
    #write an empty identifier file
    with open(outputdir + 'frame_'+str(i)+'/identifiers.txt', 'w') as f:        
        f.write('')

#loop over all pdb files
for pdbfile in tqdm(glob.glob(inputdir + '*.pdb')):
    #read the frames from the pdb file
    pdbid = pdbfile.split('/')[-1].split('.')[0]
    parser = PDBParser()
    structure = parser.get_structure(pdbid, pdbfile)
    frames = [ m for model in structure]
    #shuffle the frames
    frames = random.shuffle(frames) 
    for model in tqdm.tqdm(frames):
        io = PDBIO()
        io.set_structure(model)
        io.save(outputdir + 'frame_'+str(i)+'/structs/'+pdbid+'.pdb')


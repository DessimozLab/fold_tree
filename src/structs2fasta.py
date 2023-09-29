#import bipython pdb parser

import glob
from Bio import PDB
with open(snakemake.output[0], 'w') as f:
    parser = PDB.PDBParser()
    print('/'.join(snakemake.input[0].split('/')[0:-1]) + '/structs/' )
    files = glob.glob('/'.join(snakemake.input[0].split('/')[0:-1])+ '/structs/*.pdb')
    for pdb in files:
        #read pdb file
        structure = parser.get_structure(pdb, pdb)
        #get chain
        chain = structure[0]['A']
        #get residues
        residues = list(chain.get_residues())
        #get sequence
        sequence = ''
        for residue in residues:
            sequence += PDB.Polypeptide.three_to_one(residue.get_resname())
        #write fasta file
            f.write('>' + pdb.split('.')[0] + '\n')
            f.write(sequence + '\n')
import json
import treescore
import pandas as pd
import toytree
import numpy as np
from scipy.stats import describe
import copy


def retspecies_tree(species_set):
    #get ncbi tree of species set
    species_set = list(species_set)
    species_set = [int(s) for s in species_set]
    species_set = ncbi.get_topology(species_set, intermediate_nodes=True)
    return species_set


def calc_scores(t , uniprot_csv , species_tree):
    uniprot_df = pd.read_csv(uniprot_csv)
    tree = toytree.tree(t)
    lineages = treescore.make_lineages(uniprot_df)
    species = treescore.get_species(uniprot_df)
    tree = treescore.label_leaves( tree , lineages , species)
    #tree = treescore.labelwRED(tree.treenode)
    overlap = treescore.getTaxOverlap(tree.treenode)
    taxscore = tree.treenode.score

    #calc descriptive stats on normalized branch lens to see if trees are balanced  
    lengths = np.array([node.dist for node in tree.treenode.traverse()])
    lengths /= np.sum(lengths)
    #calc the root first taxscore

    tree4 = copy.deepcopy(tree)
    treescore.getTaxOverlap_root(tree4.treenode)
    root_score , root_score_nr = treescore.sum_rootscore(tree4.treenode)

    species_set = set()
    #label the leaves with species and number
    for l in tree.treenode.get_leaves():
        if l.sp_num:
            l.name = l.sp_num
            species_set.add(l.sp_num.split('_')[0])


    #change to phylo tree
    ncbitree = PhyloTree(species_tree  , sp_naming_function=None)

    ncbitree.write(outfile=uniprot_csv.replace('.csv','_ncbi_tree.nwk' ) , format=1)
    etetree = PhyloTree(tree.write() , sp_naming_function=None)

    for l in etetree.get_leaves():
        l.species   = l.name.split('_')[0]
    
    recon_tree, events = etetree.reconcile(ncbitree)
    recon_dups = recon_tree.search_nodes(evoltype="D")
    recon_losses = recon_tree.search_nodes(evoltype="L")
    recon_speciations = recon_tree.search_nodes(evoltype="S")
    print( 'algo 1' )
    print('dups:', len(recon_dups))
    print('losses:', len(recon_losses))
    print('speciations:', len(recon_speciations))

    rfs = rf2species(etetree , ot_samples = 10)
    print('RFs:', rfs)


    print( 'algo 2' )
    events = etetree.get_descendant_evol_events()
    dups = etetree.search_nodes(evoltype="D")
    losses = etetree.search_nodes(evoltype="L")
    speciations = etetree.search_nodes(evoltype="S")    
    print('dups:', len(dups))
    print('losses:', len(losses))
    print('speciations:', len(speciations))



    scores = {}
    #measure the distances of leaves to root
    distances = np.array([ node.get_distance(tree.treenode) for node in tree.treenode.get_leaves() ])
    distances_norm = distances / np.mean(distances)
    scores[t] = {'score': taxscore, 'stats': describe(lengths) , 'ultrametricity':  describe(distances), 
                    'ultrametricity_norm':  describe(distances_norm) , 'root_score': root_score , 'root_score_nr': root_score_nr  , 
                    'SO_speciations': len(recon_speciations) , 'SO_dups': len(recon_dups) , 'SO_losses': len(recon_losses) ,
                    'RECON_speciations':speciations ,'RECON_dups': len(dups) , 'RECON_losses': len(losses)  }
    return scores

#calc the taxscore 
uniprot_df = pd.read_csv(snakemake.input[0])
scores = {}
stats = {}

for t in snakemake.input[1:]:
    print(t)
    scores.update(calc_scores(t , snakemake.input[0]))
print(scores)
with open(snakemake.output[0], 'w') as snakeout:
    snakeout.write( json.dumps( scores ) )




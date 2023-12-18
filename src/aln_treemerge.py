import subprocess
import os
import glob
import toytree
import tqdm
import pandas as pd


#recursive function from root to leaves

allvall = snakemake.input[0] 
treefile = snakemake.input[1]


def get_leafset( treenode ):
    """
    this function returns the leafset of a node
    """
    if treenode.is_leaf():
        return [treenode.name]
    else:
        return treenode.get_leaf_names()


#traverse tree from root to leaves
def traverse_tree( treenode , leafset = []):
    """
    this function traverses the tree from root to leaves
    """
    if treenode.is_leaf():
        leafset.append(treenode.name)
        
        return leafset
    else:
        for child in treenode.children:
            leafset = traverse_tree( child , leafset = leafset)
        return leafset

#read all vs all distances

def read_allvall( ):

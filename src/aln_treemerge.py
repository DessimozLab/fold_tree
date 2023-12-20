import subprocess
import os
import glob
import toytree
import tqdm
import pandas as pd
import subprocess
import shlex


#recursive function from root to leaves
allvall = snakemake.input[0] 
treefile = snakemake.input[1]
allvall = pd.read_table(allvall, header = None)
#this should already have the 3di alignment as well as the AA aligntment


def get_leafset( treenode ):
    """
    this function returns the leafset of a node
    """
    if treenode.is_leaf():
        return [treenode.name]
    else:
        return set(treenode.get_leaf_names())

def mergeAlign( fasta1, fasta2, outfile):
	args =  'clustalo --p1 ' + fasta1 + ' --p2 ' + fasta2 + ' --force -o ' + outfile
	args = shlex.split(args)
	p = subprocess.call(args )

def sub2fasta( sub, outfile , fastacol1='qaln' , fastacol1='taln' ):
    with open(outfile, 'w') as f:
        f.write('>' + sub['query'] + '\n')
        f.write(sub[fastacol1] + '\n')
        f.write('>' + sub['target'] + '\n')
        f.write(sub[fastacol2] + '\n')    
    return outfile

def retalns(allvsall, leafname,leafset):
    sub = allvsall[allvsall['query'] == leafname]
    sub = sub[sub['target'].isin(leafset)]
    sub = sub[sub['query'] != sub['target']]
    return sub.iloc[0]

#traverse tree from root to leaves recursively
def traverse_tree_merge( treenode , resAA , res3di, allvall ):
    """
    this function traverses a tree from root to leaves recursively
    it returns a dictionary with the distances between nodes
    """
    childalns3di = []
    childalnsAA = []
    treenode.leafset = get_leafset(treenode)
    for c in treenode.get_children():
        c.leafset = get_leafset(c)
        if not c.aln:
            if len(c.leafset) == 2:
                #cherry, just grab the right pairwise alignment
                sub = retalns(allvall, list(leafset)[0] , leafset)
                c.aln = sub2fasta(sub, c.name + '_inter.fasta')
                c.aln3di = sub2fasta(sub, c.name + '_inter.fasta' , fastacol1='qaln3di' , fastacol1='taln3di')
                childalnsAA.append(c.aln)
                childalns3di.append(c.aln3di)
            if c.aln and c.aln3di:
                #if the alignment is already there, then we don't need to do anything
                childalnsAA.append(c.aln)
                childalns3di.append(c.aln3di)
            else:
                if c.is_leaf():    
                    #if the node is a leaf, then we need to add it to the alignment
                    sub = retalns(allvall, leafname, treenode.leafset)
                    c.aln = sub2fasta(sub, leafname + '_inter.fasta')
                    c.aln3di = sub2fasta(sub, leafname + '_inter.fasta' , fastacol1='qaln3di' , fastacol1='taln3di')
                    childalnsAA.append(c.aln)
                    childalns3di.append(c.aln3di)
                else:
                    #if the node is not a leaf, then we need to recursively call the function
                    childalnsAA.append(traverse_tree(c, resAA , res3di, allvall)[0])
                    childalns3di.append(traverse_tree(c, resAA , res3di, allvall)[1])        
        treenode.aln = mergeAlign( treenode.name + '_inter.fasta' , childalnsAA[0], childalnsAA[1]))            
        treenode.aln3di = mergeAlign( treenode.name + '_inter.3di.fasta' , childalns3di[0], childalns3di[1]))
    return treenode.aln, treenode.aln3di

def remove_redundant( alignment ):
    """
    this function removes redundant sequences from an alignment
    """
    aln = SeqIO.parse(alignment, 'fasta')
    seqs = []
    for s in aln:
        seqs.append(s)
    seqs = list(set(seqs))
    with open(alignment, 'w') as f:
        for s in seqs:
            f.write('>' + s.id + '\n')
            f.write(str(s.seq) + '\n')
    return alignment



#remove all alns except the final merged one
def cleanup( filedir ):
    """
    this function removes all alns except the final merged one
    """
    for f in glob.glob(filedir + '*inter.fasta'):
        os.remove(f)
    


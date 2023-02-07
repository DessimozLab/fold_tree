#foldssek cluster
import subprocess ,shlex
import numpy as np
from scipy.spatial.distance import cdist
import toytree
import pandas as pd


#smooth distmat with MDS
def MDS_smooth(distmat):
    mds = MDS(n_components=int(distmat.shape[0]/2) )#, dissimilarity='precomputed'  )
    distmat = mds.fit_transform(1-distmat )
    distmat = cdist(distmat,distmat, 'minkowski', p=1.5 )
    return distmat

def runargs(args):
    args = shlex.split( args)
    p = subprocess.run( args )
    return p
    
def runFoldseekdb(folder , outfolder):
    args = 'foldseek createdb '+  folder + ' '+ outfolder+'structblobDB '
    p = runargs(args)
    return outfolder+'structblobDB '

def runFoldseek_allvall(dbpath , outfolder):
    args = 'foldseek search '+  dbpath +' '  + dbpath + ' ' + outfolder+'aln tmp -a'
    p = runargs(args)
    args = 'foldseek createtsv '+  dbpath +' '  + dbpath  + ' ' + outfolder+'aln '  + outfolder +'aln_score.tsv'
    p = runargs(args)
    return outfolder +'aln_score.tsv'

def kernelfun(AA,BB, AB):
    return AA + BB - 2*AB

def runFastme( fastmepath , clusterfile ):
    args =  fastmepath +  ' -i ' + clusterfile + ' -o ' + clusterfile+'_tree.txt -n '
    p = runargs(args)
    return clusterfile+'_tree.txt'

def distmat_to_txt( identifiers , distmat, outfile):
    #write out distmat in phylip compatible format
    outstr = str(len(identifiers)) + '\n'
    for i,pdb in enumerate(identifiers):
        namestr = pdb.replace('.','')
        outstr += namestr+ ' ' + np.array2string( distmat[i,:], formatter={'float_kind':lambda x: "%.2f" % x}).replace('[', '').replace(']', '').replace('\n', '')  + '\n'
    with open(outfile , 'w') as handle:
        handle.write(outstr)
        handle.close()
    return outfile


def postprocess(t, delta=10**-10 ):
    #make negative branch lengths a small delta instead
    tre = toytree.tree(t)

    for n in tre.treenode.traverse():
        if n.dist< 0:
            n.dist = delta
    with open(t + '.PP.nwk' , 'w') as treeout:
        treeout.write(tre.write())
    return t + 'PP.nwk'


def structblob2tree(input_folder, logfolder):
    dbpath = runFoldseekdb(input_folder, logfolder)
    alnres = runFoldseek_allvall(dbpath , logfolder)
    res = pd.read_table(alnres, header = None ,delim_whitespace=True)
    res[0] = res[0].map(lambda x :revmapper[x.replace('.pdb', '')])
    res[1] = res[1].map(lambda x :revmapper[x.replace('.pdb', '')])
    ids = list( set(list(res[0].unique()) + list(res[1].unique())))
    pos = { protid : i for i,protid in enumerate(ids)}
    self_dists = res[res[0] == res[1]]
    self_distmap = dict(zip(self_dists[0] , self_dists[2] ) )    
    kernel_distmat = np.zeros((len(pos), len(pos)))
    distmat = np.zeros((len(pos), len(pos)))
    for idx,row in res.iterrows():
        kernel_distmat[pos[row[0]] , pos[row[1]]] = kernelfun(self_distmap[row[0]] , self_distmap[row[1]] , row[2])
        distmat[pos[row[0]] , pos[row[1]]]= row[3]

    distmat_txt = distmat_to_txt( ids , kernel_distmat , logfolder , prefix ='standard' )
    out_tree = runFastme( 'fastme' , distmat_txt )

    
    distmat_txt = distmat_to_txt( ids , distmat , logfolder , prefix ='rawscore' )
    out_tree = runFastme( 'fastme' , distmat_txt )
    
    return out_tree

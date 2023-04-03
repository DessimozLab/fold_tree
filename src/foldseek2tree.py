#foldssek cluster
import subprocess ,shlex
import numpy as np
from scipy.spatial.distance import cdist
import statsmodels
import toytree
import pandas as pd
import re
from scipy.stats import chi2


def consensustree(treelist):
    '''Get a consensus tree from a list of tree files
    
    Parameters
    ----------
    treelist : list
        list of tree files

    '''
    #creat a multitree object from the list of tree files
    treelist = [toytree.tree(i , format = 0 ) for i in treelist]
    mt = toytree.mtree(treelist)
    #get the consensus tree
    ct = mt.get_consensus_tree( )
    return ct

#smooth distmat with MDS
def MDS_smooth(distmat):
    mds = MDS(n_components=int(distmat.shape[0]/2) )#, dissimilarity='precomputed'  )
    distmat = mds.fit_transform(1-distmat )
    distmat = cdist(distmat,distmat, 'minkowski', p=1.5 )
    return distmat

def Tajima_dist(kn_ratio,bfactor=1,iter = 100):
    taj = np.add.reduce([ (kn_ratio**(np.ones(kn_ratio.shape)*i) )/ (bfactor**(i-1)*i) for i in range(1,iter) ] )
    #set diagonal to 0
    np.fill_diagonal(taj, 0)
    return taj

def clock_test(distmats, ntriplets = 1000):
    '''test if a distance matrix is clocklike using a xi squared test for n random triplets'''
    #get n random triplets without replacement
    triplets = np.random.choice(distmat.shape[0], size=(ntriplets, 3), replace=False)
    #get the distances between the triplets
    res = {}
    for mattype in distmats:
        distmat = distmats[mattype]

        dists = np.array( [ sorted[distmat[i,j],distmat[i,k],distmat[j,k]] for i,j,k in triplets ] )
        x2 = (dists[:,0] - dists[:,1])**2 / dists[:,0] + dists[:,1]
        pvals = chi2.cdf(x2, 1)

        #apply bonferroni correction
        corrected = statsmodels.stats.multitest.multipletests(pvals, alpha=0.05, method='hs', maxiter=1, is_sorted=False, returnsorted=False)
        res[mattype] = corrected

    return res

def runargs(args):
    '''run a command line command
    
    Parameters
    ----------
    args : str
        command line command
    '''
    
    args = shlex.split( args)
    p = subprocess.run( args )
    return p
    
def runFoldseekdb(folder , outfolder , foldseekpath = '../foldseek/bin/foldseek'):
    '''run foldseek createdb
    
    parameters
    ----------
    folder : str
        path to folder with pdb files
    outfolder : str 
        path to output folder
    

    '''
    args = foldseekpath + ' createdb '+  folder + ' '+ outfolder+'structblobDB '
    p = runargs(args)
    return outfolder+'structblobDB '

def runFoldseek_allvall( structfolder , outfolder , foldseekpath = '../foldseek/bin/foldseek' , maxseqs = 3000):
    '''
    run foldseek search and createtsv
    
    parameters
    ----------
    dbpath : str
        path to foldseek database
    outfolder : str 
        path to output folder
    maxseqs : int   
        maximum number of sequences to compare to

    '''
    
    foldseekpath + " easy-search " + structfolder + " "+ structfolder +" "+ outfolder+"/allvall.csv " +  structfolder+"/tmp --format-output 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,lddtfull,alntmscore' --exhaustive-search --alignment-type 2" 

    return outfolder +'aln_score.tsv'

def runFoldseek_allvall_EZsearch(infolder , outpath , foldseekpath = '../foldseek/bin/foldseek'):
    '''
    run foldseek easy-search
    
    parameters
    ----------
    infolder : str
        path to folder with pdb files
    outpath : str
        path to output folder
    foldseekpath : str  
        path to foldseek binary

        '''
    
    args = foldseekpath + ' easy-search ' + infolder + ' ' + infolder +' '+ outpath + " tmp --format-output 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,lddtfull,alntmscore' --format-mode 3 --exhaustive-search "
    p = runargs(args)
    return outpath

def kernelfun(AA,BB, AB):
    return AA + BB - 2*AB

def runFastme( fastmepath , clusterfile ):
    '''run fastme
    
    parameters
    ----------
    fastmepath : str
        path to fastme binary
    clusterfile : str
        path to all vs all distance matrix in fastme format
    '''

    args =  fastmepath +  ' -i ' + clusterfile + ' -o ' + clusterfile+'_tree.txt -n '
    p = runargs(args)
    return clusterfile+'_tree.txt'

def distmat_to_txt( identifiers , distmat, outfile):
    '''
    write out a distance matrix in fastme format

    Parameters
    ----------
    identifiers : list
        list of identifiers for your proteins
    distmat : np.array  
        distance matrix
    outfile : str   
        path to output file

    '''

    #write out distmat in phylip compatible format
    outstr = str(len(identifiers)) + '\n'
    for i,pdb in enumerate(identifiers):
        outstr += pdb + ' ' + ' '.join( [ "{:.4f}".format(d) for d in list( distmat[i,: ] )  ]  ) + '\n'
    with open(outfile , 'w') as handle:
        handle.write(outstr)
        handle.close()
    return outfile
   
def postprocess(t, outree, delta=0 ):
    '''
    postprocess a tree to make sure all branch lengths are positive
    
    Parameters
    ----------
    t : str
        path to tree file
    delta : float
        small number to replace negative branch lengths with'''
    #make negative branch lengths a small delta instead
    with open(t) as treein:
        treestr = ' '.join( [ i.strip() for i in treein ] )

    tre = toytree.tree(treestr , format = 0 )
    print(tre)

    for n in tre.treenode.traverse():
        if n.dist< 0:
            n.dist = delta
    tre.write( outree, tree_format = 0 )

    return outree


def structblob2tree(input_folder, logfolder):
    '''run structblob pipeline

    Parameters
    ----------
    input_folder : str
        path to folder with pdb files
    logfolder : str 
        path to output folder
    '''
    
    dbpath = runFoldseekdb(input_folder, logfolder)
    alnres = runFoldseek_allvall(dbpath , logfolder)
    res = pd.read_table(alnres, header = None ,delim_whitespace=True)
    res[0] = res[0].map(lambda x :x.replace('.pdb', ''))
    res[1] = res[1].map(lambda x :x.replace('.pdb', ''))

    ids = list( set(list(res[0].unique()) + list(res[1].unique())))
    pos = { protid : i for i,protid in enumerate(ids)}
    distmat = np.zeros((len(pos), len(pos)))
    for idx,row in res.iterrows():
        distmat[pos[row[0]] , pos[row[1]]]= row[3]    
    distmat = 1-distmat
    distmat += distmat.T
    distmat /= 2
    
    distmat_txt = distmat_to_txt( ids , distmat , input_folder + 'fastme_mat.txt'  )
    out_tree = runFastme( 'fastme' , distmat_txt )
    return out_tree

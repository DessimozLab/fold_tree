#foldssek cluster
import subprocess ,shlex
import numpy as np
from scipy.spatial.distance import cdist
import statsmodels
import toytree
import pandas as pd
import re
import os
from scipy.stats import chi2
import argparse

def consensustree(treelist):
	'''get a consensus tree from a list of tree files
	
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

def Tajima_dist( kn_ratio,bfactor=19/20, iter = 100 ):
	taj = np.add.reduce([ (kn_ratio**(np.ones(kn_ratio.shape)*i) )/ (bfactor**(i-1)*i) for i in range(1,iter) ] )
	#set diagonal to 0
	np.fill_diagonal(taj, 0)
	return taj

def simple_logdist( kn_ratio,bfactor=19/20, iter = 100 ):
	d = -bfactor* np.log(1- ( kn_ratio / bfactor ) )
	#set diagonal to 0
	np.fill_diagonal(d, 0)
	return d



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
	
	args = foldseekpath + ' easy-search ' + infolder + ' ' + infolder +' '+ outpath + " tmp --format-output 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,lddtfull,alntmscore' --exhaustive-search "
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

def runQuicktree(   clusterfile , quicktreepath= 'quicktree' ):
	'''
	run quicktree

	parameters
	----------
	clusterfile : str
		path to all vs all distance matrix in fastme format
	quicktreepath : str 
		path to quicktree binary

	'''
	args = quicktreepath + ' -i m ' + clusterfile +' > ' + clusterfile + '.struct_tree.nwk'
	p = runargs(args)
	return clusterfile + '.struct_tree.nwk'


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

	tre = toytree.tree(treestr )
	print(tre)

	for n in tre.treenode.traverse():
		if n.dist< 0:
			n.dist = delta
	tre.write( outree, tree_format = 0 )

	return outree



def read_dbfiles3di(folder, name = 'outdb'):
    #find positions 
    threeDidb = folder+ name +'_ss'
    threeDiseq = [ l.strip().replace('\x00','') for l in open(threeDidb)]
    
    
    lookup = folder+name+'.lookup'
    ids = [ l.split()[1].strip() for l in open(lookup)]
    
    AADB = folder+name
    AAs = [ l.strip().replace('\x00','') for l in open(AADB)]

    mapper3di = dict(zip(ids,threeDiseq))
    mapperAA = dict(zip(ids,AAs))
    
    return mapper3di, mapperAA


def structblob2tree(input_folder, outfolder, overwrite = False, fastmepath = 'fastme', quicktreepath = 'quicktree' , foldseekpath = '../foldseek/foldseek' , delta = 0.0001):
	'''run structblob pipeline for a folder of pdb files without snakemake

	Parameters
	----------
	input_folder : str
		path to folder with pdb files
	logfolder : str 
		path to output folder
	'''
	
	#check if the foldseek output is already there
	if os.path.exists(outfolder + 'res.m8') and overwrite == False:
		print('found foldseek output, skipping foldseek')
		alnres = outfolder + 'res.m8'
	else:
		alnres = runFoldseek_allvall_EZsearch(input_folder , outfolder + 'res.m8', foldseekpath = foldseekpath)
	
	res = pd.read_table(alnres , header = None )
	res[0] = res[0].map(lambda x :x.replace('.pdb', ''))
	res[1] = res[1].map(lambda x :x.replace('.pdb', ''))
	res.columns = 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,lddtfull,alntmscore'.split(',')
	ids = list( set(list(res['query'].unique()) + list(res['target'].unique())))
	pos = { protid : i for i,protid in enumerate(ids)}
	kernels = ['fident', 'alntmscore', 'lddt']
	matrices = { k:np.zeros((len(pos), len(pos))) for k in kernels }
	print(res)

	#calc kernel for tm, aln score, lddt
	for idx,row in res.iterrows():
		for k in matrices:
			matrices[k][pos[row['query']] , pos[row['target']]] += row[k]
			matrices[k][pos[row['target']] , pos[row['query']]] += row[k]

	trees = {}
	for i,k in enumerate(matrices):
		matrices[k] /= 2
		matrices[k] = 1-matrices[k]

		if k == 'fident':
			distmat = Tajima_dist( matrices[k] , bfactor = .93 , iter = 100 )
		else:
			distmat = Tajima_dist( matrices[k] )
		print(distmat,  np.amax(distmat), np.amin(distmat) )

		np.save( outfolder + k + '_distmat.npy' , distmat)
		distmat_txt = distmat_to_txt( ids , distmat , outfolder + k + '_distmat.txt' )
		out_tree = runFastme(  fastmepath = fastmepath , clusterfile = distmat_txt )
		out_tree = postprocess(out_tree, out_tree+'.pp.nwk' , delta = delta)
		trees[k] = out_tree
	return alnres, trees
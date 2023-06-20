
from . import foldseek2tree
from . import corecut 
import argparse
import numpy as np
import pandas as pd
import os

def structblob2tree(input_folder, outfolder, overwrite = False,
             fastmepath = 'fastme', quicktreepath = 'quicktree' , 
         foldseekpath = 'foldseek' , delta = 0.0001 ,
       correction = False , kernel = 'fident' , core = False 
       , hittresh = .8 , minthresh = .6):
    '''
    run fold tree pipeline for a folder of pdb files
    
    Parameters
    ----------
    input_folder : str
        path to folder with pdb files   
    outfolder : str
        path to output folder   
    overwrite : bool
        overwrite existing foldseek output  
    fastmepath : str    
        path to fastme executable
    quicktreepath : str 
        path to quicktree executable
    foldseekpath : str  
        path to foldseek executable 
    delta : float   
        small number to replace negative branch lengths with, default is .0001
    correction : str    
        correction method to use, either 'tajima' or 'none'
    kernel : str    
        kernel to use, either 'fident', 'lddt' or 'alntmscore'
    
    '''
    #check if the foldseek output is already there
    if os.path.exists(outfolder + 'res.m8') and overwrite == False:
        print('found foldseek output, skipping foldseek')
        alnres = outfolder + 'res.m8'
    else:
        alnres = foldseek2tree.runFoldseek_allvall_EZsearch(input_folder , outfolder + 'res.m8', foldseekpath = foldseekpath)
    
    if core == True:
        corecut.extract_core( alnres , resdf+'.core.csv',  hitthresh = .8 ,minthresh = .6, corefolder = input_folder+'core_structs/' , structfolder = input_folder )
        if os.path.exists(outfolder + 'res.m8') and overwrite == False:
            print('found foldseek core output, skipping foldseek')
            alnres = outfolder + 'core.res.m8'
        else:
            alnres = foldseek2tree.runFoldseek_allvall_EZsearch(input_folder , outfolder + 'core.res.m8', foldseekpath = foldseekpath)
    
    res = pd.read_table(alnres , header = None )
    res[0] = res[0].map(lambda x :x.replace('.pdb', ''))
    res[1] = res[1].map(lambda x :x.replace('.pdb', ''))
    res.columns = 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,lddtfull,alntmscore'.split(',')
    
        
    ids = list( set(list(res['query'].unique()) + list(res['target'].unique())))
    pos = { protid : i for i,protid in enumerate(ids)}
    matrices = { kernel:np.zeros((len(pos), len(pos)))  }
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
        matrices[k] = np.clip(matrices[k], 0, 1)
        
        print(matrices[k], np.amax(matrices[k]), np.amin(matrices[k]) )
        if correction:
            if kernel == 'fident':
                factor = 19/20
            else:
                factor = 1
            matrices[k] = Tajima_dist(matrices[k], factor = factor)
        np.save( input_folder + k + '_distmat.npy' , matrices[k])
        distmat_txt = foldseek2tree.distmat_to_txt( ids , matrices[k] , outfolder + k + '_distmat.txt' )
        out_tree = foldseek2tree.runFastme(  fastmepath = fastmepath , clusterfile = distmat_txt )
        out_tree = foldseek2tree.postprocess(out_tree, input_folder + 'structblob_tree.nwk' , delta = delta)
        trees[k] = out_tree
    return alnres, trees

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='run foldtree pipeline for a folder of pdb files')
    parser.add_argument('struct_dir', help='path to folder with pdb files')
    parser.add_argument('output_dir', help='output directory')
    parser.add_argument('--corecut', help='cut the core of the proteins and realign')
    parser.add_argument('--kernel', choices = ['lddt', 'fident', 'tmalign' ] ,
    default = 'fident',
    help='this is the comparison metric used to build the tree')
    parser.add_argument('--fastmepath', default='fastme', help='path to fastme binary')
    parser.add_argument('--quicktreepath', default='quicktree', help='path to quicktree binary')
    parser.add_argument('--foldseekpath', default='../foldseek/foldseek', help='path to foldseek binary')
    parser.add_argument('--delta', default=0.0001, help='small number to replace negative branch lengths with')
    parser.add_argument('--correction', help='use the -ln correction for the distance matrix')
    parser.add_argument('--overwrite', help='overwrite existing foldseek output')
    parser.add_argument( 'hittresh', default = .8, help='threshold for finding the boundaries of the core')
    parser.add_argument( 'minthresh', default = .6, help='threshold if the core is not found')
    
    args = parser.parse_args()

    if not all([args.positional1, args.positional2]):
        parser.error("Positional arguments are required.")
    flag_dict = {}
    flag_dict['corecut'] = True if args.corecut else False
    flag_dict['hittresh'] = args.hittresh
    flag_dict['minthresh'] = args.minthresh
    flag_dict['correction'] = True if args.correction else False
    flag_dict['overwrite'] = True if args.overwrite else False
    flag_dict['fastmepath'] = args.fastmepath
    flag_dict['quicktreepath'] = args.quicktreepath
    flag_dict['foldseekpath'] = args.foldseekpath
    flag_dict['delta'] = args.delta
    flag_dict['kernel'] = args.kernel
    print(f"struct dir: {args.struct_dir}, output dir: {args.output_dir}")
    structblob2tree(args.struct_dir, args.output_dir, **flag_dict)
    print('Done!')


    

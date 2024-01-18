import subprocess
import os
import glob
import toytree
import tqdm
import pandas as pd
from subprocess import PIPE, Popen
import shlex
from Bio import SeqIO
import random
import itertools
import numpy as np

def Fident(str1,str2 , verbose = False):
    #minlen= min( (len(str1),len(str2))  )
    #str1 = str1[:minlen]
    #str2 = str2[:minlen]
    str1 = np.array(list(str1))
    str2 = np.array(list(str2))            
    return len(np.where( (str1 == str2 ) & (str1 != '-' ) & (str2 != '-')  )[0]) / len(str1)

def copyaln( aln, seq):
    seqiter = iter(seq)
    newaln = ''
    for i,char in enumerate(aln):
        if char == '-':
            newaln += '-'
        else:
            newaln+=next(seqiter)
    return newaln

def read_dbfiles3di(  AADB , threeDidb):
    #find positions 
    threeDiseq = [ l.strip().replace('\x00','') for l in open(threeDidb)]
    lookup = AADB+'.lookup'
    ids = [ l.split()[1].strip() for l in open(lookup)]
    AAs = [ l.strip().replace('\x00','') for l in open(AADB)]

    mapper3di = dict(zip(ids,threeDiseq))
    mapperAA = dict(zip(ids,AAs))
    
    return mapper3di, mapperAA

def calc_fident_crossaln(row , verbose = False):
    #amino acid representations of alns using AAand3di or just 3di
    qaln_2, taln_2 = row.qaln , row.taln
    #start and stop of aln
    
    qstart_2, qend_2, tstart_2 , tend_2 = row.qstart, row.qend , row.tstart , row.tend
    #indexing starts at 1...
    
    #3di of the query and target
    structQ, structT = row['3diq'], row['3dit']
    AAq, AAt = row['AAq'], row['AAt']

    #add gaps
    t3diAA_newgaps = copyaln(taln_2, structT[tstart_2-1:tend_2]) 
    q3diAA_newgaps = copyaln(qaln_2, structQ[qstart_2-1:qend_2])
    row = pd.Series( { '3di_qaln_mode2':q3diAA_newgaps , '3di_taln_mode2':t3diAA_newgaps })
    #return columns
    return row

def get_leafset( treenode ):
    """
    this function returns the leafset of a node
    """
    if treenode.is_leaf():
        return [treenode.name]
    else:
        return treenode.get_leaf_names()

def mergeAlign( fasta1, fasta2, outfile):
    print(fasta1,fasta2)
    #args =  'clustalo --p1 ' + fasta1 + ' --p2 ' + fasta2 + ' --force -o ' + outfile
    print('merge' , fasta1 , fasta2 )
    args = 'mafft --add ' + fasta1 + ' ' + fasta2 
    args = shlex.split(args)
    p = Popen( args, stdin=PIPE, stdout=PIPE)
    aln = p.communicate()[0].decode()
    print('mafft aln',aln)
    with open( outfile , 'w') as alnout:
        alnout.write(aln)
    return outfile


def convolve_strings(str1, str2):
    # Determine the lengths of the strings
    len1, len2 = len(str1), len(str2)

    if len(str1) < len(str2):
        str1, str2 = str2, str1
        len1, len2 = len2, len1
    
    max_alignment = 0
    max_count = 0

    # Slide str2 over str1, starting with one character overlap
    # and continue until str2 is again overlapping by just one character
    for i in range(-len2 + 1, len1):
        count = 0
        for j in range(len2):
            if 0 <= i + j < len1 and str1[i + j] == str2[j]:
                count += 1

        if count > max_count:
            max_count = count
            max_alignment = i

    return max_alignment, max_count

def mergealns( aln1f, aln2f, outfile , verbose = False):
    #find sequences in common between the two alignments
    aln1 = SeqIO.parse(aln1f, 'fasta')
    aln2 = SeqIO.parse(aln2f, 'fasta')
    ids1d = {}
    ids2d = {}

    for s in aln1:
        if s.id not in ids1d:
            ids1d[s.id] = s
    for s in aln2:
        if s.id not in ids2d:
            ids2d[s.id]= s
    
    #find the intersection of the two sets
    ids1 = set(ids1d.keys())
    ids2 = set(ids2d.keys())
    commonids = ids1.intersection(ids2)
    try:
        commonid = list(commonids)[0]
    except:
        print('no common ids')
        print('ids1', ids1)
        print('ids2', ids2)
        
        raise Exception('no common ids')
    #select one common sequence
    coidx1 = -1
    coidx2 = -1
    #find the start of the common sequence in each alignment
    #select a random common sequence
    #create an index of the sequence ids in each alignment
    ids1 = {}
    ids2 = {}
    aln1 = SeqIO.parse(aln1f, 'fasta')
    aln2 = SeqIO.parse(aln2f, 'fasta')
    idlist = []
    for i,s in enumerate(aln1):
        ids1[s.id] = i
        idlist.append(s.id)    
    for i,s in enumerate(aln2):
        ids2[s.id] = i
        idlist.append(s.id)


    #transform both alignments into numpy matrices
    aln1 = SeqIO.parse(aln1f, 'fasta')
    aln2 = SeqIO.parse(aln2f, 'fasta')
    aln1 = np.array([ list(str(s.seq)) for s in aln1])
    aln2 = np.array([ list(str(s.seq)) for s in aln2])
    
    nrows1 = aln1.shape[0]
    nrows2 = aln2.shape[0]

    #generate a list of columns
    aln1 = iter([ aln1[:,i] for i in range(aln1.shape[1])])
    aln2 = iter([ aln2[:,i] for i in range(aln2.shape[1])])
    
    s1 = ids1d[commonid]
    s2 = ids2d[commonid]

    s1raw = iter(str(s1.seq))
    s2raw = iter(str(s2.seq))


    s1 = str(s1.seq).replace('-','')
    s2 = str(s2.seq).replace('-','')

    #if the common subsequence is not found start by removing the first character of the common sequence
    char1 = None
    char2 = None

    #convolution of the two sequences
    maxaln, maxcount = convolve_strings(s1,s2)
    if len(s1) < len(s2):
        if maxaln < 0:
            coidx1 = maxaln
            coidx2 = 0
        else:
            coidx1 = 0
            coidx2 = maxaln
    else:
        if maxaln < 0:
            coidx1 = 0
            coidx2 = maxaln
        else:
            coidx1 = maxaln
            coidx2 = 0
    
    if verbose == True:
        print('convolution', maxaln, maxcount)
        print('coidx1', coidx1)
        print('coidx2', coidx2)


    while np.abs(coidx2) > 0:
        char2 = next(s2raw)
        if verbose == True:
            pass
            #print('discard char', char2)
        if char2 == '-':
            pass
        else:
            coidx2 -= np.sign(coidx2)
        discard = next(aln2)
    else:
        char2 = next(s2raw)

    while np.abs(coidx1) > 0:
        char1 = next(s1raw)
        if verbose  == True:
            pass
            #print('discard char', char1)
        if char1 == '-':
            pass
        else:
            coidx1 -= np.sign(coidx1)
        discard = next(aln1)
    else:
        char1 = next(s1raw)
    
    
    newaln1 = []
    newaln2 = []

    if char1 is None:   
        char1 = next(s1raw)
    if char2 is None:
        char2 = next(s2raw)


    while True:
        try:
            if char1 == '-' and char2 == '-':
                #no pivot information, skip
                char1 = next(s1raw)
                char2 = next(s2raw)
            elif char2 == '-' and char1 != '-':
                char2 = next(s2raw)
                #create a column of gaps
                gaps = np.array(['-'] * nrows1)
                newaln1.append(gaps)
                newaln2.append(next(aln2))
                
            elif char1 == '-' and char2 != '-':
                char1 = next(s1raw)
                #create a column of gaps
                gaps = np.array(['-'] * nrows2)
                newaln1.append(next(aln1))
                newaln2.append(gaps)

              
            elif char1 == char2:
                #match. append both columns
                char1 = next(s1raw)
                char2 = next(s2raw)
                newaln1.append(next(aln1))
                newaln2.append(next(aln2))
            else:
                #mismatch. 
                #should not happen
                print('mismatch')
                print('char1', char1)
                print('char2', char2)
                print('convolution', maxaln, maxcount)

                print(''.join([c for c in s1raw]))
                print(''.join([c for c in s2raw]))
                #raise exception
                raise Exception('mismatch')
        
        except StopIteration:
            #if one of the sequences is finished, then break
            cut = min(len(newaln1), len(newaln2))
            newaln1 = newaln1[:cut]
            newaln2 = newaln2[:cut]
            newaln1 = np.vstack(newaln1).T
            newaln2 = np.vstack(newaln2).T

            break
    newaln = np.concatenate((newaln1, newaln2), axis = 0)
    #write out the new alignment
    with open(outfile, 'w') as f:
        for i in range(newaln.shape[0]):
            #print('>' + idlist[i] + '\n' + ''.join(list(newaln[i,:])) + '\n')
            f.write('>' + idlist[i] + '\n')
            f.write(''.join(list(newaln[i,:])) + '\n')
    remove_redundant( outfile )
    #if verbose == True:
    with open(outfile, 'r') as f:
        print(f.read())
        
    return outfile  

def sub2fasta( sub, outfile , fastacol1='qaln' , fastacol2='taln' ):
    with open(outfile, 'w') as f:
        f.write('>' + sub['query'] + '\n')
        f.write(sub[fastacol1] + '\n')
        f.write('>' + sub['target'] + '\n')
        f.write(sub[fastacol2] + '\n')    
    return outfile

def retalns(allvall, leafname1,leafname2):
    sub = allvall[allvall['query'].isin( leafname1)]
    sub = sub[sub['target'].isin(leafname2)]
    sub = sub[sub['query'] != sub['target']]
    if len(sub)==0:
        print(leafname1, leafname2)
        raise Exception('no sub')
    return sub.iloc[0]

#traverse tree from root to leaves recursively
def traverse_tree_merge( treenode, topleafset, allvall , alnfolder , verbose = False):
    """
    this function traverses a tree from root to leaves recursively
    it returns a dictionary with the iteratively built alignment
    """
    if verbose == True:
        print('traverse', treenode.name , treenode.is_leaf() , treenode.leafset)
    
    if treenode.is_leaf():
        #if the node is a leaf, then we need to add it to the alignment with one of the pivots in the current leafset
        sub = retalns(allvall, [treenode.name] , topleafset )
        treenode.aln = sub2fasta(sub, alnfolder + treenode.name + '_inter.fasta')
        treenode.aln3di = sub2fasta(sub, alnfolder + treenode.name + '_inter.3di.fasta' , fastacol1='3di_qaln_mode2' , fastacol2='3di_taln_mode2')
        return treenode.aln, treenode.aln3di
    
    else:
        childalns3di = {}
        childalnsAA = {}
        bridges3di = {}
        bridgesAA = {}
        #treenode.leafset = get_leafset(treenode)
        #get the intersection of the child leafsets
        treenode.leafset = get_leafset(treenode)
        children = treenode.get_children()
        
        if len(children) == 2 and children[0].is_leaf() and children[1].is_leaf():
            #treat the case of a cherry
            print('cherry', children[0].name , children[1].name)
            treenode.aln = sub2fasta( retalns(allvall, [children[0].name] , [children[1].name]) , alnfolder + treenode.name + '_inter.fasta')
            treenode.aln3di = sub2fasta( retalns(allvall, [children[0].name] , [children[1].name]) , alnfolder + treenode.name + '_inter.3di.fasta' , fastacol1='3di_qaln_mode2' , fastacol2='3di_taln_mode2')
            return treenode.aln, treenode.aln3di
        
        else:
            #not a cherry. one or both sides is a subtree

            print('not cherry', treenode.name  )
            print( 'children', [c.name for c in children])
            for c in treenode.get_children():
                #make sub aln for each child
                if verbose == True:
                    print('traverse', c.name , c.is_leaf() , c.leafset)
                if not c.aln:
                    c.aln,c.aln3di = traverse_tree_merge(c , treenode.leafset , allvall, alnfolder , verbose = verbose)
                childalnsAA[c] = c.aln
                childalns3di[c] = c.aln3di
            

            for c1,c2 in itertools.combinations(treenode.get_children(),2):
                bridgeleaf1 = get_leafset(c1)
                bridgeleaf2 = get_leafset(c2)
                bridge = retalns(allvall, bridgeleaf1 , bridgeleaf2 )
                print('bridge', bridge)
                bridgesAA[(c1,c2)] = sub2fasta(bridge, alnfolder + treenode.name + '_bridge.fasta')
                bridges3di[(c1,c2)] = sub2fasta(bridge, alnfolder + treenode.name + '_bridge.3di.fasta' , fastacol1='3di_qaln_mode2' , fastacol2='3di_taln_mode2')
                if verbose == True:
                    print('bridgeleaf1', bridgeleaf1)
                    print('bridgeleaf2', bridgeleaf2)


            #successively merge the alignments of the children
            for i, cpair in enumerate( itertools.combinations(treenode.get_children(),2)):
                c1,c2 = cpair
                if verbose == True:
                    print('merge', c1.name , c2.name)

                bridgeAA = bridgesAA[(c1,c2)]
                bridge3di = bridges3di[(c1,c2)]
                
                if i == 0:
                    try:
                        alnAA = mergealns( childalnsAA[c1], bridgeAA , alnfolder + treenode.name + '_inter.fasta' , verbose=verbose )
                        aln3di = mergealns( childalns3di[c1], bridge3di , alnfolder + treenode.name + '_inter.3di.fasta' , verbose=verbose )
                        
                        alnAA = mergealns( childalnsAA[c2], alnAA , alnfolder + treenode.name + '_inter.fasta' , verbose=verbose)
                        aln3di = mergealns( childalns3di[c2], aln3di , alnfolder + treenode.name + '_inter.3di.fasta' , verbose=verbose)
                    except:
                        print( treenode , childalnsAA , childalns3di , bridgesAA , bridges3di)
                        raise Exception('merge error 1')
                    
                else:

                    try:
                        alnAA = mergealns( alnAA, bridgeAA , alnfolder + treenode.name + '_inter.fasta' , verbose=verbose)
                        aln3di = mergealns( aln3di, bridge3di , alnfolder + treenode.name + '_inter.3di.fasta' , verbose=verbose)

                        alnAA = mergealns( childalnsAA[c1], alnAA , alnfolder + treenode.name + '_inter.fasta' , verbose=verbose)
                        aln3di = mergealns( childalns3di[c1], aln3di , alnfolder + treenode.name + '_inter.3di.fasta' , verbose=verbose)

                        alnAA = mergealns( childalnsAA[c2], alnAA , alnfolder + treenode.name + '_inter.fasta' , verbose=verbose)
                        aln3di = mergealns( childalns3di[c2], aln3di , alnfolder + treenode.name + '_inter.3di.fasta' , verbose=verbose)
                    except:
                        print( treenode , childalnsAA , childalns3di , bridgesAA , bridges3di)
                        raise Exception('merge error 2') 
            treenode.aln = alnAA
            treenode.aln3di = aln3di
        
        if verbose == True:
            #check if node is root  
            if treenode.up == None:
                print('final aln')
                print('childalnsAA', childalnsAA)
                print('childalns3di', childalns3di)
        
    return treenode.aln, treenode.aln3di

def remove_redundant( alignment ):
    """
    this function removes redundant sequences from an alignment
    """
    aln = SeqIO.parse(alignment, 'fasta')
    seqs = []
    ids = []
    for s in aln:
        if s.id not in ids:
            seqs.append(s)
            ids.append(s.id)
    
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


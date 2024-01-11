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

def stringconvolution( s1, s2):
    #convolution of the two sequences
    #select the smaller sequence

    
    if len(s1) < len(s2):
        s = s1
        temp = s2
    else:
        s = s2
        temp = s1

    #slide the smaller sequence over the larger one
    #from a position of -len(s) to len(temp) + len(s)
    conv = []
    for i in range(-len(s), len(temp) + len(s)):
        #select the part of the larger sequence that overlaps with the smaller one
        if i < 0:
            temp2 = temp[:i+len(s)]
            s2 = s[-i:]
        elif i > len(temp):
            temp2 = temp[i:]
            s2 = s[:len(temp2)]
        else:
            temp2 = temp[i:i+len(s)]
            s2 =  s[:len(temp2)]
       
        #iterate over the smaller sequence
        score = 0
        for j,char in enumerate(s2):
            try:
                if char == temp2[j]:
                    score += 1
            except IndexError:
                print('temp2', temp2)
                print('s2', s2)
                print('i', i)
                print('j', j)
                raise IndexError
        conv.append(score)
    maxaln = np.argmax(conv)
    print('conv', conv)
    print( 'nchars' , np.amax(conv) ) 
    print('maxaln', maxaln)

    #chop the larger sequence to the length of the smaller one
    if maxaln < 0:
        temp = temp[:maxaln+len(s)]
    elif maxaln > len(temp):
        temp = temp[maxaln:]
    else:
        temp = temp[maxaln:maxaln+len(s)]
    

    #chop the smaller one if needed
    if len(temp) < len(s):
        s = s[:len(temp)]

    #find the starting position for both sequences
    if len(s1) > len(s2):
        return maxaln, temp, s
    else:
        return maxaln, s, temp

def mergealns( aln1f, aln2f, outfile):
    #find sequences in common between the two alignments
    aln1 = SeqIO.parse(aln1f, 'fasta')
    aln2 = SeqIO.parse(aln2f, 'fasta')
    ids1d = {}
    ids2d = {}
    idlist1 = []

    for s in aln1:
        if s.id not in ids1d:
            ids1d[s.id] = s
            idlist1.append(s.id)
    for s in aln2:
        if s.id not in ids2d:
            ids2d[s.id]= s
            idlist1.append(s.id)
    
    #find the intersection of the two sets
    ids1 = set(ids1d.keys())
    ids2 = set(ids2d.keys())
    commonids = ids1.intersection(ids2)
    #select one common sequence
    coidx1 = -1
    coidx2 = -1

    #find the start of the common sequence in each alignment
    #select a random common sequence
    commonid = random.choice(list(commonids))

    s1 = ids1d[commonid]
    s2 = ids2d[commonid]

    s1 = str(s1.seq).replace('-','')
    s2 = str(s2.seq).replace('-','')
    #find the start of the common sequence in each alignment
    coidx1 = s1.find(s2)
    coidx2 = s2.find(s1)

    #if the common subsequence is not found start by removing the first character of the common sequence
    if coidx1 == -1 and coidx2 == -1:
        #convolution of the two sequences
        maxaln, subs1, subs2 = stringconvolution(s1,s2)

        #find the start of the common sequence in each alignment
        if len(s1) > len(s2):
            coidx1 = maxaln
            coidx2 = len(s2) - (len(s1) - maxaln)
        else:
            coidx1 = len(s1) - (len(s2) - maxaln)
            coidx2 = maxaln
    
    print('commonid', commonid)
    print('coidx1', coidx1)
    print('coidx2', coidx2)



    #create an index of the sequence ids in each alignment
    ids1 = {}
    ids2 = {}
    aln1 = SeqIO.parse(aln1f, 'fasta')
    aln2 = SeqIO.parse(aln2f, 'fasta')
    for i,s in enumerate(aln1):
        ids1[s.id] = i
    
    for i,s in enumerate(aln2):
        ids2[s.id] = i
    

    #transform both alignments into numpy matrices
    aln1 = SeqIO.parse(aln1f, 'fasta')
    aln2 = SeqIO.parse(aln2f, 'fasta')
    aln1 = np.array([ list(str(s.seq)) for s in aln1])
    aln2 = np.array([ list(str(s.seq)) for s in aln2])

    #find the index of the common sequence in each alignment

    idx1 = ids1[commonid]
    idx2 = ids2[commonid]

    #merge the two alignments
    #remove gaps from the common sequence
    commonseq1 = aln1[idx1]
    commonseq1 = iter(commonseq1)
    commonseq2 = aln2[idx2]
    commonseq2 = iter(commonseq2)
    char1 = next(commonseq1)
    char2 = next(commonseq2)
    #iterate over the characters of the common sequence
    newaln1 = []
    newaln2 = []
    i = 0
    j = 0
    while True:
        try:
            if char1 == '-' and char2 == '-':
                #no pivot information, skip
                char1 = next(commonseq1)
                char2 = next(commonseq2)
                i += 1
                j += 1
            
            elif char2 == '-':
                char2 = next(commonseq2)
                j += 1
                #create a column of gaps
                gaps = np.array(['-'] * len(aln1))
                newaln1.append(gaps)
                newaln2.append(aln2[:,j])

            elif char1 == '-':
                char1 = next(commonseq1)
                i += 1
                #create a column of gaps
                gaps = np.array(['-'] * len(aln1))
                newaln1.append(gaps)
                newaln2.append(aln2[:,j])

            elif char1 == char2:
                char1 = next(commonseq1)
                char2 = next(commonseq2)
                i += 1
                j += 1
                newaln1.append(aln1[:,i])
                newaln2.append(aln2[:,j])


        except StopIteration:
            break
    
    newaln1 = np.array(newaln1).T
    newaln2 = np.array(newaln2).T
    #concatenate the two alignments
    newaln = np.concatenate((newaln1, newaln2), axis = 1)
    #write out the new alignment
    with open(outfile, 'w') as f:
        for i,s in enumerate(newaln):
            f.write('>' + idlist1[i] + '\n')
            f.write(''.join(s) + '\n')
    
    return outfile  

def sub2fasta( sub, outfile , fastacol1='qaln' , fastacol2='taln' ):
    with open(outfile, 'w') as f:
        f.write('>' + sub['query'] + '\n')
        f.write(sub[fastacol1] + '\n')
        f.write('>' + sub['target'] + '\n')
        f.write(sub[fastacol2] + '\n')    
    return outfile

def retalns(allvall, leafname,leafset):
    sub = allvall[allvall['query'] == leafname]
    sub = sub[sub['target'].isin(leafset)]
    sub = sub[sub['query'] != sub['target']]
    return sub.iloc[0]

#traverse tree from root to leaves recursively
def traverse_tree_merge( treenode, topleafset, allvall , alnfolder ):
    """
    this function traverses a tree from root to leaves recursively
    it returns a dictionary with the iteratively built alignment
    """
    
    if treenode.is_leaf():  
        #if the node is a leaf, then we need to add it to the alignment with one of the pivots in the current leafset
        sub = retalns(allvall, treenode.name , topleafset )
        treenode.aln = sub2fasta(sub, alnfolder + treenode.name + '_inter.fasta')
        treenode.aln3di = sub2fasta(sub, alnfolder + treenode.name + '_inter.3di.fasta' , fastacol1='3di_qaln_mode2' , fastacol2='3di_taln_mode2')
    else:
        childalns3di = []
        childalnsAA = []
        treenode.leafset = get_leafset(treenode) 
        for c in treenode.get_children():
            print('traverse', c.name , c.is_leaf() , c.leafset)
            if not c.aln:
                c.aln,c.aln3di = traverse_tree_merge(c , treenode.leafset , allvall, alnfolder)
            childalnsAA.append(c.aln)
            childalns3di.append(c.aln3di)
        treenode.aln = mergealns( childalnsAA[0], childalnsAA[1] ,  alnfolder + treenode.name + '_inter.fasta')            
        treenode.aln3di = mergealns( childalns3di[0], childalns3di[1] , alnfolder + treenode.name + '_inter.3di.fasta' )
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


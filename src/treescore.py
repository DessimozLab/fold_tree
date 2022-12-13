import toytree
#taxonomy overlap score
def getTaxOverlap(node):
	if node.is_leaf() == True:
		node.add_feature( 'score' ,  0 )
		node.add_feature( 'size' ,  0 )
		return node.lineage
	else:
		lengths = []
		total = 0
		sets = [] 
		scores = []
		for i,c in enumerate(node.get_children()):
			sets.append( getTaxOverlap(c))
			scores.append(c.score)
			total += c.score
		sets = [s for s in sets if s]
		if len(sets)> 0 :
			for i,cset in enumerate(sets):
				if i == 0:
					nset = cset
				else:
					nset = nset.intersection(cset)
				lengths.append(len(cset))
			score = len(nset) + total
			node.add_feature( 'size' ,  abs( len(nset) - max(lengths) ) )
		else:
			nset = None 
			node.add_feature( 'size' ,  0 )
		node.add_feature( 'lineage' ,  nset )
		#only in the case of a leaf with no label
		node.add_feature( 'score' ,  total )
	return nset 


def make_lineages(uniprot_df):
	return dict(zip(uniprot_df['query'] 
		, uniprot_df['Taxonomic lineage (Ids)'].map( lambda x : set( x.split(',') ) ) ))

def label_leaves( tree , leaf_lineages):
	#takes a pandas dataframe with lineage info from uniport
	for n in tree.treenode.iter_leaves():
		if n.name in leaf_lineages:
			n.add_feature( 'lineage' ,   leaf_lineages[n.name] )
		else:
			n.add_feature( 'lineage' ,   None )
	return tree


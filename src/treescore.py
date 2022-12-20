import toytree
import toyplot

#taxonomy overlap score

def standard_treedraw( tre, sizes= None , colors= None ,fixed_order=None, fixed_position=None , ts = None,  save_file = False  , tiplabels = None):
	if tiplabels is None:
		tiplabels = tre.get_tip_labels()
	canvas, axes, mark = tre.draw(  
		ts = ts,
		node_sizes=sizes,
		node_colors=colors,
		tip_labels_align=True,
		scalebar=True,
		fixed_order=fixed_order, 
		fixed_position=fixed_position,
		tip_labels=tiplabels,
		tip_labels_style={
			"fill": "#262626",
			"font-size": "9px"}
	)
	if save_file:
		toyplot.svg.render(canvas, save_file)


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
			total += c.score
		sets = [s for s in sets if s]
		if len(sets)> 0:
			for i,cset in enumerate(sets):
				if i == 0:
					nset = cset
				else:
					nset = nset.intersection(cset)
				lengths.append(len(cset))
			score = len(nset) + total
			node.add_feature( 'score' ,  score )
			#show the biggest loss in lineage length
			node.add_feature( 'size' ,  abs( len(nset) - max(lengths) ) )
		else:
			nset = None
			
			node.add_feature( 'size' ,  0 )
			node.add_feature( 'score' ,  0 )
		node.add_feature( 'lineage' ,  nset )
		#only in the case of a leaf with no label
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


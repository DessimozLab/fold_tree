import toytree
import toyplot

#taxonomy overlap score

def standard_treedraw( tre, sizes= None , colors= None ,fixed_order=None, fixed_position=None , 
ts = None,  save_file = False  , tiplabels = None ,  layout='c', edge_type='p' ):
	
	'''
	standard treedraw function
	
	Parameters
	----------
	tre : toytree.tree
		tree to draw
	sizes : list
		list of node sizes
	colors : list
		list of node colors
	fixed_order : list
		list of node names in order to draw
	fixed_position : list
		list of node positions
	ts : treestlye : str
		style of tree to draw (default = 'c')
	save_file : str
		path to save file
	tiplabels : list
		list of tip labels to draw
	layout : str
		layout of tree to draw (default = 'c')
	edge_type : str
		type of edge to draw (default = 'p')
	'''

	
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
		layout=layout,
		edge_type=edge_type,
		tip_labels_style={
			"fill": "#262626",
			"font-size": "9px"}
	)
	if save_file:
		toyplot.svg.render(canvas, save_file)



def exp_score( v , exp = 1.5):
	return v**exp

def frac_score( v ):
	return v+1

#taxonomy overlap score
def getTaxOverlap(node , treelen  = None , scorefun = frac_score):
	
	"""
    Calculate the taxonomy overlap score for the given node in a phylogenetic tree.
    
    The taxonomy overlap score is defined as the number of taxonomic labels shared by all the leaf nodes
    descended from the given node, plus the sum of the scores of all its children. If a leaf node has no
    taxonomic label, it is not counted towards the score. The function also calculates the size of the
    largest loss in lineage length, defined as the difference between the length of the set of taxonomic
    labels shared by all the leaf nodes and the length of the longest set of taxonomic labels among the
    children of the node.
    
    The function adds the following features to the node object:
    - 'score': the taxonomy overlap score.
    - 'size': the largest loss in lineage length.
    - 'lineage': the set of taxonomic labels shared by all the leaf nodes descended from the node.
    
    Parameters:
    node (Toytree.): The node in a phylogenetic tree.
    
    Returns:
    set: The set of taxonomic labels shared by all the leaf nodes descended from the node, or `None` if
    the node has no children with taxonomic labels.
    """

	if node.is_leaf() == True:
		node.add_feature( 'score' ,  0 )
		node.add_feature( 'size' ,  0 )
		node.add_feature( 'score_x_frac' , 0)
		return node.lineage
	else:
		lengths = []
		total = 0
		redtotal = 0
		fractotal = 0
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
					nset = cset.intersection(nset)
				lengths.append(len(cset))
			#add the number of unique lineages
			score = len(nset) + total
			#add the number of unique lineages weighted by the fraction of the tree
			
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

#taxonomy overlap score
def getTaxOverlap_root(node , leaf_lineages = None):
	
	"""
    Calculate the taxonomy overlap score from the root down for the given node in a phylogenetic tree.
    
	start with the total set of all clades from leaves
	use the sets from the leaf to root approach and accumlate score as the total number 
	of shared elements or frac of shared elements

    The function adds the following features to the node object:
    - 'root_score': the taxonomy overlap score.

    Parameters:
    node (Toytree.): The node in a phylogenetic tree.
    
    Returns:
    set: The set of taxonomic labels shared by all the leaf nodes descended from the node, or `None` if
    the node has no children with taxonomic labels.
    """
	if node.is_root() == True:
		leaf_lineages = [ n.lineage for n in node.get_leaves()]
		leaf_lineages = [ l for l in leaf_lineages if l ]
		leaf_lineages = [ item for sublist in leaf_lineages for item in sublist ]
		leaf_lineages = set(leaf_lineages)
		if node.lineage:
			node.add_feature( 'root_score' ,  len(node.lineage) )
		else:
			node.add_feature( 'root_score' ,  0 )
		for i,c in enumerate(node.get_children()):
			getTaxOverlap_root(c , leaf_lineages = leaf_lineages)
	else:
		if node.lineage:
			total = node.up.root_score + len(node.lineage)
		else:
			total = node.up.root_score
		node.add_feature( 'root_score' ,  total )
		for i,c in enumerate(node.get_children()):
			getTaxOverlap_root(c ,  leaf_lineages = leaf_lineages)


def sum_rootscore(node):
	return sum([n.root_score for n in node.get_leaves()])


def make_lineages(uniprot_df):
	return dict(zip(uniprot_df['query'] 
		, uniprot_df['Taxonomic lineage (Ids)'].map( lambda x :  set( x.split(',') ) ) ) )

def label_leaves( tree , leaf_lineages):
	"""
    Adds lineage information to the leaves of a tree.
    
    Parameters:
    tree (toytree.tree.TreeNode): A tree object from the toytree package.
    leaf_lineages (dict): A dictionary mapping leaf names to lineage information.
    
    Returns:
    toytree.tree.TreeNode: The input tree object with the added lineage information.
    
    Examples:
    >>> tree = toytree.tree('''((a, b), c);''')
    >>> leaf_lineages = {'a': 'Eukaryota', 'b': 'Eukaryota'}
    >>> label_leaves(tree, leaf_lineages)
    toytree.tree.TreeNode
    """
	#takes a pandas dataframe with lineage info from uniprot
	for n in tree.treenode.iter_leaves():
		if n.name in leaf_lineages:
			n.add_feature( 'lineage' ,   leaf_lineages[n.name] )
		else:
			n.add_feature( 'lineage' ,   None )
	return tree

def compute_sum_dist_to_desc_leaves(node, sum_d = 0, n_leaves = 0, n_internal_nodes = 0):
    """
    Compute the sum of distances to the descendant leaves and the number of descendant leaves for each node.
    
	"""

    n_internal_nodes = 0
    if node.is_leaf():
        n_leaves = 1
        sum_d = node.get_distance(node.up)
        node.n_desc_leaves = 0
        node.sum_dist_to_desc_leaves = 0
        node.n_internal_nodes = 0
    else:
        n_leaves = 0
        sum_d = 0
        for child in node.children:
            res = compute_sum_dist_to_desc_leaves(child, sum_d, n_leaves, n_internal_nodes)
            sum_d += res[0] 
            n_leaves += res[1]
            n_internal_nodes += res[2]
        n_internal_nodes += 1
        node.sum_dist_to_desc_leaves = sum_d
        node.n_desc_leaves = n_leaves
        node.n_internal_nodes = n_internal_nodes
        if node.up:
            sum_d += node.get_distance(node.up) * n_leaves
    if node.up:
        return sum_d, n_leaves, n_internal_nodes
    
def compute_red_score(node, red = 0, level_from_root = 0):
    """
    Compute the RED score for each node.
    
    """
    if node.is_leaf():
        red = 1
    elif not node.up:
        red = 0
    else:
        p = red
        d = node.get_distance(node.up)
        u = (node.sum_dist_to_desc_leaves + (d * node.n_desc_leaves)) / node.n_desc_leaves
        red = p + (d/u) * (1-p)
    node.level = level_from_root
    node.red = red
    for child in node.children:
        compute_red_score(child, red, level_from_root + 1)
        
def labelwRED(tree):
    compute_sum_dist_to_desc_leaves(tree)
    compute_red_score(tree)
    return tree
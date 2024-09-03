import toytree
import toyplot
import re
from scipy.stats import describe
import pandas as pd
from ete3 import NCBITaxa
####### Astral functions ########

ncbi = NCBITaxa()

def parse_astral(string , index = 0):
	"""

	q1,q2,q3: these three values show quartet support (as defined in the description of -t 1) for the main topology, the first alternative, and the second alternative, respectively.
	f1, f2, f3: these three values show the total number of quartet trees in all the gene trees that support the main topology, the first alternative, and the second alternative, respectively.
	pp1, pp2, pp3: these three show the local posterior probabilities (as defined in the description of -t 4) for the main topology, the first alternative, and the second alternative, respectively.

	"""
	#'[pp1=0.284784;pp2=0.426968;pp3=0.288248;f1=0.000000;f2=0.263158;f3=0.008772;q1=0.000000;q2=0.967742;q3=0.032258]'
	string = string.replace('[','').replace(']','')
	string = string.split(';')
	string = {s.split('=')[0]:float(s.split('=')[1]) for s in string}
	#to dataframe
	string = pd.DataFrame(string , index = [index])
	return string

def retastral_support( astral_file):
	with open(astral_file , 'r') as f:
		nstring = f.read()
	#find all quoted strings with single quotes
	quoted = re.findall(r"'(.*?)'", nstring)
	dfs = [ parse_astral(q,i) for i,q in enumerate(quoted)]
	dfs = pd.concat(dfs)
	#make sure pp1 does not contain nan
	return dfs

def verify_map(mapping, tree):
	t = ete3.Tree(tree)
	for l in t.get_leaves():
		if l.name not in mapping.keys():
			print('missing', l.name)
			return False
	return True 

def return_astral_score(astrolout, logfile):
	#get dup loss and speciation events from astral logfile
	with open(logfile , 'r') as f:
		for l in f:
			if '#Duploss' in l:
				duploss = int(l.split(' ')[-1])
	astraldf = retastral_support(astrolout)
	#count nans 
	nancount = astraldf['pp1'].isna().sum()
	#get summary stats
	astral_stas = describe(astraldf['pp1'].dropna())
	return astraldf, { 'mean': astral_stas.mean , 'std': astral_stas.variance ,'skewness':astral_stas.skewness, 'nancount': nancount , 'duploss': duploss }


def run_astral( phylo , st,  mapping, output ,root  ,outfmt = '2', astralpath = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/software/ASTER-Linux/bin/astral-pro3' ):
	logfile = output+'.log'
	cmd = astralpath + ' -c {st} -a {mapping} -u {outfmt} -i {phylo} -o {output} -C -T -E --root {root} 2>{logfile}'.format(phylo = phylo , st=st, mapping = mapping , output = output , root=root, logfile = logfile , outfmt = outfmt)
	print(cmd)
	subprocess.run(cmd , shell = True)
	return output , logfile

def prepare_astral_input(uniprot_df , speciestreeout, mapperout):
	finalset = pd.read_csv(uniprot_df)
	finalset['species'] = finalset['Taxonomic lineage (Ids)'].map(lambda x: x.split(',')[-1].split('(')[0].strip())
	finalset['species'] = finalset['species'].map(lambda x: x.split('_')[0])
	mapper = dict(zip(finalset['Entry'], finalset['species']))
	with open(mapperout , 'w') as f:
		for k,v in mapper.items():
			f.write('{}\t{}\n'.format(k,int(v)))
	
	#get ncbi tree of species set
	species_set = list(finalset.species.unique())
	species_set = [int(s) for s in species_set]
	st = ncbi.get_topology(species_set, intermediate_nodes=False)
	st.name = 'root'
	#write with internal node names

	st.write(outfile=speciestreeout , format= 1 )
	return mapper , uniprot_df.replace('.csv','_speciesmap.txt' )  , uniprot_df.replace('.csv','_ncbi_tree.nwk' ) 



######## TCS functions ########


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

def node_degree(node ):
	""" add a degree feature to the node object with is the distance to the root"""
	if node.is_root() == True:
		node.add_feature( 'degree' ,  0 )
		for i,c in enumerate(node.get_children()):
			node_degree(c )
	else:
		node.add_feature( 'degree' ,  node.up.degree + 1 )
		for i,c in enumerate(node.get_children()):
			node_degree(c )
	return node
	
def node_degree_inv(node , max_degree = 0 , exp = 1.5):
	#get the max degree of all the nodes
	if node.is_root() == True:
		max_degree = max([n.degree for n in node.get_leaves()])
		print('max degree' , max_degree)
		if node.lineage:
			node.add_feature( 'inv_degree' , len(node.lineage) * 1 )
		else:
			node.add_feature( 'inv_degree' , 0 )
			
		for i,c in enumerate(node.get_children()):
			node_degree_inv(c  , max_degree= max_degree , exp= exp)
	else:
		if node.lineage:
			node.add_feature( 'inv_degree' ,  len(node.lineage)* (( max_degree - node.degree) / max_degree) **exp )
		else:
			node.add_feature( 'inv_degree' ,  0 )
		for i,c in enumerate(node.get_children()):
			node_degree_inv(c , max_degree= max_degree , exp = exp )
	return sum([n.inv_degree for n in node.traverse() if n.is_leaf() == False ])

def degree_score(node , exp):
	node = node_degree(node)
	dgscore = node_degree_inv(node , exp )
	print( 'degree score' , dgscore	 , 'exp' , exp)
	return dgscore 

#get weighted score
def lineage_score(node):
	clades = {}
	for c in node.traverse():
		if c.is_leaf() and c.lineage:
			for l in c.lineage:
				if l in clades:
					clades[l] += 1
				else:
					clades[l] = 1	
	#add the weighted score
	print('clades' , clades)
	score = sum([ sum ( [ clades[l]  for l in n.lineage] )  for n in node.traverse() if n.is_leaf()	== False and n.lineage ] )
	print('lineage score' , score)
	return score

#get weighted score
def lineage_score_woutredundant(node):
	clades = {}
	for c in node.traverse():
		if c.is_leaf() and c.lineage:
			for l in c.lineage:
				if l in clades:
					clades[l] += 1
				else:
					clades[l] = 1	
	#zero the entries that are in all leaves
	for k,v in clades.items():
		if v == len([n for n in node.get_leaves() if n.lineage]):
			clades[k] = 0

	#add the weighted score
	print('clades' , clades)
	score = sum([ sum ( [ clades[l]  for l in n.lineage] )  for n in node.traverse() if n.is_leaf()	== False and n.lineage ] )
	print('lineage score' , score)
	return score

#get weighted score
def lineage_score_woutredundant_IC(node):
	clades = {}
	for c in node.traverse():
		if c.is_leaf() and c.lineage:
			for l in c.lineage:
				if l in clades:
					clades[l] += 1
				else:
					clades[l] = 1
	#zero the entries that are in all leaves
	nleaves = len([n for n in node.get_leaves() if n.lineage])
	for k,v in clades.items():
		if v == len([n for n in node.get_leaves() if n.lineage]):
			clades[k] = 0
		else:
			#weigh by IC
			p = v/nleaves
			clades[k] = p * np.log2(p)
	#add the weighted score
	print('clades' , clades)
	score = sum([ sum ( [ clades[l]  for l in n.lineage] )  for n in node.traverse() if n.is_leaf()	== False and n.lineage ] )
	print('lineage score' , score)
	return score

def lineage_score_tax_degree(node,uniprot_df):
	clades = {}
	for c in node.traverse():
		if c.is_leaf() and c.lineage:
			for l in c.lineage:
				if l in clades:
					clades[l] += 1
				else:
					clades[l] = 1	
	#zero the entries that are in all leaves
	for k,v in clades.items():
		if v == len([n for n in node.get_leaves() if n.lineage]):
			clades[k] = 0
	taxa_degree = {}
	for lineage in list(uniprot_df['Taxonomic lineage (Ids)'] ):
		for i,tax in enumerate(lineage.split(',')):
			if tax in taxa_degree:
				taxa_degree[tax] = i
			elif tax in taxa_degree and taxa_degree[tax] > i:
				taxa_degree[tax] = i
			else:
				taxa_degree[tax] = i
	#get the max degree of all the taxa
	max_degree = max([taxa_degree[tax] for tax in taxa_degree])
	for tax in taxa_degree:
		taxa_degree[tax] = max_degree - taxa_degree[tax] + 1

	print('taxa degree' , taxa_degree)
	
	#add the weighted score
	print('clades' , clades)
	score = sum([ sum ( [ clades[l]*taxa_degree[l]  for l in n.lineage] )  for n in node.traverse() if n.is_leaf()	== False and n.lineage ] )
	print('lineage score' , score)
	return score


#taxonomy overlap score
def getTaxOverlap_root(node , clades = None):
	
	"""
	Calculate the taxonomy overlap score from the root down for the given node in a phylogenetic tree.
	
	start with the total set of all clades from leaves
	use the sets from the leaf to root approach and accumlate score as the total number 
	of shared elements

	The function adds the following features to the node object:
	- 'root_score': the taxonomy overlap score.

	Parameters:
	node (Toytree.): The node in a phylogenetic tree.
	
	Returns:
	set: The set of taxonomic labels shared by all the leaf nodes descended from the node, or `None` if
	the node has no children with taxonomic labels.
	"""
	if node.is_root() == True:
		clades = {}

		#weight of clades set to 1
		for c in node.traverse():
			if c.is_leaf() and c.lineage:
				for l in c.lineage:
					if l in clades:
						clades[l] += 1
					else:
						clades[l] = 1
		
		#zero the entries that are in all leaves			
		for k,v in clades.items():
			if v == len([n for n in node.get_leaves() if n.lineage]):
				clades[k] = 0
			else:
				clades[k] = 1
			
		if node.lineage:
			node.add_feature( 'root_score' ,  len(node.lineage) )
			node.add_feature( 'root_score_nr' ,  sum([ clades[l] for l in node.lineage]) )
		else:
			node.add_feature( 'root_score' ,  0 )
			node.add_feature( 'root_score_nr' ,  0 )

		for i,c in enumerate(node.get_children()):
			getTaxOverlap_root(c , clades = clades )
	else:
		if node.lineage:
			if node.is_leaf():
				total = node.up.root_score
			total = node.up.root_score + len(node.lineage)
			total_nr = node.up.root_score_nr + sum([ clades[l] for l in node.lineage])
		else:
			total = node.up.root_score
			total_nr = node.up.root_score_nr

		node.add_feature( 'root_score' ,  total )
		node.add_feature( 'root_score_nr' ,  total_nr )

		for i,c in enumerate(node.get_children()):
			getTaxOverlap_root(c ,  clades = clades )

def sum_rootscore(node):
	#find the taxa shared between all leaves
	print('root scoring')
	rootscore = sum([n.root_score for n in node.get_leaves()])
	print('root score' , rootscore)
	rootscore_nr = sum([n.root_score_nr for n in node.get_leaves()])	
	print('root score_nr' , rootscore_nr)

	return rootscore , rootscore_nr

def make_lineages(uniprot_df):
	return dict(zip(uniprot_df['query'] 
		, uniprot_df['Taxonomic lineage (Ids)'].map( lambda x :  set( x.split(',') ) ) ) )

def get_species(uniprot_df):
	return dict(zip(uniprot_df['query'] 
		, uniprot_df['Taxonomic lineage (Ids)'].map( lambda x :   x.split(',')[-1].split('(' )[0].strip() ) ) )

def label_leaves( tree , leaf_lineages , species_map):
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
	species_count = {}
	#takes a pandas dataframe with lineage info from uniprot
	for n in tree.treenode.iter_leaves():
		if n.name in leaf_lineages:
			n.add_feature( 'lineage' ,   leaf_lineages[n.name] )
			if species_map[n.name] in species_count:
				species_count[species_map[n.name]] += 1
			else:
				species_count[species_map[n.name]] = 1
			spcount = species_count[species_map[n.name]]
			#pad so that it is a 3 digit number
			spcount = str(spcount).zfill(3)
			n.add_feature( 'sp_num' ,   species_map[n.name]+ '_' + spcount)
		else:
			n.add_feature( 'lineage' ,   None )
			n.add_feature( 'sp_num' ,   None )
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



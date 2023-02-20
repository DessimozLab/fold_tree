import itertools
import dendropy
import sys

####small parsimony functions ##########
def process_node_smallpars_1(node,allowed_symbols , transition_dict ):

    """
    This function takes in a node, a set of allowed symbols, and a transition dictionary as input. It starts from the leaves of the node and generates a set of symbols for each node as it traverses up the tree.
    The function uses the following steps:

    If the node's symbol attribute is None, it recursively calls itself for each child of the node.
    It initializes the symbol and scores attributes for the node to be an empty set and dictionary, respectively.
    It finds the intersection of symbols for all child nodes and assigns it to the current node's symbols attribute. If the intersection is empty, it finds the union of symbols for all child nodes.
    It then iterates over the allowed symbols and checks if the current symbol is in the node's symbols. If it is not, it adds 1 to the minimum score of the current symbol from the child nodes and assigns it to the current node's scores attribute for that symbol. If the symbol is in the node's symbols, it assigns the minimum score of the current symbol from the child nodes to the current node's scores attribute for that symbol.
    Args:
    node: a node object containing child nodes and attributes for symbols and scores.
    allowed_symbols: a set of symbols that are allowed in the tree.
    transition_dict : a dictionary containing transition scores, if any, to be used in scoring the tree.

    Returns:
    None
    """
    #go from leaves up and generate character sets
    if node.symbols is None:
        for child in node.child_nodes():
            if child.symbols is None:
                process_node_smallpars_1(child, allowed_symbols, transition_dict)
        node.symbols = { }
        node.scores = { }
        symbols = set.intersection( * [ child.symbols for child in node.child_nodes( ) ] )
        if len(symbols) == 0:
            symbols = set.union( * [ child.symbols for child in node.child_nodes( ) ] )
        node.symbols = symbols
        for c in allowed_symbols:
            if c not in node.symbols:
                #add trnasition mat here if needed for more subtle scoring
                score = min(  [ child.scores[c] for child in node.child_nodes()])+1
            else:
                score = min(  [ child.scores[c] for child in node.child_nodes() ] )
            node.scores[c] = score

def process_node_smallpars_2(node , allowed_symbols , transition_dict , verbose = False):

    """
    This function takes in a node, a set of allowed symbols, a transition dictionary and a verbose flag as input. It assigns the most parsimonious character to each node as it traverses down the tree.
    The function uses the following steps:

    If the node's char attribute is None, it checks if the node has a parent.
    If the node has a parent, it initializes the char, event and eventype attributes for the node to be an empty dictionary. It then assigns the minimum score of the symbol from the node's scores attribute to the node's char attribute. It then compares the node's char attribute with the parent's char attribute. If they are the same, it assigns 0 to the node's event attribute, else it assigns 1. It also assigns the value of the transition from the parent's char to the node's char in the transition dictionary to the node's eventype attribute.
    If the node is the root node, it initializes the char, event and eventype attributes for the node to be an empty dictionary. It then assigns the minimum score of the symbol from the node's scores attribute to the node's char attribute and assigns 0 to the node's event attribute.
    It then recursively calls itself for each child of the node, if the child's char attribute is None.
    Args:
    node: a node object containing child nodes and attributes for char, event, eventype and scores.
    allowed_symbols: a set of symbols that are allowed in the tree.
    transition_dict : a dictionary containing transition scores, if any, to be used in scoring the tree.
    verbose: a flag to indicate whether to print debugging information.

    Returns:
    None
    """

    #assign the most parsimonious char from children
    if node.char is None:
        if node.parent_node:
            #node has parent
            node.char = {}
            node.event = {}
            node.eventype= {}
            node.char = min(node.scores, key=node.scores.get)
            if node.parent_node.char == node.char:
                node.event = 0
            else:
                if node.scores[node.parent_node.char] == node.scores[node.char] :
                    node.char = node.parent_node.char
                    node.event = 0
                else:
                    node.event = 1
                    node.eventype = transition_dict[(node.parent_node.char,node.char)]
        else:
            #root node
            node.char = {}
            node.event= {}
            node.eventype = {}
            node.char = min(node.scores, key=node.scores.get)
            node.event = 0
        #down one level
        for child in node.child_nodes():
            if child.char is None:
                process_node_smallpars_2(child , allowed_symbols, transition_dict)

def calculate_small_parsimony(t, allowed_symbols, transition_dict):
    """
    This function takes in a tree object (t), a set of allowed symbols, and a transition dictionary as input. It first calls the process_node_smallpars_1 function on the seed node of the tree to generate the symbol sets and scores for each node while traversing up the tree. It then calls the process_node_smallpars_2 function on the seed node of the tree to assign the most parsimonious character to each node while traversing down the tree.

    Args:
    t: a tree object containing the seed node and child nodes.
    allowed_symbols: a set of symbols that are allowed in the tree.
    transition_dict : a dictionary containing transition scores, if any, to be used in scoring the tree.

    Returns:
    t: the input tree object with the modified attributes for each node.
    """
    #up
    process_node_smallpars_1(t.seed_node , allowed_symbols, transition_dict )
    #down
    process_node_smallpars_2(t.seed_node, allowed_symbols, transition_dict )
    return t

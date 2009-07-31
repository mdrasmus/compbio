"""
 Birth death process and reconstructed process for trees

"""

from math import *
import random
from rasmus import stats, treelib


def prob_birth_death(genes1, genes2, t, birth, death):
    """Probability of 'genes1' genes at time 0 give rise to 'genes2' genes at
       time 't' with 'birth' and 'death' rates.
    """

    if birth == 0.0 and death == 0.0:
        if genes1 == genes2:
            return 1.0
        else:
            return 0.0
    
    l = birth
    u = death
    elut = exp((l-u)*t)
    a = u * (elut - 1.0) / (l*elut - u) # alpha
    b = l * (elut - 1.0) / (l*elut - u) # beta
    n = genes1
    i = genes2

    if genes1 < 1:
        return 0.0

    if genes2 == 0:
        return a ** n
    else:
        return sum(stats.choose(n,j) * stats.choose(n+i-j-1, n-1) *\
                   a**(n-j) * b**(i-j) * (1.0 - a - b)**j
                   for j in xrange(min(n, i)+1))
   

def birth_wait_time(t, n, T, birth, death):
    """Probability density for for next birth at time 't' given
       'n' lineages starting at time 0, evolving until time 'T' with a
       'birth' and 'death' rates for a reconstructed process.
    """
    
    l = birth
    u = death
    r = l - u
    a = u / l

    return n * r * exp(-n*r*t) * \
           ((1.0 - a * exp(-r * (T - t)))**(n-1)) / \
           ((1.0 - a * exp(-r * T))**n)


def prob_no_birth(n, T, birth, death):
    """Probability of no birth from 'n' lineages starting at time 0, 
       evolving until time 'T' with 'birth' and 'death' rates
       for a reconstructed process.
    """

    if birth == 0.0:
        return 1.0
    
    l = birth
    u = death
    r = l - u
    a = u / l

    return (1.0 - (l*(1.0 - exp(-r * T))) / \
                  (l - u * exp(-r * T))) ** n


def sample_birth_wait_time(n, T, birth, death):
    """
    Sample the next birth event from a reconstructed birthdeath process.
    Let there be 'n' lineages at time 0 that evolve until time 'T' with
    'birth' and 'death' rates.

    Conditioned that a birth will occur
    """
    
    # TODO: could make this much more efficient (use straight line instead of
    # flat line).
    
    # uses rejection sampling
    start_y = birth_wait_time(0, n, T, birth, death)
    end_y = birth_wait_time(T, n, T, birth, death)
    M = max(start_y, end_y)
    
    while True:
        t = random.uniform(0, T)
        f = birth_wait_time(t, n, T, birth, death)

        if random.uniform(0, 1) <= f / M:
            return t


def sample_birth_death_tree(T, birth, death, tree=None, node=None):
    """Simulate a reconstructed birth death tree"""
    
    # create tree if one is not given
    if tree == None:
        tree = treelib.Tree()
    
    # create starting node if one is not given
    if node == None:
        tree.make_root() 
        node = tree.root
    else:
        node2 = treelib.TreeNode(tree.new_name())
        tree.add_child(node, node2)
        node = node2
    
    def walk(tn, node):
        n = 1
        
        # determine if this is the last branch
        if random.uniform(0, 1) < prob_no_birth(n, T-tn, birth, death):
            node.dist = T - tn
        else:
            t = sample_birth_wait_time(n, T-tn, birth, death)
            node.dist = t
            
            assert tn + t < T
            
            # recurse
            node2 = treelib.TreeNode(tree.new_name())
            tree.add_child(node, node2)
            walk(tn + t, node2)

            node2 = treelib.TreeNode(tree.new_name())
            tree.add_child(node, node2)
            walk(tn + t, node2)
    walk(0.0, node)

    return tree


def sample_birth_death_gene_tree(stree, birth, death, 
                                 genename=lambda sp, x: sp + "_" + str(x),
                                 removeloss=True):
    """Simulate a gene tree within a species tree with birth and death rates"""
    
    # initialize gene tree
    tree = treelib.Tree()
    tree.make_root()
    recon = {tree.root: stree.root}
    events = {tree.root: "spec"}
    losses = set()
    
    def walk(snode, node):
        if snode.is_leaf():
            tree.rename(node.name, genename(snode.name, node.name))
            events[node] = "gene"
        else:
            for child in snode:
                # determine if loss will occur
                if random.uniform(0, 1) < prob_birth_death(1, 0, child.dist, 
                                                           birth, death):
                    continue
                
                # no loss, so simulate reconstructed tree
                sample_birth_death_tree(child.dist, birth, death, 
                                        tree=tree, node=node)
                
                # record reconciliation
                next_nodes = []
                def walk2(node):
                    node.recurse(walk2)
                    recon[node] = child
                    if node.isLeaf():
                        events[node] = "spec"
                        next_nodes.append(node)
                    else:
                        events[node] = "dup"
                walk2(node.children[-1])
                
                # recurse
                for leaf in next_nodes:
                    walk(child, leaf)
            
            # if no child for node then it is a loss
            if node.is_leaf():
                losses.add(node)
    walk(stree.root, tree.root)
    
    
    # remove lost nodes
    if removeloss:
        treelib.remove_exposed_internal_nodes(tree,
                                              set(tree.leaves()) - losses)
        treelib.remove_single_children(tree)
        
        delnodes = set()
        for node in recon:
            if node.name not in tree.nodes:
                delnodes.add(node)
        for node in delnodes:
            del recon[node]
            del events[node]
    
    return tree, recon, events



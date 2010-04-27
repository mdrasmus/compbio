"""

   Code for the DLCoal model
   (duplications, losses, and coalescence)

"""

import random
from itertools import chain, izip

from rasmus import treelib

from . import birthdeath
from . import coal



def sample_dlcoal(stree, n, duprate, lossrate, namefunc=None):
    """Sample a gene tree from the DLCoal model"""

    # generate the locus tree
    locus_tree, locus_recon, locus_events = \
                birthdeath.sample_birth_death_gene_tree(
        stree, duprate, lossrate)

    if len(locus_tree.nodes) <= 1:
        # total extinction
        coal_tree = treelib.Tree()
        coal_tree.make_root()
        coal_recon = {coal_tree.root: locus_tree.root}
        daughters = set()
    else:
        # simulate coalescence
        
        # choose daughter duplications
        daughters = set()
        for node in locus_tree:
            if locus_events[node] == "dup":
                daughters.add(node.children[random.randint(0, 1)])

        try:
            coal_tree, coal_recon = sample_multicoal_tree(locus_tree, n,
                                                      daughters=daughters,
                                                      namefunc=namefunc)
        except:
            print locus_tree.nodes
            treelib.draw_tree_names(locus_tree, maxlen=5)
            raise

    # store extra information
    extra = {"locus_tree": locus_tree,
             "locus_recon": locus_recon,
             "locus_events": locus_events,
             "coal_tree": coal_tree,
             "coal_recon": coal_recon,
             "daughters": daughters}

    return coal_tree, extra





def sample_multicoal_tree(stree, n, leaf_counts=None,
                          daughters=set(),
                          namefunc=None):
    """
    Returns a gene tree from a multi-species coalescence process
    n -- population size (int or dict)
         If n is a dict it must map from species name to population size
    """

    # initialize vector for how many genes per extant species
    if leaf_counts is None:
        leaf_counts = dict((l, 1) for l in stree.leaf_names())

    # initialize function for generating new gene names
    if namefunc is None:
        spcounts = dict((l, 1) for l in stree.leaf_names())
        def namefunc(sp):
            name = sp + "_" + str(spcounts[sp])
            spcounts[sp] += 1
            return name

    # initialize population sizes
    if isinstance(n, (int, float)):
        popsizes = dict.fromkeys(stree.nodes.keys(), n)
    elif isinstance(n, dict):
        popsizes = n
    else:
        raise Exception("n must be a int or dict.")

    # init gene counts
    counts = dict((n.name, 0) for n in stree)
    counts.update(leaf_counts)

    # init reconciliation, events
    recon = {}

    subtrees = {}

    # loop through species tree
    for snode in stree.postorder():        
        # simulate population for one branch
        k = counts[snode.name]
        if snode in daughters:
            # daughter branch, use bounded coalescent
            subtree = coal.sample_coal_tree_bounded(
                k, popsizes[snode.name], snode.dist, capped=True)
            lineages = set(subtree.root)
        elif snode.parent:            
            # non basal branch
            subtree, lineages = coal.sample_coal_tree_fixed(
                k, popsizes[snode.name], snode.dist, capped=True)
        else:
            # basal branch
            subtree = coal.sample_coal_tree(k, popsizes[snode.name])
            lineages = set(subtree.root)
        subtrees[snode] = (subtree, lineages)
        if snode.parent:
            counts[snode.parent.name] += len(lineages)
        for node in subtree:
            recon[node] = snode


    # stitch subtrees together
    tree = treelib.Tree()

    # add all nodes to total tree
    for subtree, lineages in subtrees.values():
        tree.merge_names(subtree)
        tree.remove(subtree.root)    
    
    for snode in stree:
        if not snode.is_leaf():
            subtree, lineages = subtrees[snode]

            # get lineages from child subtrees
            lineages2 = chain(*[subtrees[child][1]
                                for child in snode.children])

            # ensure leaves are randomly attached
            leaves = subtree.leaves()
            random.shuffle(leaves)

            # stitch leaves of the subtree to children subtree lineages
            for leaf, lineage in izip(leaves, lineages2):
                tree.add_child(leaf, lineage)


    # set root
    tree.root = subtrees[stree.root][0].root
    tree.add(tree.root)

    # name leaves
    for leaf in tree.leaves():
        tree.rename(leaf.name, namefunc(recon[leaf].name))
        
    return tree, recon



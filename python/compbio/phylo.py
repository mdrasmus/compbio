#
# Phylogeny functions
# Matt Rasmussen 2006-2010
# 


# python imports
import math
import os
import random
import sys


# rasmus imports
from rasmus import tablelib
from rasmus import stats
from rasmus import treelib
from rasmus import util

# compbio imports
from . import fasta


#=============================================================================
# Counting functions

def num_rooted_trees(nleaves):
    # (2n-3)!! = (2n-3)!/[2^(n-2)*(n-2)!] for n >= 2
    assert nleaves >= 1
    if nleaves < 2:
        return 1
    return stats.factorial(2*nleaves-3)/(2**(nleaves-2) * stats.factorial(nleaves-2))

def num_unrooted_trees(nleaves):
    # (2n-5)!! = (2n-5)!/[2^(n-3)*(n-3)!] for n >= 3
    assert nleaves >= 1
    if nleaves < 3:
        return 1
    return stats.factorial(2*nleaves-5)/(2**(nleaves-3) * stats.factorial(nleaves-3))


#=============================================================================
# Gene to species mapping functions


def gene2species(genename):
    """Default gene2species mapping"""
    return genename


def make_gene2species(maps):
    """Returns a function that maps gene names to species names

    maps -- a list of tuples [(gene_pattern, species_name), ... ]
    """
    
    # find exact matches and expressions
    exacts = {}
    exps = []
    for mapping in maps:
        if "*" not in mapping[0]:
            exacts[mapping[0]] = mapping[1]
        else:
            exps.append(mapping)
    
    # create mapping function
    def gene2species(gene):
        # eval expressions first in order of appearance
        for exp, species in exps:
            if exp[-1] == "*":
                if gene.startswith(exp[:-1]):
                    return species
            elif exp[0] == "*":
                if gene.endswith(exp[1:]):
                    return species
        
        if gene in exacts:
            return exacts[gene]
        
        raise Exception("Cannot map gene '%s' to any species" % gene)
    return gene2species


def read_gene2species(* filenames):
    """
    Reads a gene2species file

    Returns a function that will map gene names to species names.
    """
    
    for filename in filenames:
        maps = []
        for filename in filenames:
            maps.extend(util.read_delim(util.skip_comments(
                util.open_stream(filename))))
    return make_gene2species(maps)


#=============================================================================
# Reconciliation functions
#
    

def reconcile(gtree, stree, gene2species=gene2species):
    """
    Returns a reconciliation dict for a gene tree 'gtree' and species tree 'stree'
    """
    
    recon = {}
    
    # determine the preorder traversal of the stree
    order = {}
    def walk(node):
        order[node] = len(order)
        node.recurse(walk)
    walk(stree.root)

    
    # label gene leaves with their species
    for node in gtree.leaves():
        recon[node] = stree.nodes[gene2species(node.name)]
    
    # recurse through gene tree
    def walk(node):
        node.recurse(walk)
        
        if not node.is_leaf():
            # this node's species is lca of children species  
            recon[node] = reconcile_lca(stree, order, 
                                       util.mget(recon, node.children))
    walk(gtree.root)
    
    return recon


def reconcile_lca(stree, order, nodes):
    """Helper function for reconcile"""
    
    # handle simple and complex cases
    if len(nodes) == 1:
        return nodes[0]    
    if len(nodes) > 2:
        return treelib.lca(nodes)
    
    # 2 node case
    node1, node2 = nodes
    index1 = order[node1]
    index2 = order[node2]
    
    while index1 != index2:
        if index1 > index2:
            node1 = node1.parent
            index1 = order[node1]
        else:
            node2 = node2.parent
            index2 = order[node2]
    return node1
    

def reconcile_node(node, stree, recon):
    """Reconcile a single gene node to a species node"""
    return treelib.lca([recon[x] for x in node.children])


def label_events(gtree, recon):
    """Returns a dict with gene node keys and values indicating 
       'gene', 'spec', or 'dup'"""
    events = {}
    
    def walk(node):
        events[node] = label_events_node(node, recon)
        node.recurse(walk)
    walk(gtree.root)
    
    return events


def label_events_node(node, recon):
    if not node.is_leaf():
        if recon[node] in map(lambda x: recon[x], node.children):
            return "dup"
        else:
            return "spec"
    else:
        return "gene"


def find_loss_node(node, recon):
    """Finds the loss events for a branch in a reconciled gene tree"""
    loss = []
    
    # if not parent, then no losses
    if not node.parent:
        return loss
    
    # determine starting and ending species
    sstart = recon[node]
    send = recon[node.parent]
    
    # determine species path of this gene branch (node, node.parent)
    ptr = sstart
    spath = []
    while ptr != send:
        spath.append(ptr)
        ptr = ptr.parent
    
    # determine whether node.parent is a dup
    # if so, send (species end) is part of species path
    if label_events_node(node.parent, recon) == "dup":
        spath.append(send)
    
    # go up species path (skip starting species)
    # every node on the list is at least one loss
    for i, snode in enumerate(spath[1:]):
        for schild in snode.children:
            if schild != spath[i]:
                loss.append([node, schild])
        
    return loss


def find_loss_under_node(node, recon):
    loss = []
    snodes = {}
    internal = {}
    species1 = recon[node]

    # walk from child species to parent species
    for child in node.children:
        ptr = recon[child]
        snodes[ptr] = 1
        while ptr != species1:
            ptr = ptr.parent
            snodes[ptr] = 1
            internal[ptr] = 1

    # foreach internal node in partial speciation tree, all children
    # not in speciation are loss events
    for i in internal:
        for child in i.children:
            if child not in snodes:
                loss.append([node,child])
    return loss


def find_loss(gtree, stree, recon, node=None):
    """Returns a list of gene losses in a gene tree"""
    loss = []

    def walk(node):
        loss.extend(find_loss_node(node, recon))
        node.recurse(walk)
    if node:
        walk(node)
    else:
        walk(gtree.root)

    return loss


def count_dup(gtree, events, node=None):
    """Returns the number of duplications in a gene tree"""
    var = {"dups": 0}
    
    def walk(node):
        if events[node] == "dup":
            var["dups"] += len(node.children) - 1
        node.recurse(walk)
    if node:
        walk(node)
    else:
        walk(gtree.root)
    
    return var["dups"]


def count_loss(gtree, stree, recon, node=None):
    """Returns the number of losses in a gene tree"""
    return len(find_loss(gtree, stree, recon, node))


def count_dup_loss(gtree, stree, recon, events=None):
    """Returns the number of duplications + losses in a gene tree"""
    if events is None:
        events = label_events(gtree, recon)
    
    nloss = count_loss(gtree, stree, recon)
    ndups = count_dup(gtree, events)
    return nloss + ndups


def find_species_roots(tree, stree, recon):
    """Find speciation nodes in the gene tree that reconcile to the
       species tree root"""
    
    roots = []
    def walk(node):
        found = False
        for child in node.children:
            found = walk(child) or found
        if not found and recon[node] == stree.root:
            roots.append(node)
            found = True
        return found
    walk(tree.root)
    return roots       


def find_orthologs(gtree, stree, recon, counts=True):
    """Find all ortholog pairs within a gene tree"""

    events = label_events(gtree, recon)
    orths = []
    
    for node, event in events.items():
        if event == "spec":
            leavesmat = [x.leaves() for x in node.children]
            sp_counts = [util.hist_dict(util.mget(recon, row))
                         for row in leavesmat]
            
            for i in range(len(leavesmat)):
                for j in range(i+1, len(leavesmat)):
                    for gene1 in leavesmat[i]:
                        for gene2 in leavesmat[j]:
                            if gene1.name > gene2.name:
                                g1, g2 = gene2, gene1
                                a, b = j, i
                            else:
                                g1, g2 = gene1, gene2
                                a, b = i, j
                            
                            if not counts:
                                orths.append((g1.name, g2.name))
                            else:
                                orths.append((g1.name, g2.name,
                                              sp_counts[a][recon[g1]],
                                              sp_counts[b][recon[g2]]))
    
    return orths



def subset_recon(tree, recon, events=None):
    """Ensure the reconciliation only refers to nodes in tree"""

    # get all nodes that are walkable
    nodes = set(tree.postorder())
    for node in list(recon):
        if node not in nodes:
            del recon[node]
    if events:
        for node in list(events):
            if node not in nodes:
                del events[node]
        


#=============================================================================
# Reconciliation Input/Output

def write_recon(filename, recon):
    """Write a reconciliation to a file"""
    util.write_delim(filename, [(str(a.name), str(b.name))
                                for a,b in recon.items()])


def read_recon(filename, tree1, tree2):
    """Read a reconciliation from a file"""    
    recon = {}
    for a, b in util.read_delim(filename):
        if a.isdigit(): a = int(a)
        if b.isdigit(): b = int(b)
        recon[tree1.nodes[a]] = tree2.nodes[b]
    return recon


def write_events(filename, events):
    """Write events data structure to file"""
    util.write_delim(filename, [(str(a.name), b) for a,b in events.items()])


def read_events(filename, tree):
    """Read events data structure from file"""
    events = {}
    for name, event in util.read_delim(filename):
        if name.isdigit(): name = int(name)
        events[tree.nodes[name]] = event
    return events


def write_recon_events(filename, recon, events=None, noevent=""):
    """Write a reconciliation and events to a file"""
    
    if events is None:
        events = dict.fromkeys(recon.keys(), noevent)
    
    util.write_delim(filename, [(str(a.name), str(b.name), events[a])
                                for a,b in recon.items()])

def read_recon_events(filename, tree1, tree2):
    """Read a reconciliation and events data structure from file"""
    
    recon = {}
    events = {}
    for a, b, event in util.read_delim(filename):
        if a.isdigit(): a = int(a)
        if b.isdigit(): b = int(b)
        node1 = tree1.nodes[a]
        recon[node1] = tree2.nodes[b]       
        events[node1] = event
    return recon, events


#============================================================================
# duplication loss counting
#


def init_dup_loss_tree(stree):
    # initalize counts to zero
    def walk(node):
        node.data['dup'] = 0
        node.data['loss'] = 0
        node.data['appear'] = 0
        node.data['genes'] = 0
        node.recurse(walk)
    walk(stree.root)


def count_dup_loss_tree(tree, stree, gene2species, recon=None):
    """count dup loss"""

    if recon is None:
        recon = reconcile(tree, stree, gene2species)
    events = label_events(tree, recon)
    losses = find_loss(tree, stree, recon)
    
    dup = 0
    loss = 0
    appear = 0
    
    # count appearance    
    recon[tree.root].data["appear"] += 1
    appear += 1
    
    # count dups
    for node, event in events.iteritems():
        if event == "dup":
            recon[node].data['dup'] += 1
            dup += 1
        elif event == "gene":
            recon[node].data['genes'] += 1

    # count losses
    for gnode, snode in losses:
        snode.data['loss'] += 1
        loss += 1
    
    return dup, loss, appear


def count_ancestral_genes(stree):
    """count ancestral genes"""
    def walk(node):
        if not node.is_leaf():
            counts = []
            for child in node.children:
                walk(child)
                counts.append(child.data['genes'] 
                              - child.data['appear']
                              - child.data['dup'] 
                              + child.data['loss'])
            assert util.equal(* counts), str(counts)
            node.data['genes'] = counts[0]
    walk(stree.root)


def count_dup_loss_trees(trees, stree, gene2species):
    """
    Returns new species tree with dup,loss,appear,genes counts in node's data
    """
    
    stree = stree.copy()
    init_dup_loss_tree(stree)

    for tree in trees:
        count_dup_loss_tree(tree, stree, gene2species)
    count_ancestral_genes(stree)

    return stree


def write_event_tree(stree, out=sys.stdout):
    labels = {}
    for name, node in stree.nodes.iteritems():
        labels[name] = "[%s]\nD=%d,L=%d;\nG=%d;" % \
                       (str(name),
                        node.data['dup'], node.data['loss'],
                        node.data['genes'])

    treelib.draw_tree(stree, labels=labels, minlen=15, spacing=4,
                      labelOffset=-3,
                      out=out)


def dup_consistency(tree, recon, events):
    """
    Calculate duplication consistency scores for a reconcilied tree

    See Vilella2009 (Ensembl)
    """
    
    spset = {}
    def walk(node):
        for child in node.children:
            walk(child)
        if node.is_leaf():
            spset[node] = set([recon[node]])
        else:
            spset[node] = (spset[node.children[0]] |
                           spset[node.children[1]])
    walk(tree.root)
    
    conf = {}
    for node in tree:
        if events[node] == "dup":
            conf[node] = (len(spset[node.children[0]] &
                              spset[node.children[1]]) /
                          float(len(spset[node])))

    return conf


#=============================================================================
# tree rooting


def recon_root(gtree, stree, gene2species = gene2species, 
               rootby = "duploss", newCopy = True,
	       keepName = False, returnCost = False,
	       dupcost = 1, losscost = 1):
    """
    Reroot a tree by minimizing the number of duplications/losses/both
    Note that rootby trumps dupcost/losscost
    """
    
    # assert valid inputs
    assert rootby in ["dup", "loss", "duploss"], "unknown rootby value '%s'" % rootby
    assert dupcost >= 0 and losscost >= 0

    # make a consistent unrooted copy of gene tree
    if newCopy:
        gtree = gtree.copy()
        
    if len(gtree.leaves()) == 2:
        return

    if keepName:
        oldroot = gtree.root.name
    treelib.unroot(gtree, newCopy=False)
    treelib.reroot(gtree, 
                   gtree.nodes[sorted(gtree.leaf_names())[0]].parent.name, 
                   onBranch=False, newCopy=False)
   
    
    # make recon root consistent for rerooting tree of the same names
    # TODO: there is the possibility of ties, they are currently broken
    # arbitrarily.  In order to make comparison of reconRooted trees with 
    # same gene names accurate, hashOrdering must be done, for now.
    hash_order_tree(gtree, gene2species)
    
    # get list of edges to root on
    edges = []
    def walk(node):
        edges.append((node, node.parent))
        if not node.is_leaf():
            node.recurse(walk)
            edges.append((node, node.parent))
    for child in gtree.root.children:
        walk(child)

    
    # try initial root and recon    
    treelib.reroot(gtree, edges[0][0].name, newCopy=False)
    if keepName:
        gtree.rename(gtree.root.name, oldroot)
    recon = reconcile(gtree, stree, gene2species)
    events = label_events(gtree, recon)
    
    # find reconciliation that minimizes dup/loss
    minroot = edges[0]
    rootedge = sorted(edges[0])
    cost = 0
    if rootby in ["dup", "duploss"] and dupcost != 0:
        cost += count_dup(gtree, events) * dupcost
    if rootby in ["loss", "duploss"] and losscost != 0:
        cost += count_loss(gtree, stree,  recon) * losscost
    mincost = cost
    
    
    # try rooting on everything
    for edge in edges[1:]:
        if sorted(edge) == rootedge:
            continue
        rootedge = sorted(edge)
        
        node1, node2 = edge
        if node1.parent != node2:
            node1, node2 = node2, node1
        assert node1.parent == node2, "%s %s" % (node1.name, node2.name)
        
        # uncount cost
        if rootby in ["dup", "duploss"] and dupcost != 0:
            if events[gtree.root] == "dup":
                cost -= dupcost
            if events[node2] == "dup":
                cost -= dupcost
        if rootby in ["loss", "duploss"] and losscost != 0:
            cost -= len(find_loss_under_node(gtree.root, recon)) * losscost
            cost -= len(find_loss_under_node(node2, recon)) * losscost
        
        # new root and recon
        treelib.reroot(gtree, node1.name, newCopy=False, keepName=keepName)
        
        recon[node2] = reconcile_node(node2, stree, recon)
        recon[gtree.root] = reconcile_node(gtree.root, stree, recon)
        events[node2] = label_events_node(node2, recon)
        events[gtree.root] = label_events_node(gtree.root, recon)
        
        if rootby in ["dup", "duploss"] and dupcost != 0:
            if events[node2] ==  "dup":
                cost += dupcost
            if events[gtree.root] ==  "dup":
                cost += dupcost
        if rootby in ["loss", "duploss"] and losscost != 0:
            cost += len(find_loss_under_node(gtree.root, recon)) * losscost
            cost += len(find_loss_under_node(node2, recon)) * losscost
        
        # keep track of min cost
        if cost < mincost:
            mincost = cost
            minroot = edge

    
    # root tree by minroot
    if edge != minroot:
        node1, node2 = minroot
        if node1.parent != node2:
            node1, node2 = node2, node1
        assert node1.parent == node2
        treelib.reroot(gtree, node1.name, newCopy=False, keepName=keepName)
    
    if returnCost:
        return gtree, mincost
    else:
        return gtree


def midroot_recon(tree, stree, recon, events, params, generate):

    node1, node2 = tree.root.children

    specs1 = []
    specs2 = []
    
    # find nearest specs/genes
    def walk(node, specs):
        if events[node] == "dup":
            for child in node.children:
                walk(child, specs)
        else:
            specs.append(node)
    #walk(node1, specs1)
    #walk(node2, specs2)
    specs1 = node1.leaves()
    specs2 = node2.leaves()
    
    def getDists(start, end):
        exp_dist = 0
        obs_dist = 0

        sstart = recon[start]
        send = recon[end]
        while sstart != send:
            exp_dist += params[sstart.name][0]
            sstart = sstart.parent

        while start != end:
            obs_dist += start.dist
            start = start.parent

        return exp_dist, obs_dist / generate
    
    diffs1 = []
    for spec in specs1:
        if events[tree.root] == "spec":
            exp_dist1, obs_dist1 = getDists(spec, tree.root)
        else:
            exp_dist1, obs_dist1 = getDists(spec, node1)
        diffs1.append(obs_dist1 - exp_dist1)        

    diffs2 = []
    for spec in specs2:
        if events[tree.root] == "spec":
            exp_dist2, obs_dist2 = getDists(spec, tree.root)
        else:
            exp_dist2, obs_dist2 = getDists(spec, node2)
        diffs2.append(obs_dist2 - exp_dist2)
    
    totdist = (node1.dist + node2.dist) / generate

    left = node1.dist - stats.mean(diffs1)
    right =  totdist - node2.dist + stats.mean(diffs2)
    
    #print diffs1, diffs2    
    #print stats.mean(diffs1), stats.mean(diffs2)
    
    mid = util.clamp((left + right) / 2.0, 0, totdist)
    
    node1.dist = mid * generate
    node2.dist = (totdist - mid) * generate



def stree2gtree(stree, genes, gene2species):
    """Create a gene tree with the same topology as the species tree"""
    
    tree = stree.copy()
    
    for gene in genes:
        tree.rename(gene2species(gene), gene)
    return tree


#=============================================================================
# relationships
# encoded using gene name tuples


def get_gene_dups(tree, events):
    """Returns duplications as gene name tuples"""
    return set(tuple(sorted([tuple(sorted(child.leaf_names()))
                                   for child in node.children]))
               for node, kind in events.iteritems()
               if kind == "dup")

def get_speciations(tree, events):
    """Returns speciations as gene name tuples"""
    return set(tuple(sorted([tuple(sorted(child.leaf_names()))
                                   for child in node.children]))
               for node, kind in events.iteritems()
               if kind == "spec")


def get_gene_losses(tree, stree, recon):
    """Returns losses as gene name, species name tuples"""
    return set((loss[0].name, loss[1].name)
               for loss in find_loss(tree, stree, recon))
         

def get_orthologs(tree, events):
    """Returns orthologs as gene name pairs"""
    
    specs = [sorted([sorted(child.leaf_names())
                     for child in node.children])
             for node in events
             if events[node] == "spec"]
    
    return set(tuple(sorted((a, b)))
               for x in specs
               for a in x[0]
               for b in x[1])


#=============================================================================
# Tree hashing
#

def hash_tree_compose(child_hashes, node=None):
    return "(%s)" % ",".join(child_hashes)

def hash_tree_compose_names(child_hashes, node=None):
     return "(%s)%s" % (",".join(child_hashes), node.name)

def hash_tree(tree, smap=lambda x: x, compose=hash_tree_compose):
    def walk(node):
        if node.is_leaf():
            return smap(node.name)
        else:
            child_hashes = map(walk, node.children)
            child_hashes.sort()
            return compose(child_hashes, node)
    
    if isinstance(tree, treelib.Tree) or hasattr(tree, "root"):
        return walk(tree.root)
    elif isinstance(tree, treelib.TreeNode):
        return walk(tree)
    else:
        raise Exception("Expected Tree object")

def hash_tree_names(tree, smap=lambda x: x, compose=hash_tree_compose_names):
    return hash_tree(tree, smap, compose)

def hash_order_tree(tree, smap = lambda x: x):
    def walk(node):
        if node.is_leaf():
            return smap(node.name)
        else:
            child_hashes = map(walk, node.children)
            ind = util.sortindex(child_hashes)
            child_hashes = util.mget(child_hashes, ind)
            node.children = util.mget(node.children, ind)
            return hash_tree_compose(child_hashes)
    walk(tree.root)



#=============================================================================
# Branch length distributions
#


def mapRefTree(trees, reftree, refmapfunc):
    collect = util.Dict(1, [])
    
    for tree in trees:
        nodemap = refmapfunc(tree, reftree)
        
        for name, node in tree.nodes.iteritems():
            collect[nodemap[name]].append(node)
    
    return collect


def findBranchLengths(collect):
    return util.mapdict(collect, valfunc = lambda nodes: 
                        map(lambda node: node.dist, nodes))

def findTreeLengths(collect):
    totals = map(sum, util.map2(lambda x: x.dist, zip(* collect.values())))
    return totals
    



#=============================================================================
# Branch length analysis
#


def get_species_inorder(tree):
    return [node.name for node in tree.inorder()]


def get_branch_lens(trees, stree, gene2species=gene2species):
    # determine species nanes
    species = map(str, stree.nodes.keys())
    species.remove(str(stree.root.name))
    
    # make rates table
    rates = tablelib.Table(headers=species)
    
    # loop through trees
    for tree in trees:
        if isinstance(tree, str):
            tree = treelib.read_tree(tree)
        recon = reconcile(tree, stree, gene2species)
        events = label_events(tree, recon)
        
        # skip trees with duplications or with extremly long branch lengths
        assert "dup" not in events.values()
        
        row = {}
        for node in tree.nodes.values():
            row[str(recon[node].name)] = node.dist
        rates.append(row)
    
    return rates


def find_branch_distrib(trees, stree, gene2species = gene2species, 
                      relative = True):
    """Older version of getBranchLens()
    
       Will probably be deprecated soon.
    """
    
    lengths = util.Dict(1, [])
    used = []

    for tree in trees:
        recon = reconcile(tree, stree, gene2species)
        events = label_events(tree, recon)
        
        # skip trees with duplications or with extremly long branch lengths
        if "dup" in events.values():
            used.append(False)
            continue
        else:
            used.append(True)
        
        for node in tree.nodes.values():
            if relative:
                # find total length of tree
                totalLength = 0
                for node in tree.nodes.values():
                    totalLength += node.dist
            
                lengths[recon[node]].append(node.dist/totalLength)
            else:                
                lengths[recon[node]].append(node.dist)
    
    
    return lengths, used



def getRelBranchLens(rates, species=None):
    if species == None:
        species = rates.headers
    
    nonspecies = set(rates.headers) - set(species)
    
    relrates = rates.new()
    
    for row in rates:
        row2 = {}
        tot = sum(util.mget(row, species))
        
        for sp in species:
            row2[sp] = row[sp] / tot
        
        # copy over non-species data
        for key in nonspecies:
            row2[key] = row[key]
        
        relrates.append(row2)
    
    return relrates
        

def getBranchZScores(rates, params):
    zscores = rates.new()
    
    # determine column to species mapping
    col2species = {}
    for species in params:
        col2species[str(species)] = species
    
    for row in rates:
        row2 = {}
        for key, val in row.iteritems():
            if key in col2species:
                # compute zscore
                mu, sigma = params[col2species[key]]
                row2[key] = (val - mu) / sigma
            else:
                # not a branch, copy value unchanged
                row2[key] = val
        
        zscores.append(row2)
    
    return zscores


#=============================================================================
# add implied speciation nodes to a gene tree



def add_spec_node(node, snode, tree, recon, events):
    """
    insert new speciation node above gene node 'node' from gene tree 'tree'

    new node reconciles to species node 'snode'.  Modifies recon and events
    accordingly
    """
    
    newnode = treelib.TreeNode(tree.new_name())
    parent = node.parent
    
    # find index of node in parent's children
    nodei = parent.children.index(node)
    
    # insert new node into tree
    tree.add_child(parent, newnode)
    parent.children[nodei] = newnode
    parent.children.pop()
    tree.add_child(newnode, node)
    
    # add recon and events info
    recon[newnode] = snode
    events[newnode] = "spec"

    return newnode

def remove_spec_node(node, tree):
    """
    removes speciation node 'node' from gene tree 'tree'

    Modifies recon and events accordingly
    """
    assert len(node.children) == 1
    
    # remove node from tree
    tree.add_child(node.parent, node.children[0])
    tree.remove(node)
    
    # remove recon and events info
    del recon[node]
    del events[node]


def add_implied_spec_nodes(tree, stree, recon, events):
    """
    adds speciation nodes to tree that are implied but are not present
    because of gene losses
    """
    
    added_nodes = []

    for node in list(tree):
        # process this node and the branch above it

        # handle root node specially
        if node.parent is None:
            # ensure root of gene tree properly reconciles to
            # root of species tree
            if recon[node] == stree.root:                            
                continue            
            tree.root = treelib.TreeNode(tree.new_name())
            tree.add_child(tree.root, node)
            recon[tree.root] = stree.root
            events[tree.root] = "spec"
            added_nodes.append(tree.root)

        # determine starting and ending species
        sstart = recon[node]
        send = recon[node.parent]

        # the species path is too short to have implied speciations
        if sstart == send:
            continue

        parent = node.parent

        # determine species path of this gene branch (node, node->parent)
        snode = sstart.parent

        while snode != send:
            added_nodes.append(add_spec_node(node, snode, tree, recon, events))
            node = node.parent
            snode = snode.parent
        
    
        # determine whether node.parent is a dup
        # if so, send (a.k.a. species end) is part of species path
        if events[parent] == "dup":
            added_nodes.append(add_spec_node(node, send, tree, recon, events))

    return added_nodes


def remove_implied_spec_nodes(tree, added_nodes):
    """
    removes speciation nodes from tree
    """
    for node in added_nodes:
        remove_spec_node(tree, node)


#=============================================================================
# local rearrangements


def perform_nni(tree, node1, node2, change=0, rooted=True):
    """Proposes a new tree using Nearest Neighbor Interchange
       
       Branch for NNI is specified by giving its two incident nodes (node1 and 
       node2).  Change specifies which  subtree of node1 will be swapped with
       the uncle.  See figure below.

         node2
        /     \
      uncle    node1
               /  \
         child[0]  child[1]
    
    special case with rooted branch and rooted=False:
    
              node2
             /     \
        node2'      node1
       /     \     /     \
      uncle   * child[0] child[1]
    
    """
    
    if node1.parent != node2:
        node1, node2 = node2, node1  
    
    # try to see if edge is one branch (not root edge)
    if not rooted and treelib.is_rooted(tree) and \
       node2 == tree.root:
        # special case of specifying root edge
        if node2.children[0] == node1:
            node2 = node2.children[1]
        else:
            node2 = node2.children[0]
        
        # edge is not an internal edge, give up
        if len(node2.children) < 2:
            return
        
    if node1.parent == node2.parent == tree.root:
        uncle = 0
        
        if len(node2.children[0].children) < 2 and \
           len(node2.children[1].children) < 2:
            # can't do NNI on this branch
            return
    else:   
        assert node1.parent == node2
    
        # find uncle
        uncle = 0 
        if node2.children[uncle] == node1:
            uncle = 1
    
    # swap parent pointers
    node1.children[change].parent = node2
    node2.children[uncle].parent = node1
    
    # swap child pointers
    node2.children[uncle], node1.children[change] = \
        node1.children[change], node2.children[uncle]



def propose_random_nni(tree):
    """
    Propose a random NNI rearrangement
    """

    nodes = tree.nodes.values()

    # find edges for NNI
    while True:
        node1 = random.sample(nodes, 1)[0]
        if not node1.is_leaf() and node1.parent is not None:
            break
    
    node2 = node1.parent
    #a = node1.children[random.randint(0, 1)]
    #b = node2.children[1] if node2.children[0] == node1 else node2.children[0]
    #assert a.parent.parent == b.parent

    return node1, node2, random.randint(0, 1)




def perform_spr(tree, subtree, newpos):
    """
    Proposes new topology using Subtree Pruning and Regrafting (SPR)
    
        a = subtree
        e = newpos

        BEFORE
                ....
            f         d
           /           \
          c             e
         / \           ...
        a   b
       ... ...

        AFTER

            f         d
           /           \
          b             c
         ...           / \
                      a   e
                     ... ...

        Requirements:
        1. a (subtree) is not root or children of root
        2. e (newpos) is not root, a, descendant of a, c (parent of a), or 
           b (sibling of a)
        3. tree is binary

"""
    a = subtree
    e = newpos

    c = a.parent
    f = c.parent
    bi = 1 if c.children[0] == a else 0
    b = c.children[bi]
    ci = 0  if f.children[0] == c else 1
    d = e.parent
    ei = 0 if d.children[0] == e else 1

    d.children[ei] = c
    c.children[bi] = e
    f.children[ci] = b
    b.parent = f
    c.parent = d
    e.parent = c




def propose_random_spr(tree):
    """
    What if e == f  (also equivalent to NNI) this is OK

    BEFORE
    
          d
         / \
        e  ...
       / \
      c  ...         
     / \           
    a   b
   ... ...

    AFTER
          d
         / \
        c
       / \
      a   e
     ... / \
        b  ...
       ...
       
  What if d == f  (also equivalent to NNI) this is OK
  
    BEFORE
          
        f
       / \
      c   e
     / \  ...
    a   b
   ... ...

    AFTER
          
        f
       / \
      b   c  
     ... / \ 
        a   e
       ... ...  

    Requirements:
    1. a (subtree) is not root or children of root
    2. e (newpos) is not root, a, descendant of a, c (parent of a), or 
       b (sibling of a)
    3. tree is binary
    """

    assert len(tree.nodes) >= 5, "Tree is too small"

    # find subtree (a) to cut off (any node that is not root or child of root)
    nodes = tree.nodes.values()
    while True:
        a = random.sample(nodes, 1)[0]
        if (a.parent is not None and a.parent.parent is not None):
            break
    subtree = a
    
    # find sibling (b) of a
    c = a.parent
    bi = 1 if c.children[0] == a else 0
    b = c.children[bi]
    
    # choose newpos (e)
    e = None
    while True:
        e = random.sample(nodes, 1)[0]
        
        # test if e is a valid choice
        if e.parent is None or e == a or e == c or e == b:
            continue
        
        # also test if e is a descendent of a
        under_a = False
        ptr = e.parent
        while ptr is not None:
            if ptr == a:
                under_a = True
                break
            ptr = ptr.parent
        if under_a:
            continue
        
        break
    newpos = e

    return subtree, newpos




#=============================================================================
# reconciliation rearrangements

def change_recon_up(recon, node, events=None):
    """
    Move the mapping of a node up one branch
    """
        
    if events is not None and events[node] == "spec":
        # promote speciation to duplication
        # R'(v) = e(R(u))
        events[node] = "dup"
    else:
        # R'(v) = p(R(u))
        recon[node] = recon[node].parent


def change_recon_down(recon, node, schild, events=None):
    """
    Move the mapping of a node down one branch
    """

    if events is not None and recon[node] == schild:
        events[node] = "spec"
    else:
        recon[node] = schild


def can_change_recon_up(recon, node, events=None):
    """Returns True is recon can remap node one 'step' up"""

    if events is not None and events[node] == "spec" and not node.is_leaf():
        # promote speciation to duplication
        return True
    else:
        # move duplication up one branch
        rnode = recon[node]
        prnode = rnode.parent

        # rearrangement is valid if
        return (not node.is_leaf() and 
            prnode is not None and #  1. there is parent sp. branch
            (node.parent is None or # 2. no parent to restrict move
             rnode != recon[node.parent] # 3. not already matching parent
             ))


def enum_recon(tree, stree, depth=None,
               step=0, preorder=None,
               recon=None, events=None,
               gene2species=None):
    """
    Enumerate reconciliations between a gene tree and a species tree
    """
    
    if recon is None:
        recon = reconcile(tree, stree, gene2species)
        events = label_events(tree, recon)

    if preorder is None:
        preorder = list(tree.preorder())

    # yield current recon
    yield recon, events

    if depth is None or depth > 0:
        for i in xrange(step, len(preorder)):
            node = preorder[i]
            if can_change_recon_up(recon, node, events):
                schild = recon[node]
                change_recon_up(recon, node, events)
            
                # recurse
                depth2 = depth - 1 if depth is not None else None
                for r, e in enum_recon(tree, stree, depth2,
                                       i, preorder,
                                       recon, events):
                    yield r, e
            
                change_recon_down(recon, node, schild, events)



    
'''
class EnumRecon (object):
    """
    Enumerate reconciliations between a gene tree and species tree
    """

    def __init__(self, tree, stree, depth=1,
                 step=-1, preorder=None,
                 recon=None, events=None,
                 gene2species=None):
        self.tree = tree
        self.stree = stree
        self.depth = 1
        self.step = step

        if recon:
            self.recon = recon
            self.events = events
        else:
            self.recon = reconcile(tree, stree, gene2species)
            self.events = label_events(tree, self.recon)

        if preorder:
            self.preorder = preorder
        else:
            self.preorder = list(tree.preorder())
            
        
        self.sprev = None


    def __iter__(self):
        return self

    
    def next(self):
        if self.step >= 0:
            # perform a rearrangement
            if self.step >= len(self.preorder):
                # no more mappings to move up, take 1 step back
                change_recon_down(self.recon, node, schild, self.events)
            
            node = self.preorder[self.step]
            if can_change_recon_up(self.recon, node, self.events):
                self.sprev = self.recon[node]
                change_recon_up(self.recon, node, self.events)
        self.step += 1
        
        return self.recon
'''



#=============================================================================
# tree search

class TreeSearch (object):

    def __init__(self, tree):
        self.tree = tree

    def __iter__(self):
        return self

    def set_tree(self, tree):
        self.tree = tree

    def get_tree(self):
        return self.tree

    def propose(self):
        raise

    def revert(self):
        raise

    def reset(self):
        pass

    def next(self):
        return self.propose()


class TreeSearchNni (TreeSearch):

    def __init__(self, tree):
        TreeSearch.__init__(self, tree)
        self.set_tree(tree)

    def set_tree(self, tree):
        self.tree = tree
        self.node1 = None
        self.node2 = None
        self.child = None

    def propose(self):
        self.node1, self.node2, self.child = propose_random_nni(self.tree)
        perform_nni(self.tree, self.node1, self.node2, self.child)
        return self.tree

    def revert(self):
        if self.node1 is not None:
            perform_nni(self.tree, self.node1, self.node2, self.child)
        return self.tree

    def reset(self):
        self.node1 = None
        self.node2 = None
        self.child = None


class TreeSearchSpr (TreeSearch):

    def __init__(self, tree):
        TreeSearch.__init__(self, tree)
        self.set_tree(tree)

    def set_tree(self, tree):
        self.tree = tree
        self.node1 = None
        self.node2 = None

    def propose(self):

        # choose SPR move
        self.node1, node3 = propose_random_spr(self.tree)
        
        # remember sibling of node1
        p = self.node1.parent
        self.node2 = (p.children[1] if p.children[0] == self.node1
                      else p.children[0])

        # perform SPR move
        perform_spr(self.tree, self.node1, node3)
        return self.tree

    def revert(self):
        if self.node1 is not None:
            perform_spr(self.tree, self.node1, self.node2)
        return self.tree

    def reset(self):
        self.node1 = None
        self.node2 = None        


class TreeSearchMix (TreeSearch):

    def __init__(self, tree):
        TreeSearch.__init__(self, tree)
        self.total_weight = 0.0
        self.last_propose = 0
        self.methods = []
        self.set_tree(tree)

    def set_tree(self, tree):
        self.tree = tree
	for method in self.methods:
            method[0].set_tree(tree)

    def add_proposer(self, proposer, weight):
        self.total_weight += weight
        self.methods.append((proposer, weight))

    def propose(self):
        # randomly choose method
        choice = random.random() * self.total_weight
        s = self.methods[0][1]
        i = 0
        while i < len(self.methods)-1 and s < choice:
            i += 1
            s += self.methods[i][1]

        # make proposal
        self.last_propose = i
        self.tree = self.methods[i][0].propose()
	return self.tree

    def revert(self):
        self.tree = self.methods[self.last_propose][0].revert()
	return self.tree

    def reset(self):
        for method in self.methods:
            method[0].reset()


class TreeSearchUnique (TreeSearch):
    """
    Propose unique tree topologies
    """

    def __init__(self, tree, search, tree_hash=None, maxtries=5,
                 auto_add=True):
        TreeSearch.__init__(self, tree)
        self.search = search
        self.seen = set()
        self._tree_hash = tree_hash if tree_hash else hash_tree
        self.maxtries = maxtries
        self.auto_add = auto_add
        self.set_tree(tree)

    def set_tree(self, tree):
        self.tree = tree
        self.search.set_tree(tree)


    def propose(self):

        for i in xrange(self.maxtries):
            if i > 0:
                self.search.revert()
            tree = self.search.propose()
            top = self._tree_hash(tree)
            if top not in self.seen:
                #util.logger("tried", i, len(self.seen))
                break
        else:            
            #util.logger("maxtries", len(self.seen))
	    pass

        if self.auto_add:
            self.seen.add(top)
        self.tree = tree
        return self.tree
        

    def revert(self):
        self.tree = self.search.revert()
        return self.tree


    def reset(self):
        self.seen.clear()
        self.search.reset()


    def add_seen(self, tree):
        top = self._tree_hash(tree)
        self.seen.add(top)
        

class TreeSearchPrescreen (TreeSearch):

    def __init__(self, tree, search, prescreen, poolsize):
        TreeSearch.__init__(self, tree)
        self.search = TreeSearchUnique(tree, search, auto_add=False)
        self.prescreen = prescreen
        self.poolsize = poolsize
        self.oldtree = None
        self.set_tree(tree)


    def set_tree(self, tree):
        self.tree = tree
        self.search.set_tree(tree)


    def propose(self):

        # save old topology
        self.oldtree = self.tree.copy()

        pool = []
        best_score = self.prescreen(self.tree)
        total = -util.INF

        # TODO: add unique filter

        # make many subproposals
        self.search.reset()
        for i in xrange(self.poolsize):
            self.search.propose()
            score = self.prescreen(self.tree)
            tree = self.tree.copy()

            # save tree and logl
            pool.append((tree, score))
            total = stats.logadd(total, score)

            if score > best_score:
                # make more proposals off this one
                best_score = score
            else:
                self.search.revert()

        # propose one of the subproposals 
        choice = random.random()
        partsum = -util.INF

        for tree, score in pool:
            partsum = stats.logadd(partsum, score)
            if choice < math.exp(partsum - total):
                # propose tree i
                treelib.set_tree_topology(self.tree, tree)
                break


        self.search.add_seen(self.tree)


    def revert(self):
        if self.oldtree:
            treelib.set_tree_topology(self.tree, self.oldtree)


    def reset(self):
        self.oldtree = None
        self.search.reset()




#=============================================================================
# Phylogenetic reconstruction: Neighbor-Joining
#

def neighborjoin(distmat, genes, usertree=None, binary=False):
    """Neighbor joining algorithm"""
    
    tree = treelib.Tree()
    leaves = {}
    dists = util.Dict(2, None)
    restdists = {}
    
    
    # initialize distances
    for i in range(len(genes)):
        r = 0
        for j in range(len(genes)):
            dists[genes[i]][genes[j]] = distmat[i][j]
            r += distmat[i][j]
        restdists[genes[i]] = r / (len(genes) - 2)
        
    # initialize leaves
    for gene in genes:
        tree.add(treelib.TreeNode(gene))
        leaves[gene] = 1
    
    # if usertree is given, determine merging order
    merges = []
    newnames = {}
    if usertree != None:
        def walk(node):
            if not node.is_leaf():
                assert len(node.children) == 2, \
                    Exception("usertree is not binary")
            
                for child in node:
                    walk(child)
                merges.append(node)
                newnames[node] = len(merges)
            else:
                newnames[node] = node.name
        walk(usertree.root)
        merges.reverse()
    
    # join loop
    while len(leaves) > 2:
        # search for closest genes
        if not usertree:
            low = util.INF
            lowpair = (None, None)
            leaveslst = leaves.keys()

            for i in range(len(leaves)):
                for j in range(i+1, len(leaves)):
                    gene1, gene2 = leaveslst[i], leaveslst[j]
                    dist = dists[gene1][gene2] - restdists[gene1] \
                                               - restdists[gene2]
                    if dist < low:
                        low = dist
                        lowpair = (gene1, gene2)
        else:
            node = merges.pop()
            lowpair = (newnames[node.children[0]],
                       newnames[node.children[1]])
        
        # join gene1 and gene2
        gene1, gene2 = lowpair
        parent = treelib.TreeNode(tree.new_name())
        tree.add_child(parent, tree.nodes[gene1])
        tree.add_child(parent, tree.nodes[gene2])
        
        # set distances
        tree.nodes[gene1].dist = (dists[gene1][gene2] + restdists[gene1] - 
                                  restdists[gene2]) / 2.0
        tree.nodes[gene2].dist = dists[gene1][gene2] - tree.nodes[gene1].dist
        
        # gene1 and gene2 are no longer leaves
        del leaves[gene1]
        del leaves[gene2]
        
        gene3 = parent.name
        r = 0
        for gene in leaves:
            dists[gene3][gene] = (dists[gene1][gene] + dists[gene2][gene] -
                                  dists[gene1][gene2]) / 2.0
            dists[gene][gene3] = dists[gene3][gene]
            r += dists[gene3][gene]
        leaves[gene3] = 1
        
        if len(leaves) > 2:
            restdists[gene3] = r / (len(leaves) - 2)

    gene1, gene2 = leaves.keys()
    if binary:
        # join the last two genes at the root and split the remaining dist evenly
        parent = treelib.TreeNode(tree.new_name())
        tree.add_child(parent, tree.nodes[gene1])
        tree.add_child(parent, tree.nodes[gene2])
        tree.nodes[gene1].dist = dists[gene1][gene2] / 2.0
        tree.nodes[gene2].dist = dists[gene1][gene2] / 2.0
        tree.root = parent
    else:
        # join the last two genes into a tribranch
        if type(gene1) != int:
            gene1, gene2 = gene2, gene1
        tree.add_child(tree.nodes[gene1], tree.nodes[gene2])
        tree.nodes[gene2].dist = dists[gene1][gene2]
        tree.root = tree.nodes[gene1]

    # root tree according to usertree    
    if usertree != None and treelib.is_rooted(usertree):
        roots = set([newnames[usertree.root.children[0]],
                     newnames[usertree.root.children[1]]])
        newroot = None
        for child in tree.root.children:
            if child.name in roots:
                newroot = child
        
        assert newroot != None
        
        treelib.reroot(tree, newroot.name, newCopy=False)
    
    return tree




#=============================================================================
# Phylogenetic reconstruct: Least Square Error

def least_square_error(tree, distmat, genes, forcePos=True, weighting=False):
    """Least Squared Error algorithm for phylogenetic reconstruction"""
    
    # use SCIPY to perform LSE
    import scipy
    import scipy.linalg
    
    def makeVector(array):
        """convience function for handling different configurations of scipy"""
        if len(array.shape) == 2:
            if array.shape[0] == 1:
                return array[0]
            else:
                return scipy.transpose(array)[0]
        else:
            return array
            
    
    if treelib.is_rooted(tree):
        rootedge = sorted([x.name for x in tree.root.children])
        treelib.unroot(tree, newCopy=False)
    else:
        rootedge = None        
    
    # create pairwise dist array
    dists = []
    for i in xrange(len(genes)):
        for j in xrange(i+1, len(genes)):
            dists.append(distmat[i][j])
    
    # create topology matrix
    topmat, edges = make_topology_matrix(tree, genes)
    
    # setup matrix and vector
    if weighting:
        topmat2 = scipy.array([[util.safediv(x, math.sqrt(dists[i]), 0) 
                                for x in row]
                               for i, row in enumerate(topmat)])
        paths = scipy.array(map(math.sqrt, dists))
    else:
        topmat2 = scipy.array(topmat)
        paths = scipy.array(dists)

    
    # solve LSE
    edgelens, resids, rank, singlars = scipy.linalg.lstsq(topmat2, paths)
    
    # force non-negative branch lengths
    if forcePos:
        edgelens = [max(float(x), 0) for x in makeVector(edgelens)]
    else:
        edgelens = [float(x) for x in makeVector(edgelens)]
    
    # calc path residuals (errors)
    paths2 = makeVector(scipy.dot(topmat2, edgelens))
    resids = (paths2 - paths).tolist()
    paths = paths.tolist()
    
    # set branch lengths
    set_branch_lengths_from_matrix(tree, edges, edgelens, paths, resids, 
                                   topmat=topmat, rootedge=rootedge)
    
    return util.Bundle(resids=resids, 
                       paths=paths, 
                       edges=edges, 
                       topmat=topmat)


def make_topology_matrix(tree, genes):

    # find how edges split vertices
    network = treelib.tree2graph(tree)
    splits = find_all_branch_splits(network, set(genes))
    edges = splits.keys()

    # create topology matrix
    n = len(genes) 
    ndists = n*(n-1) / 2
    topmat = util.makeMatrix(ndists, len(edges))
    
    vlookup = util.list2lookup(genes)
    n = len(genes)
    for e in xrange(len(edges)):
        set1, set2 = splits[edges[e]]
        for gene1 in set1:
            for gene2 in set2:
                i, j = util.sort([vlookup[gene1], vlookup[gene2]])
                index = i*n-i*(i+1)/2+j-i-1
                topmat[index][e] = 1.0
    
    return topmat, edges


def set_branch_lengths_from_matrix(tree, edges, edgelens, paths, resids, 
                                   topmat=None, rootedge=None):
    # recreate rooting branches
    if rootedge != None:
        # restore original rooting
        if tree.nodes[rootedge[0]].parent == tree.nodes[rootedge[1]]:
            treelib.reroot(tree, rootedge[0], newCopy=False)
        else:
            treelib.reroot(tree, rootedge[1], newCopy=False)
    
        # find root edge in edges
        for i in xrange(len(edges)):
            if sorted(edges[i]) == rootedge:
                break
                
        edges[i] = [rootedge[0], tree.root.name]
        edges.append([rootedge[1], tree.root.name])
        edgelens[i] /= 2.0
        edgelens.append(edgelens[i])
        resids[i] /= 2.0
        resids.append(resids[i])
        paths[i] /= 2.0
        paths.append(paths[i])
        
        if topmat != None:
            for row in topmat:
                row.append(row[i])
    
    # set branch lengths
    for i in xrange(len(edges)):
        gene1, gene2 = edges[i]
        if tree.nodes[gene2].parent == tree.nodes[gene1]:
            gene1, gene2 = gene2, gene1
        tree.nodes[gene1].dist = edgelens[i]




def tree2distmat(tree, leaves):
    """Returns pair-wise distances between leaves of a tree"""
    
    # TODO: not implemented efficiently
    mat = []
    for i in range(len(leaves)):
        mat.append([])
        for j in range(len(leaves)):
            mat[-1].append(treelib.find_dist(tree, leaves[i], leaves[j]))
    
    return mat



#============================================================================
# branch splits
#

def find_all_branch_splits(network, leaves):
    # find vertice and edge visit history
    start = network.keys()[0]

    openset = [start]
    closedset = {}
    
    vhistory = []
    ehistory = []
    elookup = util.Dict(1, [])
    
    
    while len(openset) > 0:
        vertex = openset.pop()
        
        vhistory.append(vertex)
        
        if len(vhistory) > 1:
            edge = tuple(util.sort(vhistory[-2:]))        
            ehistory.append(edge)
            elookup[edge].append(len(ehistory) - 1)
        
        # skip closed vertices
        if vertex in closedset:
            continue
        
        for v in network[vertex].keys():
            if v not in closedset:
                openset.append(vertex)            
                openset.append(v)
        

        # close new vertex
        closedset[vertex] = 1
    
    
    # use histories to define half each split
    splits = {}
    for edge in elookup:
        set1 = {}
        
        start, end = elookup[edge]
        for i in range(start+1, end+1):
            if vhistory[i] in leaves:
                set1[vhistory[i]] = 1
        
        # fill in other half of splits using complement
        set2 = {}
        for v in leaves:
            if v not in set1:
                set2[v] = 1
        
        if edge[0] == vhistory[start]:
            splits[edge] = [set2, set1]
        else:
            splits[edge] = [set1, set2]
        
    
    return splits


def find_branch_splits(tree):
    splits = find_all_branch_splits(treelib.tree2graph(tree),
                                    tree.leaf_names())
    splits2 = {}
    
    for edge, sets in splits.iteritems():
        # skip external edges
        if len(sets[0]) == 1 or len(sets[1]) == 1:
            continue
        
        s = tuple(sorted([tuple(sorted(i.keys())) for i in sets]))
        splits2[edge] = s
    
    # if tree is rooted, remove duplicate edge
    if treelib.is_rooted(tree):
        edge1 = tuple(sorted([tree.root.name, tree.root.children[0].name]))
        edge2 = tuple(sorted([tree.root.name, tree.root.children[1].name]))
        if edge1 > edge2:
            edge1, edge2 = edge2, edge1
        if edge1 in splits2 and edge2 in splits2:
            del splits2[edge1]
    
    return splits2


def find_splits(tree):
    """Find branch splits for a tree"""
    
    allLeaves = set(tree.leaf_names())

    # find descendants
    descendants = []
    def walk(node):
        if node.is_leaf():
            descendants.append(set([node.name]))
        else:
            s = set()
            for child in node.children:
                s.update(walk(child))
            descendants.append(s)
        return descendants[-1]
    for child in tree.root.children:
        walk(child)

    # left child's descendants immediately defines
    # right child's descendants (by complement)
    if len(tree.root.children) == 2:
        descendants.pop()

    # build splits list
    splits = []
    for leaves in descendants:
        if len(leaves) > 1:
            set1 = tuple(sorted(leaves))
            set2 = tuple(sorted(allLeaves - leaves))
            if len(set1) > len(set2):
                set1, set2 = set2, set1
                
            splits.append((set1, set2))
    
    return splits


def split_string(split, leaves=None, leafDelim=" ", splitDelim="|"):
    """
    Returns a string representing a split

    If leaves are specified, leaf names will be displayed in that order.
    """

    if leaves is not None:
        lookup = util.list2lookup(leaves)
        split = (sorted(split[0], key=lambda x: lookup[x]),
                 sorted(split[0], key=lambda x: lookup[x]))

    return leafDelim.join(split[0]) + splitDelim + leafDelim.join(split[1])


def split_bit_string(split, leaves=None, char1="*", char2=".", nochar=" "):
    """Returns a bit string representation of a split"""

    if leaves is None:
        leaves = split[0] + split[1]
    set1, set2 = map(set, split)

    chars = []
    for leaf in leaves:
        if leaf in set1:
            chars.append(char1)
        elif leaf in set2:
            chars.append(char2)
        else:
            chars.append(nochar)

    return "".join(chars)
    
    


def robinson_foulds_error(tree1, tree2):
    """Returns the normalized Robinson Foulds Error metric, e.g. (A+B)/max(C,D),
       where A is the number of partitions implied by tree1 but not by tree2,
       B is the number of partitions implied by tree2 but not by tree1,
       C is the total number of partitions implied by tree1, and
       D is the total number of partitions implied by tree2.
    """
    splits1 = find_branch_splits(tree1)
    splits2 = find_branch_splits(tree2)

    overlap = set(splits1.values()) & set(splits2.values())
    
    #assert len(splits1) == len(splits2)

    denom = float(max(len(splits1), len(splits2)))
    
    if denom == 0.0:
        return 0.0
    else:
        return 1 - (len(overlap) / denom)


#=============================================================================
# file functions

def phylofile(famdir, famid, ext):
    """Creates a filename using my gene family format

    famdir/famid/famid.ext
    """
    return os.path.join(famdir, famid, famid + ext)



#=============================================================================
# visualization

def view_tree(tree, options = "-t 1"):
    tmpfile = util.tempfile(".", "vistree", ".tree")
    tree.write(tmpfile)
    os.system("vistree.py -n %s %s" % (tmpfile, options))
    os.remove(tmpfile)
viewTree = view_tree



#=============================================================================
# miscellaneous code that is not used frequently



def findRootedSubtrees(tree):
    """tree is unrooted"""
    
    trees = []
    
    # convert tree to graph
    mat = treelib.tree2graph(tree)
    
    donei = {}
    
    # loop over branches, and collect rooted trees on each node of branch
    for i in mat:
        donei[i] = 1
        for j in mat[i]:
            if j not in donei:
                tree1 = treelib.graph2tree(mat, i, closedset={j:1})
                tree2 = treelib.graph2tree(mat, j, closedset={i:1})
                treelib.removeSingleChildren(tree1)
                treelib.removeSingleChildren(tree2)
                trees.append(tree1)
                trees.append(tree2)
    
    return trees


def findOrthoNeighbors(parts, hits):

    from . import blast

    lookup = {}
    for i in xrange(len(parts)):
        for item in parts[i]:
            lookup[item] = i

    mat = util.Dict(1, (0, None))
    
    # find unidirectional best hits at partition level
    for hit in hits:
        try:
            part1 = lookup[blast.query(hit)]
            part2 = lookup[blast.subject(hit)]
            score = blast.bitscore(hit)
        except:
            # skip genes not in parts
            continue
        
        # don't count hits within a cluster
        if part1 == part2:
            continue
        
        if score > mat[part1][0]:
            mat[part1] = (score, part2)
        
        if score > mat[part2][0]:
            mat[part2] = (score, part1)
    
    # find best bidirectional hits
    nbrs = []
    touched = {}
    for part in xrange(len(parts)):
        other = mat[part][1]
        if mat[other][1] == part and \
           part not in touched and \
           other not in touched:
            nbrs.append((part, other))
            touched[part] = 1
    
    return nbrs



def partition_tree(tree, stree, gene2species):
    recon = reconcile(tree, stree, gene2species)
    sroots = find_species_roots(tree, stree, recon)
    
    # extend subroots
    sroots = treelib.maxDisjointSubtrees(tree, sroots)
    
    # partition
    trees = []
    for sroot in sroots:
        trees.append(treelib.subtree(tree, sroot))
        trees[-1].root.parent = None
    
    return trees
partitionTree = partition_tree


def findSpeciesSets(stree):
    """
    returns a mapping for each species tree node 
    to a set of modern day species
    """
    
    setmap = {}
    def walk(node):
        node.recurse(walk)
        if node.is_leaf():
            setmap[node] = {node.name:1}
        else:
            setmap[node] = {}
            for child in node.children:
                setmap[node].update(setmap[child])
    walk(stree.root)
    
    return setmap



def jukesCantorCorrection(dist):
    """Applies the Jukes Cantor correction to distances
       
       Only valid for distances less than .75 sub/site
    """
    assert (dist < .75)
    return - (3/4.0) * log(1 - (4/3.) * dist)


'''
#=============================================================================
# Sequence Distance Estimation


def getSeqPairDist(seq1, seq2, infile=None, outfile=None):
    aln = fasta.FastaDict()
    aln["0"] = seq1
    aln["1"] = seq2
    
    if os.path.isfile("infile"):
        raise Exception("infile already exists")
    
    # force PHYLIP to ask for outfile
    if not os.path.exists("outfile"):
        file("outfile", "w").close()
        madePhylip = True
    else:
        madePhylip = False
    
    
    
    # write file
    if infile == None:
        infile = util.tempfile(".", "tmp_in", ".align")
        madeInfile = True
    else:
        madeInfile = False
    if outfile == None:    
        outfile = util.tempfile(".", "tmp_out", ".dist")
        madeOutfile = True
    else:
        madeOutfile = False
    
    if os.path.exists(outfile):
        args = "%s\nf\n%s\nr\ny\n" % (infile, outfile)
    else:
        args = "%s\nf\n%s\ny\n" % (infile, outfile)
    
    phylip.writePhylipAlign(file(infile, "w"), aln)
    phylip.execPhylip("dnadist", args, verbose=False)
    labels, distmat = phylip.read_dist_matrix(outfile)

    if madePhylip:
        os.remove("outfile")
    
    if madeInfile:
        os.remove(infile)
    if madeOutfile:
        os.remove(outfile)
    
    return distmat[0][1]


def getGaplessDistMatrix(aln):
    infile = util.tempfile("/tmp/", "tmp_in", ".align")
    outfile = util.tempfile("/tmp/", "tmp_out", ".dist")
    
    # force PHYLIP to ask for outfile
    if not os.path.exists("outfile"):
        file("outfile", "w").close()
        madeOutfile = True
    else:
        madeOutfile = False
    
    distmat = util.makeMatrix(len(aln), len(aln), 0.0)
    keys = aln.keys()
    
    for i in xrange(0, len(aln)):
        for j in xrange(i+1, len(aln)):
            distmat[i][j] = getSeqPairDist(aln[keys[i]], aln[keys[j]], 
                                           infile=infile, outfile=outfile)
            distmat[j][i] = distmat[i][j]
    
    
    if madeOutfile:
        os.remove("outfile")
    
    os.remove(infile)
    os.remove(outfile)

    return distmat

'''

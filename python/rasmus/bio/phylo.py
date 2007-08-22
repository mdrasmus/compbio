#
# Phylogeny functions
#
# 


# python imports
import copy
import math
import os
import sys


# rasmus imports
from rasmus import cluster
from rasmus import graph
from rasmus import tablelib
from rasmus import stats
from rasmus import treelib
from rasmus import util

from rasmus.bio import blast
from rasmus.bio import fasta
from rasmus.bio import phylip


# I was thinking about moving these functions over here... not sure yet
from rasmus.genomeutil import gene2species, makeGene2species, readGene2species



#=============================================================================
# Phylogeny functions
#
    

def viewTree(tree, options = "-t 1"):
    tmpfile = util.tempfile(".", "vistree", ".tree")
    tree.write(tmpfile)
    os.system("vistree.py -n %s %s" % (tmpfile, options))
    os.remove(tmpfile)
    

def tree2distmat(tree, leaves):
    """get pair-wise distances between leaves of a tree"""
    
    # TODO: not implemented efficiently
    mat = []
    for i in range(len(leaves)):
        mat.append([])
        for j in range(len(leaves)):
            mat[-1].append(treelib.findDist(tree, leaves[i], leaves[j]))
    
    return mat


def reconcile(gtree, stree, gene2species = gene2species):
    recon = {}
    #depths = stree.findDepths(stree.root)
    
    # label gene leaves with their species
    for node in gtree.leaves():
        recon[node] = stree.nodes[gene2species(node.name)]
    
    # recurse through gene tree
    def walk(node):
        node.recurse(walk)
        
        if not node.isLeaf():
            # this node's species is lca of children species        
            recon[node] = treelib.lca(util.mget(recon, node.children))
    walk(gtree.root)
    
    return recon


def reconcileNode(node, stree, recon):
    return treelib.lca(util.mget(recon, node.children))


def labelEvents(gtree, recon):
    events = {}
    
    def walk(node):
        events[node] = labelEventsNode(node, recon)
        node.recurse(walk)
    walk(gtree.root)
    
    return events

def labelEventsNode(node, recon):
    if not node.isLeaf():
        if recon[node] in map(lambda x: recon[x], node.children):
            return "dup"
        else:
            return "spec"
    else:
        return "gene"

def findLossUnderNode(node, recon):
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


def findLossNode(node, recon):
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
    if labelEventsNode(node.parent, recon) == "dup":
        spath.append(send)
    
    # go up species path (skip starting species)
    # every node on the list is at least one loss
    for i, snode in enumerate(spath[1:]):
        for schild in snode.children:
            if schild != spath[i]:
                loss.append([node, schild])
        
    return loss


def countDupNode(node, events):
    if events[node] == "dup":
        return len(node.children) - 1
    else:
        return 0


def findLoss(gtree, stree, recon, node=None):
    #depths = stree.findDepths()
    loss = []

    def walk(node):
        loss.extend(findLossNode(node, recon))
        node.recurse(walk)
    if node:
        walk(node)
    else:
        walk(gtree.root)

    return loss


def countDup(gtree, events, node=None):    
    var = {"dups": 0}
    
    def walk(node):
        var["dups"] += countDupNode(node, events)
        node.recurse(walk)
    if node:
        walk(node)
    else:
        walk(gtree.root)
    
    return var["dups"]
        

def countDupLoss(gtree, stree, recon, events=None):
    if events is None:
        events = labelEvents(gtree, recon)
    
    nloss = len(findLoss(gtree, stree, recon))
    ndups = countDup(gtree, events)
    return nloss + ndups
    


def reconRoot(gtree, stree, gene2species = gene2species, 
               rootby = "duploss", newCopy=True):
    # make a consistent unrooted copy of gene tree
    if newCopy:
        gtree = gtree.copy()
    treelib.unroot(gtree, newCopy=False)
    treelib.reroot(gtree, 
                   gtree.nodes[util.sort(gtree.leafNames())[0]].parent.name, 
                   onBranch=False, newCopy=False)
    
    
    # make recon root consistent for rerooting tree of the same names
    # TODO: there is the possibility of ties, they are currently broken
    # arbitrarily.  In order to make comparison of reconRooted trees with 
    # same gene names accurate, hashOrdering must be done, for now.
    hashOrderTree(gtree, gene2species)
    
    # get list of edges to root on
    edges = []
    def walk(node):
        edges.append((node, node.parent))
        if not node.isLeaf():
            node.recurse(walk)
            edges.append((node, node.parent))
    for child in gtree.root.children:
        walk(child)
    
    
    # try initial root and recon    
    treelib.reroot(gtree, edges[0][0].name, newCopy=False)
    recon = reconcile(gtree, stree, gene2species)
    events = labelEvents(gtree, recon)     
    
    # find reconciliation that minimizes loss
    minroot = edges[0]
    rootedge = sorted(edges[0])
    if rootby == "dup": 
        cost = countDup(gtree, events)
    elif rootby == "loss":
        cost = len(findLoss(gtree, stree, recon))
    elif rootby == "duploss":
        cost = countDupLoss(gtree, stree, recon, events)
    else:
        raise "unknown rootby value '%s'"  % rootby
    mincost = cost
    
    
    # try rooting on everything
    for edge in edges[1:-1]:
        if sorted(edge) == rootedge:
            continue
        rootedge = sorted(edge)
        
        node1, node2 = edge
        if node1.parent != node2:
            node1, node2 = node2, node1
        assert node1.parent == node2, "%s %s" % (node1.name, node2.name)
        
        # uncount cost
        if rootby in ["dup", "duploss"]:
            if events[gtree.root] == "dup":
                cost -= 1
            if events[node2] == "dup":
                cost -= 1
        if rootby in ["loss", "duploss"]:
            cost -= len(findLossUnderNode(gtree.root, recon))
            cost -= len(findLossUnderNode(node2, recon))
        
        # new root and recon
        treelib.reroot(gtree, node1.name, newCopy=False)        
        
        recon[node2] = reconcileNode(node2, stree, recon)
        recon[gtree.root] = reconcileNode(gtree.root, stree, recon)
        events[node2] = labelEventsNode(node2, recon)
        events[gtree.root] = labelEventsNode(gtree.root, recon)
        
        if rootby in ["dup", "duploss"]:
            if events[node2] ==  "dup":
                cost += 1
            if events[gtree.root] ==  "dup":
                cost += 1
        if rootby in ["loss", "duploss"]:
            cost += len(findLossUnderNode(gtree.root, recon))
            cost += len(findLossUnderNode(node2, recon))
        
        #print edge[0].name, edge[1].name, cost
        
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
        treelib.reroot(gtree, node1.name, newCopy=False)
    
    return gtree


def midrootRecon(tree, stree, recon, events, params, generate):

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





def partitionTree(tree, stree, gene2species):
    recon = reconcile(tree, stree, gene2species)
    sroots = findSpeciesRoots(tree, stree, recon)
    
    # extend subroots
    sroots = treelib.maxDisjointSubtrees(tree, sroots)
    
    # partition
    trees = []
    for sroot in sroots:
        trees.append(treelib.subtree(tree, sroot))
        trees[-1].root.parent = None
    
    return trees


def findSpeciesRoots(tree, stree, recon):
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
        

def findSpeciesSets(stree):
    """
    returns a mapping for each species tree node 
    to a set of modern day species
    """
    
    setmap = {}
    def walk(node):
        node.recurse(walk)
        if node.isLeaf():
            setmap[node] = {node.name:1}
        else:
            setmap[node] = {}
            for child in node.children:
                setmap[node].update(setmap[child])
    walk(stree.root)
    
    return setmap


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
    lookup = cluster.item2part(parts)
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


def stree2gtree(stree, genes, gene2species):
    """Create a gene tree with the same topology as the species tree"""
    
    tree = stree.copy()
    
    for gene in genes:
        tree.rename(gene2species(gene), gene)
    return tree



def findOrthologs(gtree, stree, recon):
    """Find all ortholog pairs within a gene tree"""

    events = labelEvents(gtree, recon)
    orths = []
    
    for node, event in events.items():
        if event == "spec":
            leavesmat = [x.leaves() for x in node.children]
            
            for i in range(len(leavesmat)):
                for j in range(i+1, len(leavesmat)):
                    for gene1 in leavesmat[i]:
                        for gene2 in leavesmat[j]:
                            if gene1.name < gene2.name:
                                orth = (gene1.name, gene2.name)
                            else:
                                orth = (gene2.name, gene1.name)
                            orths.append(orth)
    
    return orths






#=============================================================================
# Tree hashing
#

def hashTreeCompose(childHashes):
    return "(%s)" % ",".join(childHashes)

def hashTree(tree, smap = lambda x: x):
    def walk(node):
        if node.isLeaf():
            return smap(node.name)
        else:
            childHashes = map(walk, node.children)
            childHashes.sort()
            return hashTreeCompose(childHashes)
    
    if isinstance(tree, treelib.Tree) or hasattr(tree, "root"):
        return walk(tree.root)
    elif isinstance(tree, treelib.TreeNode):
        return walk(tree)
    else:
        raise Exception("Expected Tree object")


def hashOrderTree(tree, smap = lambda x: x):
    def walk(node):
        if node.isLeaf():
            return smap(node.name)
        else:
            childHashes = map(walk, node.children)
            ind = util.sortInd(childHashes)
            childHashes = util.mget(childHashes, ind)
            node.children = util.mget(node.children, ind)
            return hashTreeCompose(childHashes)
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


def traverseTreeInorder(tree, func):
    def walk(node):
        if node.isLeaf():
            func(node)
        else:
            assert len(node.children) == 2            
            walk(node.children[0])
            func(node)
            walk(node.children[1])
    walk(tree.root)


def getSpeciesInorder(tree):
    species = []
    def func(node):
        species.append(node.name)
    traverseTreeInorder(tree, func)
    return species



def getBranchLens(trees, stree, gene2species=gene2species):
    # determine species nanes
    species = map(str, stree.nodes.keys())
    species.remove(str(stree.root.name))
    
    # make rates table
    rates = tablelib.Table(headers=species)
    
    # loop through trees
    for tree in trees:
        if isinstance(tree, str):
            tree = treelib.readTree(tree)
        recon = reconcile(tree, stree, gene2species)
        events = labelEvents(tree, recon)
        
        # skip trees with duplications or with extremly long branch lengths
        assert "dup" not in events.values()
        
        row = {}
        for node in tree.nodes.values():
            row[str(recon[node].name)] = node.dist
        rates.append(row)
    
    return rates


def findBranchDistrib(trees, stree, gene2species = gene2species, 
                      relative = True):
    """Older version of getBranchLens()
    
       Will probably be deprecated soon.
    """
    
    lengths = util.Dict(1, [])
    used = []

    for tree in trees:
        #tree = reconRoot(tree, stree, gene2species, newCopy=False)
        recon = reconcile(tree, stree, gene2species)
        events = labelEvents(tree, recon)
        
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
# Phylogenetic reconstruction: Neighbor-Joining
#

def neighborjoin(distmat, genes, usertree=None):
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
            if not node.isLeaf():
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
        parent = treelib.TreeNode(tree.newName())
        tree.addChild(parent, tree.nodes[gene1])
        tree.addChild(parent, tree.nodes[gene2])
        
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
    
    # join the last two genes into a tribranch
    gene1, gene2 = leaves.keys()
    if type(gene1) != int:
        gene1, gene2 = gene2, gene1
    tree.addChild(tree.nodes[gene1], tree.nodes[gene2])
    tree.nodes[gene2].dist = dists[gene1][gene2]
    tree.root = tree.nodes[gene1]

    # root tree according to usertree    
    if usertree != None and treelib.isRooted(usertree):
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

def leastSquareError(tree, distmat, genes, forcePos=True, weighting=False):
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
            
    
    if treelib.isRooted(tree):
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
    topmat, edges = makeTopologyMatrix(tree, genes)
    
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
    setBranchLengths(tree, edges, edgelens, paths, resids, 
                     topmat=topmat, rootedge=rootedge)
    
    return util.Bundle(resids=resids, 
                       paths=paths, 
                       edges=edges, 
                       topmat=topmat)


def makeWeightMatrix(topmat, paths):
    import scipy
    
    weightmat = scipy.transpose(topmat)
    
    
    for row in weightmat:
        tot = sum(row)
        
        for i in xrange(len(row)):
            row[i] *= 1.0 #* float(paths[i])
    
    return weightmat


def makeTopologyMatrix(tree, genes):

    # find how edges split vertices
    network = treelib.tree2graph(tree)
    splits = findAllBranchSplits(network, set(genes))
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


def setBranchLengths(tree, edges, edgelens, paths, resids, 
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
    labels, distmat = phylip.readDistMatrix(outfile)

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



def jukesCantorCorrection(dist):
    """Applies the Jukes Cantor correction to distances
       
       Only valid for distances less than .75 sub/site
    """
    assert (dist < .75)
    return - (3/4.0) * log(1 - (4/3.) * dist)




#============================================================================
# branch splits
#

def findAllBranchSplits(network, leaves):
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


def findBranchSplits(tree):
    splits = findAllBranchSplits(treelib.tree2graph(tree), tree.leafNames())
    splits2 = {}
    
    for edge, sets in splits.iteritems():
        # skip external edges
        if len(sets[0]) == 1 or len(sets[1]) == 1:
            continue
        
        s = tuple(sorted([tuple(sorted(i.keys())) for i in sets]))
        splits2[edge] = s
    
    # if tree is rooted, remove duplicate edge
    if treelib.isRooted(tree):
        edge1 = tuple(sorted([tree.root.name, tree.root.children[0].name]))
        edge2 = tuple(sorted([tree.root.name, tree.root.children[1].name]))
        if edge1 > edge2:
            edge1, edge2 = edge2, edge1
        if edge1 in splits2 and edge2 in splits2:
            del splits2[edge1]
    
    return splits2


def robinsonFouldsError(tree1, tree2):
    splits1 = findBranchSplits(tree1)
    splits2 = findBranchSplits(tree2)

    overlap = set(splits1.values()) & set(splits2.values())
    
    #assert len(splits1) == len(splits2)

    return 1 - (len(overlap) / float(max(len(splits1), len(splits2))))



#============================================================================
# duplication loss counting
#


def initDupLossTree(stree):
    # initalize counts to zero
    def walk(node):
        node.data['dup'] = 0
        node.data['loss'] = 0
        node.data['appear'] = 0
        node.data['genes'] = 0
        node.recurse(walk)
    walk(stree.root)


# count dup loss
def countDupLossTree(tree, stree, gene2species):
    recon = reconcile(tree, stree, gene2species)
    events = labelEvents(tree, recon)
    losses = findLoss(tree, stree, recon)
    
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


# count ancestral genes
def countAncestralGenes(stree):
    def walk(node):
        if not node.isLeaf():
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


def writeEventTree(stree, out=sys.stdout):
    labels = {}
    for name, node in stree.nodes.iteritems():
        labels[name] = "[%s]\nD=%d,L=%d;\nG=%d;" % \
                       (str(name),
                        node.data['dup'], node.data['loss'],
                        node.data['genes'])

    treelib.drawTree(stree, labels=labels, minlen=15, spacing=4, labelOffset=-3,
                     out=out)


#=============================================================================
# Old code 
#

"""

# OLD inefficent version of reconRoot()

def reconRoot2(gtree, stree, gene2species = gene2species, 
              rootby = "duploss"):
    # find reconciliation that minimizes loss
    mincost = util.INF
    minroot = None
    minrecon = None
    
    # make an unrooted copy of gene tree
    # TODO: this can be simplified (root on node.parent)
    gtree = treelib.reroot(gtree, util.sort(gtree.leafNames())[0])
    gtree = treelib.unroot(gtree)
    
    # make recon root consistent for rerooting tree of the same names
    # TODO: there is the possibility of ties, they are currently broken
    # arbitrarily.  In order to make comparison of reconRooted trees with 
    # same gene names accurate, hashOrdering must be done, for now.
    hashOrderTree(gtree, gene2species)
    
    # determine graph and possible roots
    mat = treelib.tree2graph(gtree)
    newroots = util.sort(gtree.nodes.keys())
    newroots.remove(gtree.root.name)

    # try rooting on everything
    for root in newroots:
        gtree2 = treelib.reroot(gtree, root, mat)
        recon = reconcile(gtree2, stree, gene2species)
        
        if rootby == "dup":
            events = labelEvents(gtree2, recon)        
            cost = countDup(gtree2, events)
        elif rootby == "loss":
            cost = len(findLoss(gtree2, stree, recon))
        elif rootby == "duploss":
            cost = countDupLoss(gtree2, stree, recon)
        else:
            raise "unknown rootby value '%s'"  % rootby
        
        # keep track of min loss
        if cost < mincost:
            mincost = cost
            minroot = root
            minrecon = recon
    
    # handle the case where no rerooting was attempted (nleaves == 1)
    if minroot == None:
        return gtree
    
    # root tree by minroot
    return treelib.reroot(gtree, minroot)
"""



#============================================================================
# old orthology stuff
#

"""   
def orthologs(gene, otherGenome, gtree, stree, recon):
    genome = recon[gtree.nodes[gene]]
    genome2 = stree.nodes[otherGenome]
    
    # determine additional labels for trees
    setmap = findSpeciesSets(stree)
    events = labelEvents(gtree, recon)
    
    # ascend until target genome is found
    ptr = gtree.nodes[gene]
    while ptr != gtree.root and \
          genome2.name not in setmap[recon[ptr]]:
        ptr = ptr.parent
    
    # the first node that has genome2 is a duplication node, then
    # there are no orthologs
    if events[ptr] == "dup":
        return []
    
    # orthologs are leaves of ptr that are in otherGenome
    nodes = filter(lambda x: recon[x] == genome2, gtree.leaves(ptr))
    return [x.name for x in nodes]


def paralogs(gene, gtree, recon):
    # go up gene tree until other genomes are encountered
    ptr = gtree.nodes[gene]
    while ptr != gtree.root and recon[ptr.parent] == recon[ptr]:
        ptr = ptr.parent
    
    genes = gtree.leafNames(ptr)
    genes.remove(gene)
    return genes



def findBootTreeHomology(tree, stree, homology, niters, 
                         gene2species=gene2species):
    genomes = stree.leafNames()
    trees = partitionTree(tree, stree)

    util.log("processing tree of %d parts" % len(trees))
    
    for tree2 in trees:
        tree2 = reconRoot(tree2, stree)
        recon = reconcile(tree2, stree)
        
        for gene in tree2.leafNames():
            # count orthologs
            for genome in genomes:
                if genome == gene2species(gene):
                    continue

                ogenes = orthologs(gene, genome, tree2, 
                                             stree, recon)
                for gene2 in ogenes:
                    homology.incOrtholog(gene, gene2, .5 / niters)

            # count paralogs
            pgenes = paralogs(gene, tree2, recon)
            for gene2 in pgenes:
                homology.incParalog(gene, gene2, .5 / niters)
    util.toc()
    
    return trees


def findSimplePartHomology(genes, homology,
                         gene2species=gene2species):
    parts = util.Dict(1, [])
    for gene in genes:
        parts[gene2species(gene)].append(gene)
    
    # count orthologs
    for a in xrange(len(parts)):
        for b in xrange(a+1, len(parts)):
            for gene1 in parts.values()[a]:
                for gene2 in parts.values()[b]:
                    homology.incOrtholog(gene1, gene2, 1)
    
    # count paralogs
    for part in parts.values():
        for a in xrange(len(part)):
            for b in xrange(a+1, len(part)):
                homology.incParalog(part[a], part[b], 1)


def homology2orthologSets(homology, cutoff = 0):
    mat = util.Dict(2)
    genes = homology.getGenes()
    
    for gene1 in genes:
        for gene2 in homology.getOrthologs(gene1):
            if homology.getOrthologBootstrap(gene1, gene2) > cutoff:
                mat[gene1][gene2] = 1
        for gene2 in homology.getParalogs(gene1):
            if homology.getParalogBootstrap(gene1, gene2) > cutoff:
                mat[gene1][gene2] = 1
    
    comps = graph.connectedComponents(genes, lambda x: mat[x].keys())
    comps = filter(lambda x: len(x) > 1, comps)
    
    return comps


class Homology:
    def __init__(self):
        self.orths = {}
        self.paras = {}
    
    def incOrtholog(self, gene1, gene2, inc=1):
        val = self.orths.setdefault(gene1, {}).setdefault(gene2, 0)
        self.orths[gene1][gene2] = val + inc
    
        val = self.orths.setdefault(gene2, {}).setdefault(gene1, 0)
        self.orths[gene2][gene1] = val + inc

    
    def incParalog(self, gene1, gene2, inc=1):
        val = self.paras.setdefault(gene1, {}).setdefault(gene2, 0)
        self.paras[gene1][gene2] = val + inc
        
        val = self.paras.setdefault(gene2, {}).setdefault(gene1, 0)
        self.paras[gene2][gene1] = val + inc
    
    def addOrtholog(self, gene1, gene2, bootstrap=1.0):
        self.orths.setdefault(gene1, {})[gene2] = bootstrap
        self.orths.setdefault(gene2, {})[gene1] = bootstrap
    
    def addParalog(self, gene1, gene2, bootstrap=1.0):
        self.paras.setdefault(gene1, {})[gene2] = bootstrap
        self.paras.setdefault(gene2, {})[gene1] = bootstrap
        
    def getOrthologs(self, gene1):
        if gene1 in self.orths:
            return self.orths[gene1].keys()
        else:
            return []

    def getParalogs(self, gene1):
        if gene1 in self.paras:
            return self.paras[gene1].keys()
        else:
            return []
    
    def getOrthologBootstrap(self, gene1, gene2):
        return self.orths[gene1][gene2]
    
    def getParalogBootstrap(self, gene1, gene2):
        return self.paras[gene1][gene2]
    
    def getGenes(self):
        return util.unique(self.orths.keys() + self.paras.keys())
    
    
    def getParts(self, cutoff=0):
        def getNeighbors(gene):
            return filter(lambda x: self.getOrthologBootstrap(gene, x) >= cutoff,
                          self.getOrthologs(gene)) + \
                   filter(lambda x: self.getParalogBootstrap(gene, x) >= cutoff,
                          self.getParalogs(gene))
        
        comps = graph.connectedComponents(self.getGenes(), getNeighbors)
        
        return comps
    
    def write(self, out = sys.stdout):
        genes = self.getGenes()
        genes.sort()
        
        print >>out, "<?xml version='1.0' encoding='ISO-8859-1'?>"
        print >>out, "<homology>"

        for gene1 in genes:
            print >>out, "<gene name='%s'>" % gene1

            keys = self.getOrthologs(gene1)
            keys.sort()
            for gene2 in keys:
                print >>out, "  <homolog name='%s' type='ortholog' bootstrap='%f' />" % \
                    (gene2, self.getOrthologBootstrap(gene1, gene2))

            keys = self.getParalogs(gene1)
            keys.sort()
            for gene2 in keys:
                print >>out, "  <homolog name='%s' type='paralog' bootstrap='%f' />" % \
                    (gene2, self.getParalogBootstrap(gene1, gene2))
            print >>out, "</gene>"
            print >>out
        
        print >>out, "</homology>"
    

    def read(self, filename):
        infile = util.openStream(filename)

        # Create a parser
        parser = make_parser()

        # Tell the parser we are not interested in XML namespaces
        #parser.setFeature(feature_namespaces, 0)

        # Tell the parser to use our handler
        handler = HomologyHandler(self)
        parser.setContentHandler(handler)

        # Parse the input
        parser.parse(infile)



class HomologyHandler(xml.sax.handler.ContentHandler):
    def __init__(self, homology):
        self.homology = homology
        self.gene1 = None
        self.boostrap = None
        self.kind = None
        self.elm = ""
    
    def startElement(self, name, attrs):
        if name == "gene":
            self.gene1 = str(attrs["name"])
        elif name == "homolog":
            if attrs["type"] == "ortholog":
                self.homology.addOrtholog(self.gene1, str(attrs["name"]), 
                                          float(attrs["bootstrap"]))
            elif attrs["type"] == "paralog":
                self.homology.addParalog(self.gene1, str(attrs["name"]), 
                                         float(attrs["bootstrap"]))
        self.elm = name
    
    def endElement(self, name):
        self.elm = ""
    
    def characters(self, text):
        pass

"""


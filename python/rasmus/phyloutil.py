# python imports
import copy
import math
import os
import sys

# xml support
from xml.sax import make_parser
import xml.sax.handler

# rasmus imports
import blast
import cluster
import ensembl
import genomeutil
import graph
import util
import treelib




def viewTree(tree, options = "-t 1"):
    tmpfile = util.tempfile(".", "vistree", ".tree")
    tree.write(tmpfile)
    os.system("vistree.py -n %s %s" % (tmpfile, options))
    os.remove(tmpfile)
    


def reconcile(gtree, stree, gene2species = genomeutil.gene2species):
    recon = {}
    depths = stree.findDepths(stree.root)
    
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


def findLossNode(node, recon):
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


def countDupNode(node, events):
    if events[node] == "dup":
        return len(node.children) - 1
    else:
        return 0


def findLoss(gtree, stree, recon, node=None):
    depths = stree.findDepths()
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
    


def reconRoot(gtree, stree, gene2species = genomeutil.gene2species, 
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
            cost -= len(findLossNode(gtree.root, recon))
            cost -= len(findLossNode(node2, recon))
        
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
            cost += len(findLossNode(gtree.root, recon))
            cost += len(findLossNode(node2, recon))
        
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


def reconRoot2(gtree, stree, gene2species = genomeutil.gene2species, 
              rootby = "duploss"):
    # find reconciliation that minimizes loss
    mincost = util.INF
    minroot = None
    minrecon = None
    
    # make an unrooted copy of gene tree
    # TODO: this can be simplified (root on node.parent)
    gtree = treelib.reroot(gtree, util.sort(gtree.leaveNames())[0])
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
    return walk(tree.root)

def hashOrderTree(tree, smap = lambda x: x):
    def walk(node):
        if node.isLeaf():
            return smap(node.name)
        else:
            childHashes = map(walk, node.children)
            ind = util.sortInd(childHashes)
            childHashes = util.permute(childHashes, ind)
            node.children = util.permute(node.children, ind)
            return hashTreeCompose(childHashes)
    walk(tree.root)



def findBranchDistrib(trees, stree, gene2species = genomeutil.gene2species, 
                      relative = True):
    lengths = util.Dict(1, [])
    trees2 = []

    for tree in trees:
        tree = reconRoot(tree, stree, gene2species)
        recon = reconcile(tree, stree, gene2species)
        events = labelEvents(tree, recon)
        
        # skip trees with duplications or with extremly long branch lengths
        if "dup" in events.values() or \
            max(x.dist for x in tree.nodes.values()) > 2:
            continue
        else:
            trees2.append(tree)
        
        # find total length of tree
        totalLength = 0
        for node in tree.nodes.values():
            totalLength += node.dist
        
        for node in tree.nodes.values():
            if relative:
                lengths[recon[node]].append(node.dist/totalLength)
            else:                
                lengths[recon[node]].append(node.dist)
    
    
    return trees2, lengths


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
        
        # dont count hits within a cluster
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
    tree = stree.copy()
    
    for gene in genes:
        tree.rename(gene2species(gene), gene)
    return tree




def findAllOrthologs(gtree, stree, recon):
    events = labelEvents(gtree, recon)
    orths = []
    
    for node, event in events.items():
        if event == "spec":
            leavesmat = [x.leaves() for x in node.children]
            
            for i in range(len(leavesmat)):
                for j in range(i+1, len(leavesmat)):
                    for gene1 in leavesmat[i]:
                        for gene2 in leavesmat[j]:
                            orths.append(util.sort([gene1.name, gene2.name]))
    
    return orths




def neighborjoin(distmat, genes):
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
    
    # join loop
    while len(leaves) > 2:
        # search for closest genes
        low = util.INF
        lowpair = (None, None)
        leaveslst = leaves.keys()

        for i in range(len(leaves)):
            for j in range(i+1, len(leaves)):
                gene1, gene2 = leaveslst[i], leaveslst[j]
                dist = dists[gene1][gene2] - restdists[gene1] - restdists[gene2]
                
                if dist < low:
                    low = dist
                    lowpair = (gene1, gene2)
        
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
    
    return tree



def neighborjoinSparse(dists):
    """dists must be symmetric"""
    
    tree = treelib.Tree()
    leaves = {}
    restdists = {}
    genes = dists.keys()
    
    
    # initialize rest distances
    for gene1 in genes:
        r = 0
        n = 0
        for gene2 in genes:
            edges = dists[gene1]
            if gene2 in edges:
                r += edges[gene2]
                n += 1
        restdists[gene1] = r / float(max(n-2, 1))
    
    
    # initialize leaves
    for gene in genes:
        tree.add(treelib.TreeNode(gene))
        leaves[gene] = 1
    
    # join loop
    while True:
        # search for closest genes
        low = util.INF
        lowpair = (None, None)

        for gene1 in leaves:
            for gene2 in dists[gene1]:
                if gene2 not in leaves:
                    continue
                if gene1 in dists and gene2 in dists[gene1]:
                    dist = dists[gene1][gene2] - restdists[gene1] \
                                               - restdists[gene2]
                    
                    if dist < low:
                        low = dist
                        lowpair = (gene1, gene2)
        
        #print lowpair
        
        if lowpair == (None, None):
            # split partial tree into mutliple trees
            trees = []
            for leaf in leaves:
                trees.append(treelib.subtree(tree, tree.nodes[leaf]))
            
            return trees
        elif len(leaves) == 2:
            break
        
        # join gene1 and gene2
        gene1, gene2 = lowpair
        parent = treelib.TreeNode(tree.newName())
        tree.addChild(parent, tree.nodes[gene1])
        tree.addChild(parent, tree.nodes[gene2])
        
        
        
        # set distances
        tree.nodes[gene1].dist = (dists[gene1][gene2] + restdists[gene1] - 
                                  restdists[gene2]) / 2.0
        tree.nodes[gene2].dist = dists[gene1][gene2] - tree.nodes[gene1].dist
        
        #print "----------"
        #for leaf in leaves:
        #    tree.writeNewickNode(tree.nodes[leaf])
        
        # gene1 and gene2 are no longer leaves
        del leaves[gene1]
        del leaves[gene2]
        
        gene3 = parent.name
        r = 0
        n = 0
        for gene in leaves:
            if gene1 in dists and \
               gene in dists[gene1] and \
               gene2 in dists and \
               gene in dists[gene2]:
                d = (dists[gene1][gene] + dists[gene2][gene] -
                     dists[gene1][gene2]) / 2.0
                dists[gene3][gene] = d
                dists[gene][gene3] = d
                r += d
                n += 1
        leaves[gene3] = 1
        
        restdists[gene3] = r / float(max(n - 2, 1))
    
    # join the last two genes into a tribranch
    gene1, gene2 = leaves.keys()
    if type(gene1) != int:
        gene1, gene2 = gene2, gene1
    tree.addChild(tree.nodes[gene1], tree.nodes[gene2])
    tree.nodes[gene2].dist = dists[gene1][gene2]
    tree.root = tree.nodes[gene1]
    
    return [tree]



#============================================================================
# duplication loss counting
#============================================================================


def initDupLossTree(stree):
    # initalize counts to zero
    def walk(node):
        node.data['dup'] = 0
        node.data['loss'] = 0
        node.data['genes'] = 0
        node.recurse(walk)
    walk(stree.root)


# count dup loss
def countDupLossTree(tree, stree, gene2species):
    reconRoot(tree, stree, gene2species, newCopy=False)
    recon = reconcile(tree, stree, gene2species)
    events = labelEvents(tree, recon)
    losses = findLoss(tree, stree, recon)

    # count dups
    for node, event in events.iteritems():
        if event == "dup":
            recon[node].data['dup'] += 1
        elif event == "gene":
            recon[node].data['genes'] += 1

    # count losses
    for gnode, snode in losses:
        snode.data['loss'] += 1


# count ancestral genes
def countAncestralGenes(stree):
    def walk(node):
        if not node.isLeaf():
            counts = []
            for child in node.children:
                walk(child)
                counts.append(child.data['genes'] -
                              child.data['dup'] + child.data['loss'])
            assert util.equal(* counts), str(counts)
            node.data['genes'] = counts[0]
    walk(stree.root)



#============================================================================
# old orthology stuff
#============================================================================
            
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
    
    genes = gtree.leaveNames(ptr)
    genes.remove(gene)
    return genes



def findBootTreeHomology(tree, stree, homology, niters, 
                         gene2species=genomeutil.gene2species):
    genomes = stree.leaveNames()
    trees = partitionTree(tree, stree)

    util.log("processing tree of %d parts" % len(trees))
    
    for tree2 in trees:
        tree2 = reconRoot(tree2, stree)
        recon = reconcile(tree2, stree)
        
        for gene in tree2.leaveNames():
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
                         gene2species=genomeutil.gene2species):
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




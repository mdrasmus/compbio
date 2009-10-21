import random
import math
import StringIO


from rasmus import treelib
from rasmus import stats
from rasmus import matrixlib
from rasmus.bio import phylo


import Spidir
from Spidir.Debug import *

#=============================================================================
# Tree search
#


def proposeNni(tree, node1, node2, change=0):
    """Proposes a new tree using Nearest Neighbor Interchange
       
       Branch for NNI is specified by giving its two incident nodes (node1 and 
       node2).  Change specifies which  subtree of node1 will be swapped with
       the uncle.  See figure below.

         node2
        /     \
      uncle    node1
               /  \
         child[0]  child[1]
    
    special case with rooted branch:
    
              node2
             /     \
        node2'      node1
       /     \     /     \
      uncle   * child[0] child[1]
    
    """
    
    if node1.parent != node2:
        node1, node2 = node2, node1  
    
    # try to see if edge is one branch (not root edge)
    if treelib.is_rooted(tree) and \
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
        assert node1.parent == node2, debugBranch(tree, node1, node2)
    
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


def edge2branch(tree, node1, node2):
    # test for root edge
    if node1.parent == node2.parent == tree.root:
        return node1, node2
    else:
        # ensure node2 is parent of node1
        if node1.parent != node2:
            node1, node2 = node2, node1
        assert node1.parent == node2, debugBranch(tree, node1, node2)
        
        return node1, node2

def isRootEdge(tree, node1, node2):
    return node1.parent == node2.parent == tree.root


def getNniUncle(tree, node1, node2):
    print >>sys.stderr, "n>>", node1.name, node2.name

    if node1.parent != node2:
        node1, node2 = node2, node1

    if node1.parent == node2.parent == tree.root:
        # root edge
        return node2.children[0]
    else:
        assert node1.parent == node2, debugBranch(tree, node1, node2)
    
        # find uncle
        uncle = 0 
        if node2.children[uncle] == node1:
            uncle = 1

        return node2.children[uncle]


def debugBranch(tree, node1, node2):
    treelib.drawTreeNames(tree, minlen=5, maxlen=5, out=sys.stderr)
    print >>sys.stderr, "branch:", node1.name, node2.name
    


def splitTree(tree, node1, node2):
    """Distroys original tree"""

    # ensure node2 is parent of node1
    if node1.parent != node2:
        node1, node2 = node2, node1
    assert node1.parent == node2
        
    # create new tree object
    tree2 = treelib.Tree(nextname = tree.new_name())
    tree2.root = node1  
    tree2.copyData(tree)
    
    # add nodes
    def walk(node):
        tree2.add(node)
        node.recurse(walk)
    walk(tree2.root)

    # break subtree at node1 off of original tree
    tree1 = tree
    tree1.removeTree(node1)
    
    # remove single child node
    if node2 != tree1.root:
        node3 = node2.parent
        tree1.remove(node2)
        for n in node2.children:
            tree1.add_child(node3, n)
    else:
        # special case
        assert len(node2.children) == 1
        tree1.remove(node2)
        tree1.root = node2.children[0]
        tree1.root.parent = None
            
    tree2.root.parent = None
    
    #treelib.drawTreeNames(tree1, minlen=5, maxlen=5)
    #treelib.drawTreeNames(tree2, minlen=5, maxlen=5)
    
    #treelib.assertTree(tree1)
    #treelib.assertTree(tree2)
    
    return tree1, tree2


def regraftTree(tree, subtree, node1, node2):
    """
    Add a subtree to an edge in an existing tree
    Note: Tree distances become invalid after add
    """

    # ensure node2 is parent of node1
    if node1.parent != node2:
        node1, node2 = node2, node1
    assert node1.parent == node2    

    # add new node on edge
    newNode = treelib.TreeNode(tree.new_name())
    tree.add_child(newNode, node1)
    node2.children[node2.children.index(node1)] = newNode
    newNode.parent = node2
    
    tree.addTree(newNode, subtree)
    
    #debug("regraft>>")
    #treelib.drawTreeNames(tree, minlen=5, maxlen=5)

    #treelib.assertTree(tree)

    

def proposeTree(conf, tree):
    tree2 = tree.copy()
    
    #if random.random() < conf["rerootprob"]:
    #    nodes = tree.nodes.values()
    #    newnode = nodes[random.randint(0, len(nodes)-1)]
    #    tree2 = treelib.reroot(tree2, newnode.name)
    
    # find edges for NNI
    nodes = tree2.nodes.values()
    nodes = filter(lambda x: not x.isLeaf() and 
                             x != tree2.root, nodes)
    edges = [(node, node.parent) for node in nodes]
    edge = edges[int(random.random() * len(edges))]
    
    proposeNni(tree2, edge[0], edge[1], int(round(random.random())))
    return tree2


def proposeTreeWeighted(tree):
    """Nodes in tree must have logl in their data dict"""
    
    tree2 = tree.copy()
    
    # find edges for NNI
    nodes = tree2.nodes.values()
    nodes = filter(lambda x: not x.isLeaf() and 
                             x != tree2.root, nodes)
    edges = [(node, node.parent) for node in nodes]
    
    # create weights
    weights = []    
    for edge in edges:
        weights.append(edge[0].data["error"])
    
    if sum(weights) == 0:
        l = float(len(weights))
        weights = [1./l for x in weights]
    
    # sample by weight
    edge = edges[stats.sample(weights)]
    
    proposeNni(tree2, edge[0], edge[1], int(round(random.random())))
    return tree2


def proposeTree2(conf, tree,  distmat, labels, 
                  stree, gene2species, params, visited):
    tree2 = tree
    nniprobs = [.70, .20, .10]
    
    pick = stats.sample(nniprobs)
    
    for i in xrange(pick+1):
        tree2 = proposeTree(conf, tree2)

    if random.random() < conf["rerootprob"]:
        phylo.reconRoot(tree2, stree, gene2species, newCopy=False)
    return tree2


def proposeTree3(conf, tree,  distmat, labels, 
                  stree, gene2species, params, visited):
    toplogl = tree.data["logl"]
    toptree = tree.copy()
    
    tree = tree.copy()
    
    nodes = tree.nodes.values()
    nodes.remove(tree.root)
    weights = [1 for x in nodes] #[x.data["error"] for x in nodes]
    badgene = nodes[stats.sample(weights)]
    
    
    # detemine distance from badgene to everyone else
    dists = util.Dict(default=-util.INF)
    def walk(node, dist):
        dists[node.name] = dist
        for child in node.children:
            walk(child, dist + child.dist)
    walk(badgene, 0)
    seen = set([badgene])
    node = badgene.parent
    dist = badgene.dist
    while node != None:        
        for child in node.children:
            if child not in seen:
                walk(child, dist)
        seen.add(node)
        dist +=  node.dist
        node = node.parent
    
    tree1, tree2 = splitTree(tree, badgene, badgene.parent)
    
    names = tree1.nodes.keys()
    names.remove(tree1.root.name)
    #names.sort(key=lambda x: dists[x])
    random.shuffle(names)
    
    
    for name in names[:min(len(names), conf["regraftloop"])]:
        tree = tree1.copy()
        node = tree.nodes[name]
        
        #print "p3>>", node.name, node.parent.name
        regraftTree(tree, tree2.copy(), node, node.parent)
        
        thash = phylo.hash_tree(tree)
        
        if thash not in visited:        
            Spidir.setTreeDistances(conf, tree, distmat, labels)
            logl = Spidir.treeLogLikelihood(conf, tree, stree, gene2species, params)
        addVisited(conf, visited, tree, gene2species, thash)
        logl, tree, count = visited[thash]
        
        if logl > toplogl:
            toplogl = logl
            toptree = tree
            
            # try returning immediately
            #return toptree

    
    assert toptree != None
    
    return toptree



def printMCMC(conf, i, tree, stree, gene2species, visited):
    if isDebug(DEBUG_LOW):
        recon = phylo.reconcile(tree, stree, gene2species)
        events = phylo.labelEvents(tree, recon)            

        debug("\n=======================================")
        debug("iter:", i, " visited:", len(visited))
        drawTreeLogl(tree, events=events)
        debug()
        debug()




def addVisited(conf, visited, tree, gene2species, thash=None):
    if thash is None:
        thash = phylo.hash_tree(tree)
    
    if thash in visited:
        visited[thash][2] += 1
    else:
        visited[thash] = [tree.data["logl"], tree.copy(), 1]

    if "correcthash" in conf:
        if thash == conf["correcthash"]:
            debug("PROPOSED CORRECT TREE: visisted = ", len(visited))
            if conf["searchtest"]:
                drawTreeLogl(tree)
                sys.exit(0)
        
    if "debugtab_file" in conf:
        shash = phylo.hash_tree(tree, gene2species)

        if "correcthash" in conf:
            correct = (conf["correcthash"] == thash)
        else:
            correct = False

        conf["debugtab"].writeRow(conf["debugtab_file"],
                              {"correct": correct,
                               "logl": tree.data["logl"], 
                               "treelen": sum(x.dist for x in tree), 
                               "baserate": tree.data["baserate"], 
                               "error": tree.data["error"], 
                               "errorlogl": tree.data["errorlogl"],
                               "eventlogl": tree.data["eventlogl"], 
                               "tree": tree.getOnelineNewick(),
                               "topology": thash,
                               "species_hash": shash})



def searchHillClimb(conf, distmat, labels, stree, gene2species, params,
               initTree=None, visited=None):

    if visited == None:
        visited = {}
    
    # init with NJ    
    if initTree != None:
        tree = initTree
    else:
        #tree = bionj.bionj(labels=labels, distmat=distmat, verbose=False)
        tree = phylo.neighborjoin(distmat, labels)
        tree = phylo.reconRoot(tree, stree, gene2species)
        Spidir.setTreeDistances(conf, tree, distmat, labels)

    # init likelihood score
    logl = treeLogLikelihood(conf, tree, stree, gene2species, params)

    # store tree in visited
    addVisited(conf, visited, tree, gene2species)
    
    stuck = False
        
    for i in range(conf["hilliters"]):
        printMCMC(conf, i, tree, stree, gene2species, visited)
        
        proposals = getProposals(conf, tree, distmat, labels, 
                                 stree, gene2species, params, visited, stuck)
        
        util.printcols(map(lambda (a,(b,c),d): [a, b.name, c.name, d], proposals))
        print
        
        # determine which proposals to use
        edgeset = set()
        proposals2 = []
        for logl2, edge, change in proposals:
            if edge in edgeset:
                continue
            proposals2.append([logl2, edge, change])
            
            edgeset.add((getNniUncle(tree, edge[0], edge[1]), edge[1]))
            edgeset.add((edge[0].children[change], edge[0]))
            edgeset.add((edge[0], edge[1]))
        
        util.printcols(map(lambda (a,(b,c),d): [a, b.name, c.name, d], proposals2))
        print
        
        heat = 1.0
        start = 0
        while start < len(proposals2):
            nproposals = int(math.ceil(len(proposals2) * heat))
            
            # apply proposals
            for logl3, edge, change in proposals2[start:start+nproposals]:
                proposeNni(tree, edge[0], edge[1], change)
            tree2 = phylo.reconRoot(tree, stree, gene2species, newCopy=True)
            
            # calc likelihood
            thash = phylo.hash_tree(tree2)
            if thash not in visited:
                Spidir.setTreeDistances(conf, tree2, distmat, labels)
                logl2 = Spidir.treeLogLikelihood(conf, tree2, stree, 
                                          gene2species, params)
                stuck = False
            else:
                logl2 = visited[thash][0]
                
                Spidir.setTreeDistances(conf, tree2, distmat, labels)
                logl2 = Spidir.treeLogLikelihood(conf, tree2, stree, 
                                          gene2species, params)
                
                if nproposals == 1:
                    stuck = True
            
            addVisited(conf, visited, tree2, gene2species, thash)
            
            
            debug("logl2", logl2)
            
            if logl2 > logl:
                logl = logl2
                tree = tree2
                break
            
            if nproposals == 1:
                logl = logl2
                tree = tree2
                break
            
            heat *= .5
            
            # undo reversals
            for logl3, edge, change in util.reverse(proposals2[start:start+nproposals]):
                proposeNni(tree, edge[0], edge[1], change)
        
        debug("start:", start)
        debug("swaps:", nproposals)
        debug("heat:", heat)
        debug("stuck:", stuck)

    
    items = visited.items()
    i = util.argmaxfunc(lambda x: x[1][0], items)
    thash, (logl, tree, count) = items[i]
    return tree, logl

        
    


def getProposals(conf, tree, distmat, labels, stree, 
                 gene2species, params, visited, stuck=False):
                 
    # TODO: handle root edges
    
    # try all NNI
    # find edges for NNI
    nodes = tree.nodes.values()
    nodes = filter(lambda x: not x.isLeaf() and 
                             x != tree.root and
                             x not in tree.root.children, nodes)
    edges = [(node, node.parent) for node in nodes]
    edges.append(tuple(tree.root.children))
    
    
    treelib.drawTreeNames(tree, minlen=5, maxlen=5, out=sys.stderr)
    util.printcols(util.map2(lambda x: x.name, edges), out=sys.stderr)
    
    proposals = []
    for edge in edges:
        for change in (0,1):
            proposeNni(tree, edge[0], edge[1], change)
            tree2 = phylo.reconRoot(tree, stree, gene2species, newCopy=True)
            
            thash = phylo.hash_tree(tree2)
            if thash not in visited:
                Spidir.setTreeDistances(conf, tree2, distmat, labels)
                logl = Spidir.treeLogLikelihood(conf, tree2, stree, 
                                         gene2species, params)
                visited[thash] = [logl, tree2, 1]
                
                proposals.append([logl, edge, change])
            else:
                visited[thash][2] += 1
                logl = visited[thash][0]
                
                if not stuck:
                    proposals.append([logl, edge, change])
            
            
            
            # switch branch back
            proposeNni(tree, edge[0], edge[1], change)
    
    proposals.sort(key=lambda x: x[0], reverse=True)
    return proposals


class McmcChain:
    def __init__(self, name, state, logl, propose):
        self.name = name
        self.state = state
        self.logl = logl
        self.propose = propose
        self.relax = 0
    
    
    def step(self):
        nextState, nextLogl = self.propose(self, self.state)

        # accept/reject
        if nextLogl > self.logl or \
           nextLogl - self.logl > log(random.random()) - self.relax:
            # accept new state
            self.state = nextState
            self.logl = nextLogl



def searchRegraft(conf, distmat, labels, stree, gene2species, params,
                  initTree=None, visited=None, proposeFunc=proposeTree2):
    if visited == None:
        visited = {}
    
    # init with NJ    
    if initTree != None:
        tree = initTree
    else:
        tree = phylo.neighborjoin(distmat, labels)
        tree = phylo.reconRoot(tree, stree, gene2species)
        Spidir.setTreeDistances(conf, tree, distmat, labels)

    # init likelihood score
    logl = Spidir.treeLogLikelihood(conf, tree, stree, gene2species, params)
    
    # store tree in visited
    addVisited(conf, visited, tree, gene2species)
    
    # show initial tree
    printMCMC(conf, 0, tree, stree, gene2species, visited)    
    
    
    for i in xrange(conf["regrafts"]):
        tree = proposeTree3(conf, tree,  distmat, labels, 
                               stree, gene2species, params, visited)
        if tree.data["logl"] > logl:
            printMCMC(conf, i, tree, stree, gene2species, visited)
            logl = tree.data["logl"]
    
    return tree, tree.data["logl"]
    
    

def searchMCMC(conf, distmat, labels, stree, gene2species, params,
               initTree=None, visited=None, proposeFunc=proposeTree2):
    if visited == None:
        visited = {}
    
    
    this = util.Bundle(
        nold=0,
        toplogl = -util.INF,
        toptree = None,
        iter=0)
    
    
    # init with NJ    
    if initTree != None:
        tree = initTree
    else:
        tree = phylo.neighborjoin(distmat, labels)
        tree = phylo.reconRoot(tree, stree, gene2species)
        Spidir.setTreeDistances(conf, tree, distmat, labels)

    # init likelihood score
    this.toplogl = Spidir.treeLogLikelihood(conf, tree, stree, gene2species, params)
    this.toptree = tree
    
    # store tree in visited
    addVisited(conf, visited, tree, gene2species)
    
    # show initial tree
    printMCMC(conf, 0, tree, stree, gene2species, visited)
    
    
    # proposal function
    def propose(chain, tree):
        tree2 = proposeFunc(conf, tree,  distmat, labels, 
                            stree, gene2species, params, visited)
        
        # check visited dict
        thash = phylo.hash_tree(tree2)
        if thash in visited:
            logl, tree2, count = visited[thash]
            #this.nold += 1
        else:
            Spidir.setTreeDistances(conf, tree2, distmat, labels)
            logl = Spidir.treeLogLikelihood(conf, tree2, stree, gene2species, params)
            this.nold = 0
        
        addVisited(conf, visited, tree2, gene2species, thash)
        
        # best yet tree
        if logl > this.toplogl:
            printMCMC(conf, "%d:%d" % (chain.name, this.iter), 
                      tree2, stree, gene2species, visited)
            this.toplogl = logl
            this.toptree = tree2.copy()
            
            # move some other chains to best state
            #chains2 = sorted(chains, key=lambda x: x.logl)
            #for chain in chains2[:1]:
            #    chain.state = this.toptree.copy()
            #    chain.logl = this.toplogl
        
        # alter logl to influence search only
        #chain.relax = conf["speedup"] * this.nold
               
        return tree2, logl
        
    # init chains    
    chains = []
    for i in range(conf["nchains"]):
        chains.append(McmcChain(i, tree.copy(), this.toplogl, propose))
    
    
    # run chains
    for i in xrange(1, conf["iters"]):
        this.iter += 1
        
        for chain in chains:
            chain.step()   
   

    return this.toptree, this.toplogl



def searchGreedy(conf, distmat, labels, stree, gene2species, params, visited=None):
    if visited == None:
        visited = {}

    totalgenes = len(labels)
    ngenes = 2
    
    # create initial 2 gene tree (labels[0], labels[1])
    tree = treelib.Tree()
    tree.make_root()
    tree.add_child(tree.root, treelib.TreeNode(labels[0]))
    tree.add_child(tree.root, treelib.TreeNode(labels[1]))
    
    
    for ngenes in xrange(2, totalgenes):
        debug("adding", labels[ngenes])
        
        toplogl = -util.INF
        toptree = None
        
        distmat2 = matrixlib.submatrix(distmat, range(ngenes+1), range(ngenes+1))
        labels2  = labels[:ngenes+1]
        
        
        # place new gene on every branch
        for name in tree.nodes:
            tree2 = tree.copy()
            node = tree2.nodes[name]

            if node == tree2.root:
                newnode = treelib.TreeNode(tree2.new_name())
                tree2.add_child(newnode, tree2.root)
                tree2.root = newnode
                tree2.add_child(newnode, treelib.TreeNode(labels[ngenes]))
            else:
                parent = node.parent
                tree2.remove(node)
                newnode = treelib.TreeNode(tree2.new_name())
                tree2.add_child(parent, newnode)
                tree2.add_child(newnode, node)
                tree2.add_child(newnode, treelib.TreeNode(labels[ngenes]))
            
            #tree2 = phylo.reconRoot(tree2, stree, gene2species)
            Spidir.setTreeDistances(conf, tree2, distmat2, labels2)
            logl = Spidir.treeLogLikelihood(conf, tree2, stree, gene2species, params)

            if logl >= toplogl:
                toplogl = logl
                toptree = tree2
        tree = toptree

        # only use visited hash table if all genes are present        
        if ngenes == totalgenes:
            visited2 = visited
        else:
            # otherwise use a new temp hash table
            visited2 = {}
        
        tree, logl = searchExhaustive(conf, distmat2, labels2, 
                                      tree, stree, gene2species, params,
                                      visited=visited2,
                                      depth=conf["depth"])
            
            
        if logl >= toplogl:
            toplogl = logl
            toptree = tree
        tree = toptree
        
        
        debug()
    
    visited.update(visited2)
    
    return tree, toplogl



def searchExhaustive(conf, distmat, labels, tree, stree, gene2species, params,
                     depth=2, visited=None, visited2=None, topDepth=True,
                     toplogl=None, short=False):
    if visited == None:
        visited = {}
    if visited2 == None:
        visited2 = {}
    
    tree = tree.copy()
    
    # find initial logl
    thash = phylo.hash_tree(tree)
    if thash not in visited:
        Spidir.setTreeDistances(conf, tree, distmat, labels)
        logl = Spidir.treeLogLikelihood(conf, tree, stree, 
                                    gene2species, params)
        visited[thash] = [logl, tree.copy(), 1]
        
        drawTreeLogl(tree)
    else:
        logl = visited[thash][0]
        
    if toplogl == None:
        toplogl = [logl]
    
    
    debug(" " * (depth*2), "(%d)" % len(visited))
    sys.stdout.flush()
    
    if depth < 1:
        return tree, logl
    
    
    # try all NNI
    # find edges for NNI
    nodes = tree.nodes.values()
    nodes = filter(lambda x: not x.isLeaf() and 
                             x != tree.root and \
                             x.parent != tree.root, nodes)
    edges = [(node, node.parent) for node in nodes]
    
    for edge in edges:
        for change in (0,1):
            proposeNni(tree, edge[0], edge[1], change)
            
            thash = phylo.hash_tree(tree)
            if thash not in visited:
                Spidir.setTreeDistances(conf, tree, distmat, labels)
                logl = Spidir.treeLogLikelihood(conf, tree, stree, 
                                         gene2species, params)
                visited[thash] = [logl, tree.copy(), 1]
            else:
                logl = visited[thash][0]
            
            if logl > toplogl[0]:
                toplogl[0] = logl
                
                if short:
                    return tree, logl
                else:
                    printMCMC(conf, "N/A", tree, stree, gene2species, visited)
                
            
            if (thash not in visited2 or \
                depth > visited2[thash]) and \
                logl - toplogl[0] >= conf["eprune"]:
                visited2[thash] = depth
                
                # dig deeper
                if depth > 1:
                    tree2, logl2 = searchExhaustive(conf, distmat, labels, 
                                     tree, stree, gene2species, params,
                                     depth=depth-1, visited=visited,
                                     visited2=visited2,
                                     topDepth=False,
                                     toplogl=toplogl, short=short)
                    
                    if short and tree2 != None:
                        return tree2, logl2
                    
            
            # switch branch back
            proposeNni(tree, edge[0], edge[1], change)
    
    # debug
    if topDepth:
        items = visited.items()
        i = util.argmaxfunc(lambda x: x[1][0], items)
        
        thash, (logl, tree, count) = items[i]
        
        return tree, logl
    else:
        return None, None

#!/usr/bin/env python

# python libs
import sys, time, os, random, math

# rasmus libs
from rasmus import env, treelib
from rasmus import util
from rasmus import stats
import Spidir

from rasmus.bio import alignlib, fasta, phylo, genomeutil

# scipy libs
#from scipy import matrixmultiply as mm
#from scipy import transpose, array, dot
#from scipy.linalg import expm, logm, eig



options = [
"""\
    Simulates sequence evolution using SPIDIR's model. 
""",

    ["s:", "stree=", "stree", "<species tree>",
        {"single": True}],
    ["S:", "smap=", "smap", "<gene2species map>"],
    ["p:", "params=", "params", "<params>",
        {"single": True}],
    ["t:", "tree=", "tree", "<gene tree topology>",
        {"single": True}],
    ["f:", "famrate=", "famrate", "<family rate>",
        {"single": True,
         "default": None,
         "parser": float}],
    ["n:", "ntrees=", "ntrees", "<number of trees to produce>",
        {"single": True,
         "default": 1,
         "parser": int}], 
    ["", "start=", "start", "<starting number>",
        {"single": True,
         "default": 0,
         "parser": int}],
    ["m:", "midpoints=", "midpoints", "<midpoints file>",
        {"single": True}],
            
    "Output extensions",
    ["T:", "outtree=", "outtree", "<output prefix>",
        {"single": True,
         "default": "./"}],
    
    ["", "outtreeext=", "outtreeext", "<output extension>",
        {"single": True,
         "default": ".tree"}],

#    ["d:", "outdist=", "outdist", "<output distance filename extension>",
#        {"single": True,
#         "default": ".dist"}],
    
]


conf = util.parseOptions(sys.argv, options, quit=True)



def main(conf):
    env.addEnvPaths("DATAPATH")
    
    # read configuration    
    stree = treelib.readTree(env.findFile(conf["stree"]))
    params = Spidir.readParams(conf["params"])
    gene2species = genomeutil.readGene2species(env.findFile(conf["smap"][-1]))
    
    tree = treelib.readTree(conf["tree"])
    
    #if "midpoints" in conf:
    #    midpoints = readMidpoints(conf["midpoints"])
    
    
    # simulate
    util.tic("simulating %d trees" % conf["ntrees"])
    for i in range(conf["start"], conf["ntrees"]):
        util.log("simulating", i)
        
        # create tree
        tree2 = simTree(conf, tree, stree, params, gene2species)
        
        
        # write output
        tree2.write(conf["outtree"] + "/" + str(i) + conf["outtreeext"])
    util.toc()


def simTree(conf, tree, stree, params, gene2species):
    if conf["famrate"] is None:
        baserate = random.gammavariate(params["baserate"][0],
                                       1 / params["baserate"][1])
    else:
        baserate = conf["famrate"]

    tree = simTreeLens(conf, tree, stree, params, gene2species, baserate)
    
    return tree



def simTreeLens(conf, tree, stree, params, gene2species, baserate):
    tree2 = tree.copy()
    
    recon = phylo.reconcile(tree2, stree, gene2species)
    events = phylo.labelEvents(tree2, recon)
    
    # set midpoint randomly
    midpoints = {}
    setMidpointsRandom(tree2.root, events, recon, params, midpoints, wholeTree=True)
    
    #util.printDict(midpoints, keyfunc=lambda x: x.name)
    
    # set branch lengths
    def walk(node):
        if node.parent and \
           recon[node] == recon[node.parent]:
            startpoint = midpoints[node.parent]
        else:
            startpoint = 0.0
        
        # if not root, set branch length
        if node.parent != None:
            param = determineParam(node, recon, events, stree, params, midpoints)
            node.dist = simBranchLen(conf, baserate, param)
        
        node.recurse(walk)
    walk(tree2.root)
    
    
    return tree2



def determineParam(node, recon, events, stree, params, midpoints):
    fracs, snodes = reconBranch(node, recon, events, stree, params, midpoints)
    
    param = [0, 0]
    
    for frac, snode in zip(fracs, snodes):
        param[0] += frac * params[snode][0]
        param[1] += frac * params[snode][1]**2
    
    # convert variance to sdev
    param[1] = math.sqrt(param[1])
    
    return param



def reconBranch(node, recon, events, stree, params, midpoints):
    fracs = []
    snodes = []

    if recon[node] == recon[node.parent]:
        # we begin and end on same branch
        fracs.append(midpoints[node] - midpoints[node.parent])
        snodes.append(recon[node].name)
    else:
        # we begin and end on different branches

        # ending (most recent)
        if recon[node] != stree.root:
            fracs.append(midpoints[node])
            snodes.append(recon[node].name)

        # walk up until starting species branch
        snode = recon[node].parent
        snode2 = recon[node.parent]
        while snode != snode2:
            fracs.append(1) # full branch
            snodes.append(snode.name)
            snode = snode.parent

        # beginning branch
        if midpoints[node.parent] != 1.0 and \
           snode2 != stree.root:
            fracs.append(1.0 - midpoints[node.parent])
            snodes.append(snode2.name)
    
    return fracs, snodes


def setMidpointsRandom(node, events, recon, params, midpoints, wholeTree=True):
    def walk(node):
        # determine this node's midpoint
        if events[node] == "dup" and \
           node.parent != None:
            if recon[node] == recon[node.parent]:
                # if im the same species branch as my parent 
                # then he is my last midpoint
                startpoint = midpoints[node.parent]
            else:
                # im the first on this branch so the last midpoint is zero
                #startpoint = 0.0
                
                startpoint = chooseStartpoint(node, recon, params)
            
            # pick a midpoint uniformly after the last one
            midpoints[node] = startpoint + \
                        (random.random() * (1 - startpoint))
        else:
            # genes or speciations reconcile exactly to the end of the branch
            # gene tree roots also reconcile exactly to the end of the branch
            midpoints[node] = 1.0
    
        # recurse within dup-only subtree
        if events[node] == "dup" or wholeTree:
            node.recurse(walk)
    
    walk(node)


def chooseStartpoint(node, recon, params):
    originalLen = params[recon[node].name][0]
    child = filter(lambda x: recon[x] == recon[node], node.children)[0]
    grandChildren = child.children
    
    minlen = min(map(lambda x: params[recon[x].name][0], grandChildren))
    startpoint = max(0, 1 - (minlen / originalLen))
    
    return startpoint


def simBranchLen(conf, baserate, param):
    blen = 0
    
    while blen <= 0:
        blen = baserate * random.normalvariate(param[0], param[1])
    
    return blen


main(conf)

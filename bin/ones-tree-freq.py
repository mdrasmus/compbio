#!/usr/bin/env python

from rasmus import env
from rasmus import genomeutil
from rasmus import phyloutil
from rasmus import treelib
from rasmus import util


import sys

options = [
    ["s:", "stree=", "stree", "<species tree>",
        {"single": True}],
    ["S:", "smap=", "smap", "<gene to species mapping>",
        {"single": True}]
]


conf = util.parseOptions(sys.argv, options, quit=True)


env.addEnvPaths("DATAPATH")
gene2species = genomeutil.readGene2species(env.findFile(conf["smap"]))
stree = treelib.readTree(env.findFile(conf["stree"]))


correctSplits = phyloutil.findBranchSplits(stree)
splitLookup = util.revdict(correctSplits)

splitCounts = {}
for branch in correctSplits.itervalues():
    splitCounts[branch] = 0

ntrees = 0

for treefile in sys.stdin:
    tree = treelib.readTree(treefile.rstrip())
    ntrees += 1
    
    # map genes to species
    for leaf in tree.leaves():
        tree.rename(leaf.name, gene2species(leaf.name))
        splits = phyloutil.findBranchSplits(tree)
        
        for branch in splits.itervalues():
            if branch in splitCounts:
                splitCounts[branch] += 1


for branch, num in splitCounts.iteritems():
    branch2 = splitLookup[branch]
    
    # get proper node for branch
    node = stree.nodes[branch2[0]]
    node2 = stree.nodes[branch2[1]]
    if node2.parent == node:
        node = node2
    
    node.data["boot"] = num / float(ntrees)


stree.write()

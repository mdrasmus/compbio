#!/usr/bin/env python

from rasmus import treelib, util
from rasmus.bio import genomeutil, phylo

import sys

options = [
  ["t:", "tree=", "tree", "<newick file>",
    {"default": []}],
  ["c:", "cost=", "cost", "dup|loss|duploss",
    {"default": "duploss",
     "single": True}],
#  ["r:", "reroot=", "reroot", "<branch to root tree>",
#    {"single": True}],
] + genomeutil.options


# parse options
conf = util.parseOptions(sys.argv, options, quit=True, resthelp="<trees> ...")
genomeutil.readOptions(conf)

gene2species = conf["gene2species"]

if "stree" in conf:
    stree = conf["stree"]



for treefile in conf["REST"]:
    print "rerooting %s..." % treefile
    tree = treelib.readTree(treefile)
    phylo.reconRoot(tree, stree, gene2species, 
                        rootby=conf["cost"], newCopy=False)
    tree.write(treefile)

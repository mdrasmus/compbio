#!/usr/bin/env python

import sys

from rasmus import cluster
from rasmus import genomeutil
from rasmus import util
from rasmus import treelib


options = [
    ["r:", "relcutoff=", "relcutoff", "<relative cutoff>",
        {"single": True,
         "parser": float,
         "default": .8}],
    ["c:", "signif=", "signif", "<significance cutoff>",
        {"single": True,
         "parser": float,
         "default": 1e-3}],
    ["a", "all", "all", "",
        {"single": True,
         "help": "write all partitionings in the tree"}]
] + genomeutil.options


def main(conf):
    genomeutil.readOptions(conf)
    
    blastfiles = conf[""]
    
    tree = cluster.mergeTree(conf, 
                              conf["stree"], 
                              conf["gene2species"],
                              cluster.makeBlastFileLookup(blastfiles))
    
    if conf["all"]:
        for node in tree.nodes.values():
            if "parts" in dir(node):
                util.writeDelim("node-" + str(node.name) + ".part", node.parts)
    else:
        util.writeDelim(sys.stdout, tree.root.parts)
    
    

conf = util.parseOptions(sys.argv, options, quit=True, resthelp="<blast files>")
main(conf)

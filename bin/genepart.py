#!/usr/bin/env python

import sys

from rasmus import cluster
from rasmus import genomeutil
from rasmus import util
from rasmus import treelib


options = [
    """ note: Blast files should be in NCBI blast -m8 output format
    """,
    
    ["r:", "relcutoff=", "relcutoff", "<relative cutoff>",
        {"single": True,
         "parser": float,
         "default": .8}],
    ["c:", "signif=", "signif", "<significance cutoff>",
        {"single": True,
         "parser": float,
         "default": 1e-3}],
    ["a:", "all=", "all", "<output prefix>",
        {"single": True,
         "default": None,
         "help": "write all partitionings in the tree"}],
    ["m:", "merge=", "merge", "avg|buh",
        {"default": "avg",
         "single": True}],
    ["", "accept=", "accept", "<accept list>"]
] + genomeutil.options


def main(conf):
    genomeutil.readOptions(conf)
    
    blastfiles = conf[""]
    
    accept = set()
    for filename in conf["accept"]:
        for row in util.DelimReader(filename):
            for gene in row:
                accept.add(gene)
    conf["accept"] = accept
    
    if conf["all"] != None:
        conf["output"] = conf["all"]
    
    tree = cluster.mergeTree(conf, 
                             conf["stree"], 
                             conf["gene2species"],
                             cluster.makeBlastFileLookup(blastfiles))
    
    if conf["all"] == None:    
        util.writeDelim(sys.stdout, tree.root.parts)

    

conf = util.parseOptions(sys.argv, options, quit=True, resthelp="<blast files>")
main(conf)


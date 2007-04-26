#!/usr/bin/env python

import sys

from rasmus import genecluster
from rasmus import genomeutil
from rasmus import util
from rasmus import treelib
from rasmus import tablelib


options = [
    """ note: Blast files should be in NCBI blast -m8 output format
    """,
    
    ["r:", "relcutoff=", "relcutoff", "<relative cutoff>",
        {"single": True,
         "parser": float,
         "default": .8}],
    #["c:", "signif=", "signif", "<significance cutoff>",
    #    {"single": True,
    #     "parser": float,
    #     "default": 1e-3}],
    ["b:", "bitspersite=", "bitspersite", "<bits/site cutoff>",
        {"single": True,
         "parser": float,
         "default": 1.0}],
    ["c:", "coverage=", "coverage", "<minimum coverage>",
        {"single": True,
         "parser": float,
         "default": 0.5}],
    ["a:", "all=", "all", "<output prefix>",
        {"single": True,
         "default": None,
         "help": "write all partitionings in the tree"}],
    ["m:", "merge=", "merge", "avg|buh",
        {"default": "avg",
         "single": True}],
    ["l:", "genelens=", "genelens", "<table of gene lengths>",
        {"single": True}],
    ["", "accept=", "accept", "<accept list>"]
] + genomeutil.options


def main(conf):
    genomeutil.readOptions(conf)
    
    blastfiles = conf[""]
    
    accept = set()
    if "accept" in conf:
        for filename in conf["accept"]:
            for row in util.DelimReader(filename):
                for gene in row:
                    accept.add(gene)                    
        conf["accept"] = accept
    
    if conf["all"] != None:
        conf["output"] = conf["all"]
    
    
    genes = tablelib.readTable(conf["genelens"]).lookup("gene")
    
    tree = genecluster.mergeTree(conf, 
                             genes,
                             conf["stree"], 
                             conf["gene2species"],
                             genecluster.makeBlastFileLookup(blastfiles))
    
    if conf["all"] == None:    
        util.writeDelim(sys.stdout, tree.root.parts)

    

conf = util.parseOptions(sys.argv, options, quit=True, resthelp="<blast files>")
main(conf)


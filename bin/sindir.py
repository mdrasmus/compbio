#!/usr/bin/env python
#
# SINDIR - Species INformed DIstance-based Reconciliation
# Matt Rasmussen
# May 2006
#

# python libs
import copy
import math
import os
import random
import StringIO
import sys

# rasmus libs
from rasmus import algorithms
from rasmus import genomeutil
from rasmus import fasta
from rasmus import phylip
from rasmus import phyloutil
from rasmus import sindirlib
from rasmus import stats
from rasmus import util

    

# command line options
options = [
 ["o:", "out=", "out", "<output file>",
    {"help": "Tree output file (Newick format) (default: stdout)", 
     "single": True,
     "default": sys.stdout}],
 
 ["d:", "dist=", "dist", "<distance file>", 
    {"help": "Build a tree from a phylip distance matrix file"}],

 ["l:", "labels=", "labels", "<gene labels file>", 
    {"help": "Supply gene names for the rows and cols of a distance matrix"}],
     
 ["p:", "param=", "param", "<param file>",
    {"help": "When training this is used as an output file",
     "single": True,
     "req": True}],
 ["s:", "stree=", "stree", "<species tree>",
    {"single": True,
     "req": True,
     "help": "species tree (newick format)"}],
 ["S:", "smap=", "smap", "<gene2species map>",
    {"req": True,
     "help": "mapping of gene names to species names"}],

 #["a:", "align=", "align", "<fasta alignment>",
 #   {"help": "Build a tree from a multiple sequence alignment in fasta format"}],
     
 
 "Search options",
 ["R:", "search=", "search", "<method>",
    {"default": ["mcmc"],
     "help": """\
    tree search method.  Options are:
        mcmc   - Markov Chain Monte Carlo (default)
        greedy - Phylip-like search
        none   - no search"""}],
 ["I:", "depth=", "depth", "<NNI depth>",
    {"default": 3,
     "parser": int,
     "single": True}],
 ["i:", "iters=", "iters", "<iterations>", 
    {"default": 2000, 
     "parser": int,
     "single": True}],
 ["T:", "tree=", "tree", "<propose tree>",
    {"default": []}],     

 "Training options",
 ["t:", "traindir=", "traindir", "<training directory>",
    {"help": """\
    Run SINDIR in training mode.  """,
    "single":True}],
 ["e:", "trainext=", "trainext", "<training tree extension>",
    {"single":True}],
 ["z:", "trainstats=", "trainstats", "<training stats prefix>",
    {"single": True,
     "default": ""}],
 
 "Parameter options",
 ["A:", "accuracy=", "accuracy", "<accuracy>",
    {"default": .1,
     "parser": float,
     "single": True}],
 ["D:", "dupprob=", "dupprob", "<probability of duplication>",
    {"default": .01,
     "parser": float,
     "single": True}],
 ["L:", "lossprob=", "lossprob", "<probability of loss>",
    {"default": .01,
     "parser": float,
     "single": True}],
 
 "Miscellaneous options",
 ["V:", "debug=", "debug", "<level>",
    {"single": True,
     "default": 0, 
     "parser": int,
     "help": "set SINDIR debug level"}],

 ["P:", "paths=", "paths", "<files path>",
    {"single": True,
     "help": "colon separated paths used to search for data files",
     "default": "."}]
]

 

def main(argv):   
    # parse options
    conf = util.parseOptions(argv, options, quit=True)
    
    if conf["debug"] > 0:
        util.globalTimer().removeStream(sys.stderr)
        util.globalTimer().addStream(sys.stdout)
        print "SINDIR"
        print "configuration:"
        util.printDict(conf, justify=lambda x: "left")
        print
        print
    
    genomeutil.readOptions(conf)
    conf["specprob"] = 1.0 - conf["dupprob"]

    
    # read input
    gene2species = conf["gene2species"]
    stree = conf["stree"]    
    
    if "traindir" in conf:
        trainTree(conf, stree, gene2species)
    else:
        buildTree(conf, stree, gene2species)




def trainTree(conf, stree, gene2species):
    treefiles = util.listFiles(conf["traindir"], conf["trainext"])

    util.tic("reading trees")
    trees = []
    prog = util.ProgressBar(len(treefiles))
    for treefile in treefiles:
        prog.update()
        trees.append(algorithms.readTree(treefile))
    util.toc()
    
    params = sindirlib.learnModel(trees, stree, gene2species, conf["trainstats"])
    
    sindirlib.writeParams(conf["param"], params)
    
    


def buildTree(conf, stree, gene2species):
    params = sindirlib.readParams(conf["param"])

    if "dist" in conf:
        for i in range(len(conf["dist"])):
            distfile = conf["dist"][i]
            
            labels, distmat = phylip.readDistMatrix(distfile)
        
            # read in different labels if needed
            if "labels" in conf:
                labels = sindirlib.readLabels(conf["labels"][i])
            
            tree, logl = sindirlib.sindir(conf, distmat, labels, stree, 
                                          gene2species, params)
            tree.write(conf["out"])



#
# run as script
#
if __name__ == "__main__":
    try:
        main(sys.argv)
    except util.OptionError, e:
        print >>sys.stdout, "%s: %s" % (os.path.basename(sys.argv[0]), e)
        sys.exit(1)




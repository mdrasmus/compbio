#!/usr/bin/env python

from rasmus import env
from rasmus import fasta
from rasmus import genomeutil
from rasmus import phylip
from rasmus import phyloutil
from rasmus import util
from rasmus import treelib
import sys


options = [
 ["S:", "smap=", "smap", "<gene2speces map>"],
 ["s:", "stree=", "stree", "<species tree>",
    {"single": True}],
 ["F:", "fastaext=", "fastaext", "<fasta extension>",
    {"single": True}],
 ["T:", "treeext=", "treeext", "<tree extension>",
    {"single": True}]
]

conf = util.parseOptions(sys.argv, options, quit=True)


def main(conf):
    env.addEnvPaths("DATAPATH")

    files = conf["REST"]
    gene2species = genomeutil.readGene2species(* map(env.findFile, conf["smap"]))
    stree = treelib.readTree(env.findFile(conf["stree"]))
    
    util.tic("species tree -> gene tree")
    
    i = 0
    for f in files:
        i += 1
        util.logger("%d of %d (%s)" % (i, len(alnfiles), f))
        
        seqs = fasta.readFasta(f)
        tree = phyloutil.stree2gtree(stree, seqs.keys(), gene2species)
        tree.writeNewick(f.replace(conf["fastaext"], conf["treeext"]))

    util.toc()

main(conf)

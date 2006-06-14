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
 ["s:", "stree=", "stree", "<species tree>"],
 ["A:", "alignext=", "alignext", "<alignment extension>"],
 ["T:", "treeext=", "treeext", "<tree extension>"]
]

param = util.parseOptions(sys.argv, options, quit=True)


def main(param):
    env.addEnvPaths("DATAPATH")

    alnfiles = param[""]
    gene2species = genomeutil.readGene2species(* map(env.findFile, param["smap"]))
    stree = treelib.readTree(env.findFile(param["stree"][-1]))
    
    util.tic("make ML trees from stree")
    
    i = 0
    for f in alnfiles:
        i += 1
        util.logger("%d of %d (%s)" % (i, len(alnfiles), f))

        aln = fasta.readFasta(f)
        tree = phyloutil.stree2gtree(stree, aln.keys(), gene2species)

        logl, tree = phylip.promlTreelk(aln, tree, verbose=False)
        tree = phyloutil.reconRoot(tree, stree, gene2species)
        
        tree.writeNewick(f.replace(param["alignext"][-1], param["treeext"][-1]))

    util.toc()

main(param)

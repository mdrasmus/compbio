#!/usr/bin/env python

from rasmus.common import *
from rasmus.bio import phylo
import Spidir
import Spidir.Likelihood

tree = readTree("../test/0.nt.tree")
stree = readTree("../test/flies.stree")
gene2species = genomeutil.readGene2species("../test/flies.smap")
params = Spidir.readParams("../test/flies.nt.param")

drawTree(tree)

print sum(x.dist for x in tree)
print Spidir.estGeneRate(tree, stree, params, gene2species)
print Spidir.estGeneRate(tree, stree, params, gene2species)
print Spidir.estGeneRate(tree, stree, params, gene2species)
print Spidir.Likelihood.getBaserate(tree, stree, params, gene2species=gene2species)

conf = {"python_only": True, 
        "famprob": True}
print Spidir.treeLogLikelihood(conf, tree, stree, gene2species, params)

conf = {}
#generate = Spidir.estGeneRate(tree, stree, params, gene2species)

for generate in frange(1.5, 2.3, .05):
    print generate, Spidir.treeLogLikelihood(conf, tree, stree, gene2species, params, 
                                             baserate=generate)

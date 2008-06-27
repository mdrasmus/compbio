#!/usr/bin/env python

from rasmus.common import *
from rasmus.bio import phylo
import Spidir

import spidir

tree = readTree("test/data/0.nt.tree")
stree = readTree("test/data/flies.stree")
gene2species = genomeutil.readGene2species("test/data/flies.smap")
params = Spidir.readParams("test/data/flies.nt.param")
aln = readFasta("test/data/1.nt.align")
bgfreqs = [.258,.267,.266,.209]
tsvratio = 1.59

drawTree(tree)

print sum(x.dist for x in tree)
print Spidir.estGeneRate(tree, stree, params, gene2species)
print Spidir.estGeneRate(tree, stree, params, gene2species)
print Spidir.Likelihood.getBaserate(tree, stree, params,
                                    gene2species=gene2species)

#conf = {"python_only": True, 
#        "famprob": True}
#print Spidir.treeLogLikelihood(conf, tree, stree, gene2species, params)

util.tic("sample")
generates = spidir.sample_gene_rate(tree, stree, gene2species, params,
                                    aln, bgfreqs, tsvratio, 20000)
util.toc()

print Spidir.estGeneRate(tree, stree, params, gene2species)
print mean(generates), (percentile(generates, .050),
                        percentile(generates, .950))

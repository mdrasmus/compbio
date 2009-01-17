#!/usr/bin/env python

from rasmus.common import *
import sys
import pyspidir
reload(pyspidir)

from rasmus import treelib
from rasmus.bio import genomeutil
import Spidir

def gene2speciesArray(tree, stree, gene2species):
    ptree, nodes, nodelookup = Spidir.makePtree(tree)
    sptree, snodes, snodelookup = Spidir.makePtree(stree)

    gene2speciesarray = []
    for node in nodes:
        if node.isLeaf():
            gene2speciesarray.append(snodelookup[
                                     stree.nodes[gene2species(node.name)]])
        else:
            gene2speciesarray.append(-1)
    return gene2speciesarray
    

tree = treelib.readTree("../data/0.nt.tree")
tree = treelib.parseNewick("((dmoj_sim,(dvir_sim,dgri_sim)),(dwil_sim,(dpse_sim,(dana_sim,(dmel_sim,(dyak_sim,dere_sim))))));")
tree = treelib.parseNewick("((dwil_sim,(dvir_sim,dgri_sim)),(dmoj_sim,(dpse_sim,(dana_sim,(dmel_sim,(dyak_sim,dere_sim))))));")

stree = treelib.readTree("../data/flies.stree")
gene2species = genomeutil.readGene2species("../data/flies.smap")
params = Spidir.readParams("../data/flies.nt.param")


ptree, nodes, nodelookup = Spidir.makePtree(tree)
pstree, snodes, snodelookup = Spidir.makePtree(stree)
g2s = gene2speciesArray(tree, stree, gene2species)
mu = [float(params[snode.name][0]) for snode in snodes]
sigma = [float(params[snode.name][1]) for snode in snodes]
alpha = float(params['baserate'][0])
beta = float(params['baserate'][1])


dists = pyspidir.genbranches(ptree, pstree, g2s, mu, sigma, alpha, beta)

for node, d in zip(nodes, dists):
    node.dist = d

drawTree(tree)
print tree.root.children[0].dist



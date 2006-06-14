#!/usr/bin/env python

from rasmus import algorithms, stats, util, phylip, fasta
from rasmus import ensembl, phyloutil, genomeutil, env
import sys, os
import math, StringIO, copy, random

from compbio.tools import pp

from rasmus import sindirlib


options = [
    ["p:", "part=", "part", "<part file>"],
    ["n:", "nbrs=", "nbrs", "<neighbor file>"],
    ["s:", "stree=", "stree", "<species tree>"],
    ["S:", "smap=", "smap", "<gene2species map>"],
    ["d:", "outdir=", "outdir", "<output directory>"],
    ["o:", "outpart=", "outpart", "<output partition>"]
]


param = util.parseOptions(sys.argv, options, quit=True)


def main(param):
    env.addEnvPaths("DATAPATH")
    
    util.tic("reading input")
    gene2species = genomeutil.readGene2species(env.findFile(param["smap"][-1]))
    stree = algorithms.readTree(env.findFile(param["stree"][-1]))
    genomes = tuple(util.sort(stree.leaveNames()))
    parts = util.readDelim(param["part"][-1])
    orthonbrs = util.map2(int, util.readDelim(param["nbrs"][-1]))
    util.toc()
    
    
    # make trees
    util.tic("make trees")
    trees = makeSpecificityTest(parts, orthonbrs, stree, gene2species)
    util.toc()
    
    util.tic("write output")
    # make output partition
    parts2 = [tree.leaveNames() for tree in trees]
    util.writeDelim(param["outpart"][-1], parts2)
    
    # write all trees
    for i, tree in enumerate(trees):
        tree.writeNewick(os.path.join(param["outdir"][-1], "%d.correct.tree" % i))
    
    util.toc()
        


def makeSpecificityTest(ones, orthonbrs, stree, gene2species):
    genomes = stree.leaveNames()
    trees = []
    
    for groups in orthonbrs:
        # chose which group to take each species
        take = {}
        for genome in genomes:
            take[genome] = random.randint(0, 1)

        # create list of genes to remove
        dels = filter(lambda x: take[gene2species(x)] != 0, ones[groups[0]])
        dels += filter(lambda x: take[gene2species(x)] != 1, ones[groups[1]])
        
        # construct new tree
        tree1 = makeGeneTree(ones[groups[0]], stree, gene2species)
        tree2 = makeGeneTree(ones[groups[1]], stree, gene2species)
        tree = algorithms.Tree()
        tree.makeRoot(0)
        tree.addTree(tree.root, tree1)
        tree.addTree(tree.root, tree2)
        tree = removeLeaves(tree, dels)
        
        # save tree
        trees.append(tree)
        
    return trees


def makeGeneTree(part, stree, gene2species):
    tree = stree.copy()
    lookup = {}

    # make reverse lookup species2gene
    for gene in part:
        lookup[gene2species(gene)] = gene

    for leaf in stree.leaveNames():
        tree.rename(leaf, lookup[leaf])
    
    return tree


def removeLeaves(tree, leaves, depth=0):
    for leaf in leaves:
        tree.remove(tree.nodes[leaf])
    exposed = filter(lambda x: type(x) == int, tree.leaveNames())
    if len(exposed) > 0:
        tree = removeLeaves(tree, exposed, 1)

    if depth == 0:
        algorithms.removeSingleChildren(tree)
    
    return tree
  



main(param)

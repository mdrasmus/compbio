#
# Python module for SPIDIR
#


#
# requires rasmus library
#

from math import *
from ctypes import *

from rasmus import treelib, util
from rasmus.bio import fasta

import pyspidir
spidir = cdll.LoadLibrary("lib/libspidir.so")


def export(lib, funcname, return_type, arg_types, scope=globals()):
    """Exports a C function"""

    scope[funcname] = lib.__getattr__(funcname)
    scope[funcname].restype = return_type
    scope[funcname].argtypes = arg_types


export(spidir, "numHistories", c_int, [c_int])
export(spidir, "birthDeathCount", c_float, [c_int, c_float, c_float, c_float])
export(spidir, "makeTree", c_void_p, [c_int, POINTER(c_int)])



#=============================================================================


def make_ptree(tree):
    """Make parent tree array from tree"""

    nodes = []
    nodelookup = {}
    ptree = []
    
    def walk(node):
        for child in node.children:
            walk(child)
        nodes.append(node)
    walk(tree.root)
    
    def leafsort(a, b):
        if a.isLeaf():
            if b.isLeaf():
                return 0
            else:
                return -1
        else:
            if b.isLeaf():
                return 1
            else:
                return 0
    
    # bring leaves to front
    nodes.sort(cmp=leafsort)
    nodelookup = util.list2lookup(nodes)
    
    for node in nodes:
        if node == tree.root:
            ptree.append(-1)
        else:
            ptree.append(nodelookup[node.parent])
    
    assert nodes[-1] == tree.root
    
    return ptree, nodes, nodelookup


def ptree2ctree(ptree):
    """Makes a c++ Tree from a parent array"""
    
    pint = c_int * len(ptree)
    tree = makeTree(len(ptree), pint(* ptree))
    return tree


def tree2ctree(tree):
    """Make a c++ Tree from a treelib.Tree datastructure"""

    ptree, nodes, nodelookup = make_ptree(tree)
    return ptree2ctree(ptree)


def make_gene2species_array(stree, nodes, snodelookup, gene2species):
    gene2speciesarray = []
    for node in nodes:
        if node.isLeaf():
            gene2speciesarray.append(snodelookup[
                                     stree.nodes[gene2species(node.name)]])
        else:
            gene2speciesarray.append(-1)
    return gene2speciesarray


def parsimony(aln, tree):    
    ptree, nodes, nodelookup = make_ptree(tree)
    leaves = [x.name for x in nodes if isinstance(x.name, str)]
    seqs = util.mget(aln, leaves)
    
    dists = pyspidir.parsimony(ptree, seqs)
    
    for i in xrange(len(dists)):
        nodes[i].dist = dists[i]


def sample_gene_rate(tree, stree, gene2species, params,
                     aln,
                     bgfreq, tsvratio,
                     nsamples=100):

    ptree, nodes, nodelookup = make_ptree(tree)
    pstree, snodes, snodelookup = make_ptree(stree)
    smap = make_gene2species_array(stree, nodes, snodelookup, gene2species)

    mu = [float(params[snode.name][0]) for snode in snodes]
    sigma = [float(params[snode.name][1]) for snode in snodes]
    alpha, beta = params['baserate']

    seqs = [aln[node.name] for node in nodes if isinstance(node.name, str)]

    generates = []
    def callback(generate):
        generates.append(generate)
    
    pyspidir.sample_gene_rate(nsamples,
                              ptree, 
                              pstree,
                              smap,
                              mu,
                              sigma,
                              alpha, beta,
                              bgfreq,
                              tsvratio,
                              seqs,
                              callback)

    return generates



def est_gene_rate(tree, stree, gene2species, params,
                  aln, bgfreq, tsvratio,
                  nsamples=1000):

    generates = sample_gene_rate(tree, stree, gene2species, params,
                                 aln,
                                 bgfreq, tsvratio,
                                 nsamples)
    generates.sort()
    low = generates[int(.05 * len(generates))]
    high = generates[int(.95 * len(generates))]
    
    return sum(generates) / float(nsamples), (low, high)


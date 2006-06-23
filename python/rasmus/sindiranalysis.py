import util
import treelib
import phyloutil
import sindirlib
import stats

from common import *


def showrates(lens, low=0, high=.5, step=.001):
    p, params, resid = stats.plotdistribFit(stats.normalPdf, [1,1], 
                                            lens, low, high, step)
    p.set(xmin=low, xmax=high, main="mean=%f, sdev=%f" % tuple(params))
    
    return p
    
    
def showparams(g):
    pc(map(util.flatten, g.params.items()))


def showstree(g):
    treelib.drawTreeNames(g.stree, minlen=8)

def showcorr(data, keys):
    heatmap(corrmatrix(data), 
            colormap=ColorMap([[-1,red],
                               [-.4, yellow],
                               [0,black],
                               [.4,green],
                               [1, blue]]), 
    rlabels=keys,
    clabels=keys, 
    xmargin=100, 
    ymargin=100)



def makeGeneTree(part, stree, gene2species):
    tree = stree.copy()
    lookup = {}

    # make reverse lookup species2gene
    for gene in part:
        lookup[gene2species(gene)] = gene

    for leaf in stree.leaveNames():
        tree.rename(leaf, lookup[leaf])

    return tree



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
        dels += filter(lambda x: take[gene2species(x)] != 1,
                       ones[groups[1]])

        # construct new tree
        tree1 = makeGeneTree(ones[groups[0]], stree, gene2species)
        tree2 = makeGeneTree(ones[groups[1]], stree, gene2species)
        tree = Tree()
        tree.makeRoot(0)
        tree.addTree(tree.root, tree1)
        tree.addTree(tree.root, tree2)
        tree = removeLeaves(tree, dels)

        # save tree
        trees.append(tree)

    return trees


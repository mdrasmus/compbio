#!/usr/bin/python -i
#
# ultimate tree/alignment/distmat viewer
#


import sys

from summon.core import *
from summon import multiwindow, sumtree, matrix
import summon

from rasmus import treelib, util
from rasmus.bio import alignlib, fasta, phylip
from rasmus.vis.genomebrowser import *


options = [
    ["a:", "align=", "align", "<fasta alignment>",
        {"single": True}],
    ["t:", "tree=", "tree", "<newick file>",
        {"single": True}],
    ["d:", "distmat=", "distmat", "<phylip distance matrix>",
        {"single": True}],
    ["M:", "maxdist=", "maxdist", "<maximum distance>",
        {"single": True,
         "parser": float}]
    ]


conf = util.parseOptions(sys.argv, options)

      
def colorAlign(aln):
    if guessAlign(aln) == "pep":
        return pep_colors
    else:
        return dna_colors
        

windows = []
coords = []


# read tree
if "tree" in conf:
    tree = treelib.readTree(conf["tree"])
    vistree = sumtree.SumTree(tree, name=conf["tree"])
                              #xscale=conf["usedist"]
    vistree.show()
    leaves = tree.leafNames()
    leaves.reverse()
    
    windows.append(vistree.win)
    coords.append(max(node.x for node in tree.nodes.itervalues()) + 1.5)


# read alignment
if "align" in conf:
    view = Region("", "", "", 1, 1)
    aln = fasta.readFasta(conf["align"])
    
    original_order = aln.keys()
    
    if "tree" in conf:
        aln = aln.get(leaves)
    
    view.end = max(view.end, aln.alignlen())
    height = len(aln)
    colors = colorAlign(aln)


# read distance matrix
if "distmat" in conf:
    util.tic("reading matrix")
    label, distmat = phylip.readDistMatrix(conf["distmat"])
    util.toc()

    if "align" in conf:
        label = original_order
    
    if "tree" in conf:
        lookup = util.list2lookup(label)
        rperm = util.mget(lookup, leaves)
        cperm = util.mget(lookup, leaves)
    else:
        rperm = None
        cperm = None

    if "maxdist" in conf:
        colormap = util.rainbowColorMap(low=0, high=conf["maxdist"])
    else:
        colormap = util.ColorMap([[-1e-10, util.black],
                                  [0, util.blue],
                                  [.5, util.green],
                                  [1.0, util.yellow],
                                  [5.0, util.red],
                                  [10, util.white]])
    
    visdist = matrix.DenseMatrixViewer(distmat, colormap=colormap, 
                                rlabels=label, clabels=label, 
                                cutoff=-util.INF,
                                rperm=rperm, cperm=cperm)
    visdist.show()
    
    windows.append(visdist.win)
    coords.append(1.5)


# show alignment
if "align" in conf:
    visalign = GenomeStackBrowser(view=view)
    visalign.addTrack(RulerTrack(bottom=-height))
    visalign.addTrack(AlignTrack(aln, colorBases=colors))
    visalign.show()

    windows.append(visalign.win)
    coords.append(0)


# tie all windows by their y-coordinate
if len(windows) > 1:
    multiwindow.tie_windows(windows, tiey=True, piny=True, coordsy=coords)
    e = multiwindow.WindowEnsembl(windows, stacky=True, sameh=True)

#!/usr/bin/python -i
#
# ultimate tree/alignment/distmat viewer
#


import sys

from summon.core import *
from summon import multiwindow, sumtree, matrix
import summon

from rasmus import treelib, util
from rasmus.bio import alignlib, blast, fasta, phylip
from rasmus.vis.genomebrowser import *


options = [
    ["a:", "align=", "align", "<fasta alignment>",
        {"single": True}],
    ["t:", "tree=", "tree", "<newick file>",
        {"single": True}],
    ["d:", "distmat=", "distmat", "<phylip distance matrix>"],
    ["b:", "blast=", "blast", "<blast hit -m8 format>"],
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


def readBlastMatrix(blastfile, order=None):
    mat = matrix.Matrix()

    qnames = []
    snames = []
    qseen = set()
    sseen = set()
    
    queries = []
    subjects = []
    vals = []
    
    for hit in blast.BlastReader(blastfile):
        q = hit[0]
        s = hit[1]
        v = blast.bitscore(hit)
        
        # record label order
        if q not in qseen:
            qnames.append(q)
            qseen.add(q)
        if s not in sseen:
            snames.append(s)
            sseen.add(s)
        
        # store data
        queries.append(q)
        subjects.append(s)
        vals.append(v)

    # determine label order
    if order:
        labels = order
    else:
        labels = qnames + filter(lambda x: x not in sseen, snames)
    # set label names    
    mat.rowlabels = labels
    mat.collabels = labels
    
    
    # populate matrix        
    lookup = util.list2lookup(labels)
    
    rows, cols, vals2 = mat.rows, mat.cols, mat.vals
    for i in xrange(len(vals)):
        r = lookup[queries[i]]
        c = lookup[subjects[i]]
        rows.append(r)
        cols.append(c)
        vals2.append(vals[i])
        mat[r][c] = vals[i]
    mat.setup(len(labels), len(labels), len(vals))

    return mat




        
#=============================================================================
#

windows = []
coords = []


# read tree
if "tree" in conf:
    tree = treelib.readTree(conf["tree"])
    vistree = sumtree.SumTree(tree, name=conf["tree"])
                              #xscale=conf["usedist"]
    vistree.show()
    vistree.win.set_size(340, 500)
    vistree.win.set_position(0, 0)
    
    leaves = tree.leafNames()
    leaves.reverse()
    
    windows.append(vistree.win)
    coords.append(max(node.x for node in tree.nodes.itervalues()))


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
if "distmat" in conf or "blast" in conf:
    mats = []
    matnames = []
    
    currentMatrix = 0

    # determine row/col labels from alignment if it exists
    if "align" in conf:
        label = original_order
    else:
        label = None
    
    # reorder according to any given tree
    if "tree" in conf:
        lookup = util.list2lookup(label)
        rperm = util.mget(lookup, leaves)
        cperm = util.mget(lookup, leaves)
    else:
        rperm = []
        cperm = []
    
    
    # setup color map
    if "maxdist" in conf:
        colormap = util.rainbowColorMap(low=0, high=conf["maxdist"])
    else:
        colormap = util.ColorMap([[-1e-10, (1, .7, .8)],
                                  [0, util.black],
                                  [1e-10, util.blue],
                                  [.5, util.green],
                                  [1.0, util.yellow],
                                  [5.0, util.red],
                                  [10, util.white]])    
    
    
    # read in multiple distance matrices
    if "distmat" in conf:
        for matfile in conf["distmat"]:
            util.tic("reading matrix '%s'" % matfile)
            label, mat = phylip.readDistMatrix(matfile)
            util.toc()

            # convert distmatrix to summon Matrix
            mat2 = matrix.Matrix()
            mat2.from2DList(mat, cutoff=-util.INF)
            mat2.colormap = colormap
            mat2.rowlabels = label
            mat2.collabels = label
            mat2.rperm = rperm
            mat2.cperm = cperm
            mat2.setup()

            mats.append(mat2)
            matnames.append(matfile)
    
    # read in multiple blast hits
    if "blast" in conf:
        for blastfile in conf["blast"]:
            util.tic("reading blast file '%s'" % blastfile)
            mat = readBlastMatrix(blastfile, order=label)
            util.toc()
            
            mat.colormap = colormap
            mat.rperm = rperm
            mat.cperm = cperm
            mat.setup()
            mats.append(mat)
            matnames.append(blastfile)


    
    # allow easy switching between matrices
    def nextMatrix():
        global currentMatrix
        currentMatrix = (currentMatrix + 1) % len(mats)
        visdist.setMatrix(mats[currentMatrix])
        visdist.win.set_name(matnames[currentMatrix])
        visdist.redraw()
    
    def prevMatrix():
        global currentMatrix
        currentMatrix = (currentMatrix - 1) % len(mats)
        visdist.setMatrix(mats[currentMatrix])
        visdist.win.set_name(matnames[currentMatrix])
        visdist.redraw()
    
    
    # create matrix vis
    visdist = matrix.MatrixViewer(mats[0], conf={"style": "quads"})
    visdist.show()
    visdist.win.set_name(matnames[0])
    visdist.win.set_bgcolor(1, 1, 1)    
    visdist.win.set_size(300, 500)
    visdist.win.set_position(0, 0)
    
    windows.append(visdist.win)
    coords.append(0)
    
    
    visdist.win.set_binding(input_key("n"), nextMatrix)
    visdist.win.set_binding(input_key("p"), prevMatrix)    
    

# show alignment
if "align" in conf:
    visalign = GenomeStackBrowser(view=view)
    visalign.addTrack(RulerTrack(bottom=-height))
    visalign.addTrack(AlignTrack(aln, colorBases=colors))
    visalign.show()
    visalign.win.set_name(conf["align"])
    visalign.win.set_size(580, 500)
    visalign.win.set_position(0, 0)    

    windows.append(visalign.win)
    coords.append(-1.5)


# tie all windows by their y-coordinate
if len(windows) > 1:
    multiwindow.tie_windows(windows, tiey=True, piny=True, coordsy=coords)
    e = multiwindow.WindowEnsembl(windows, stacky=True, sameh=True)

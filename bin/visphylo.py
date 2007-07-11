#!/usr/bin/python -i
#
# Thu May 17 11:42:05 EDT 2007 
# ultimate tree/alignment/distmat viewer
#



import sys

from summon.core import *
from summon import multiwindow, sumtree, matrix
import summon

from rasmus import treelib, util
from rasmus.bio import alignlib, blast, fasta, phylip
from rasmus.vis import distmatrixvis
from rasmus.vis.genomebrowser import *


options = [
    "data files",
    ["f:", "fasta=", "fasta", "<fasta sequences>",
        {"single": True}],
    ["a:", "align=", "align", "<fasta alignment>",
        {"single": True}],
    ["t:", "tree=", "tree", "<newick file>",
        {"single": True}],
    ["d:", "distmat=", "distmat", "<phylip distance matrix>"],
    ["b:", "blast=", "blast", "<blast hit -m8 format>"],
    
    "data file extensions",
    ["F:", "fastaext=", "fastaext", "<fasta sequences extension>",
        {"default": []}],
    ["A:", "alignext=", "alignext", "<alignment extension>",
        {"default": []}],
    ["T:", "treeext=", "treeext", "<tree file extension>",
        {"default": []}],
    ["D:", "distmatext=", "distmatext", "<phylip distance matrix extension>",
        {"default": []}],
    ["B:", "blastext=", "blastext", "<blast hit -m8 format extension>",
        {"default": []}],
    
    
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

# determine data files
fastafiles = []
alignfiles = []
distmatfiles = []
blastfiles = []
treefiles = []

# add explicit filenames
if "fasta" in conf:
    fastafiles.append(conf["fasta"])
if "align" in conf:
    alignfiles.append(conf["align"])
if "distmat" in conf:
    distmatfiles.extend(conf["distmat"])
if "blast" in conf:
    blastfiles.extend(conf["blast"])
if "tree" in conf:
    treefiles.append(conf["tree"])
        

# add any implicit files with matching extension
if len(conf["REST"]) > 0:
    basename = conf["REST"][0]

    for ext in conf["fastaext"]:
        fastafiles.append(basename + ext)
    
    for ext in conf["alignext"]:
        alignfiles.append(basename + ext)
    
    for ext in conf["distmatext"]:
        distmatfiles.append(basename + ext)
    
    for ext in conf["blastext"]:
        blastfiles.append(basename + ext)
    
    for ext in conf["treeext"]:
        treefiles.append(basename + ext)




# read sequences
if len(fastafiles) > 0:
    seqs = fasta.readFasta(fastafiles[0])
else:
    seqs = None
    

# read tree
if len(treefiles) > 0:
    tree = treelib.readTree(treefiles[0])
    vistree = sumtree.SumTree(tree, name=treefiles[0])
                              #xscale=conf["usedist"]
    vistree.show()
    vistree.win.set_size(340, 500)
    vistree.win.set_position(0, 0)
    
    leaves = tree.leafNames()
    
    windows.append(vistree.win)
    coords.append(max(node.y for node in tree.nodes.itervalues()))
else:
    leaves = None

# read alignment
if len(alignfiles) > 0:
    view = Region("", "", "", 1, 1)
    aln = fasta.readFasta(alignfiles[0])
    
    original_order = aln.keys()
    
    if len(treefiles) > 0:
        aln = aln.get(leaves)
    
    view.end = max(view.end, aln.alignlen())
    height = len(aln)
    colors = colorAlign(aln)


# read distance matrix
if len(distmatfiles) > 0 or len(blastfiles) > 0:
    mats = []
    matnames = []
    
    currentMatrix = 0

    # determine row/col labels from alignment if it exists
    if len(alignfiles) > 0:
        label = original_order
    else:
        label = None
    
    # reorder according to any given tree
    if leaves != None:
        lookup = util.list2lookup(label)
        rperm = util.mget(lookup, leaves)
        cperm = util.mget(lookup, leaves)
    else:
        rperm = []
        cperm = []
    
    
    # setup color map
    if "maxdist" in conf:
        low = 0
        high = conf["maxdist"]
        colormap = util.ColorMap([[-1e-10, (1, .7, .8)],
                                  [0, util.black],
                                  [1e-10, util.blue],
                                  [.3 * high, util.green],
                                  [.7 * high, util.yellow],
                                  [     high, util.red]])
    else:
        colormap = util.ColorMap([[-1e-10, (1, .7, .8)],
                                  [0, util.black],
                                  [1e-10, util.blue],
                                  [.5, util.green],
                                  [1.0, util.yellow],
                                  [5.0, util.red],
                                  [10, util.white]])    
    
    
    # read in multiple distance matrices
    for matfile in distmatfiles:
        util.tic("reading matrix '%s'" % matfile)
        label2, mat = phylip.readDistMatrix(matfile)
        util.toc()
        
        if label == None:
            label = label2

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
    for blastfile in blastfiles:
        if label == None and len(fastafiles) > 0:
            label = seqs.keys()
            rperm = []
            cperm = []
    
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
    visdist = distmatrixvis.DistMatrixViewer(mats[0], seqs=seqs)
    visdist.show()
    visdist.win.set_name(matnames[0])
    visdist.win.set_bgcolor(1, 1, 1)    
    #visdist.win.set_size(300, 500)
    #visdist.win.set_position(0, 0)
    
    windows.append(visdist.win)
    coords.append(0)
    
    
    visdist.win.set_binding(input_key("n"), nextMatrix)
    visdist.win.set_binding(input_key("p"), prevMatrix)    
    

# show alignment
if len(alignfiles) > 0:
    visalign = GenomeStackBrowser(view=view)
    visalign.addTrack(RulerTrack(bottom=-height))
    visalign.addTrack(AlignTrack(aln, colorBases=colors))
    visalign.show()
    visalign.win.set_name(alignfiles[0])
    visalign.win.set_size(580, 500)
    visalign.win.set_position(0, 0)    

    windows.append(visalign.win)
    coords.append(-1.5)


# tie all windows by their y-coordinate
if len(windows) > 1:
    e = multiwindow.WindowEnsemble(windows, stacky=True, sameh=True,
                                   tiey=True, coordsy=coords)


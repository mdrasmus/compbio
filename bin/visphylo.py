#!/usr/bin/env python-i
#
# Thu May 17 11:42:05 EDT 2007 
# ultimate tree/alignment/distmat viewer
#


# python libs
import sys

# summon libs
from summon.core import *
from summon import matrix, sumtree
import summon

# rasmus libs
from rasmus import treelib, util
from rasmus.bio import alignlib, blast, fasta, phylip, genomeutil
from rasmus.vis import phylovis, treevis



options = [
"""\
Visualizes phylogenetic trees, distance matrices, and sequence alignments
simultaneously.

"""

    "data files",
    ["f:", "fasta=", "fasta", "<fasta sequences>",
        {"single": True}],
    ["a:", "align=", "align", "<fasta alignment>"],
    ["t:", "tree=", "tree", "<newick file>"],
    ["d:", "distmat=", "distmat", "<phylip distance matrix>"],
    ["b:", "blast=", "blast", "<blast hit -m8 format>"],
    
    ["s:", "stree=", "stree", "<species tree newick file>",
        {"single": True}],
    ["S:", "smap=", "smap", "<gene to species mapping>",
        {"single": True}],
    
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
    
    "optional",
    ["M:", "maxdist=", "maxdist", "<maximum distance>",
        {"single": True,
         "parser": float,
         "help": "set max distance in matrix colormap"}],
    ["", "treecolormap=", "treecolormap", "<tree color map file>",
        {"single": True}]
    ]


conf = util.parseOptions(sys.argv, options)

      

def readBlastMatrix(blastfile, order=None):
    """Reads a BLAST file (-m8 format) into a sparse SUMMON Matrix"""
    
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
    
    # blast hits have their own color map
    mat.colormap = util.ColorMap([[0, (0, 0, 0)],
                                  [100, (0, 0, 1)],
                                  [300, (0, 1, 0)],
                                  [500, (1, 1, 0)],
                                  [1000, (1, 0, 0)]])

    return mat




        
#=============================================================================
# determine input filenames

treefiles = []
distmatfiles = []
blastfiles = []
alignfiles = []
fastafiles = []

treefiles.extend(conf.get("tree", []))
distmatfiles.extend(conf.get("distmat", []))
blastfiles.extend(conf.get("blast", []))
alignfiles.extend(conf.get("align", []))

# add explicit filenames
if "fasta" in conf:
    fastafiles.append(conf["fasta"])


# add any implicit files with matching extension
if len(conf["REST"]) > 0:
    basename = conf["REST"][0]
    
    for ext in conf["treeext"]:
        treefiles.append(basename + ext)
    
    for ext in conf["distmatext"]:
        distmatfiles.append(basename + ext)
    
    for ext in conf["blastext"]:
        blastfiles.append(basename + ext)
    
    for ext in conf["alignext"]:
        alignfiles.append(basename + ext)
        
    for ext in conf["fastaext"]:
        fastafiles.append(basename + ext)
    
    
#=============================================================================
# read tree
trees = []
for filename in treefiles:
    trees.append(treelib.readTree(filename))



# read species information

if "smap" in conf:
    gene2species = genomeutil.readGene2species(conf["smap"])
else:
    gene2species = None

if "stree" in conf:
    stree = treelib.readTree(conf["stree"])
else:
    stree = None
    
if "treecolormap" in conf:
    treeColormap = treevis.read_tree_color_map(conf["treecolormap"])
else:
    treeColormap = None

#=============================================================================
# read distance matrix
distmats = []
distmatNames = []
distmatLabels = []
    
        
# read in multiple distance matrices
for matfile in distmatfiles:
    label, mat = phylip.readDistMatrix(matfile)
    distmats.append(mat)
    distmatLabels.append(label)
    distmatNames.append(matfile)    

# read in multiple blast hits
for blastfile in blastfiles:
    mat = readBlastMatrix(blastfile)
    distmats.append(mat)
    distmatLabels.append(label)
    distmatNames.append(blastfile)        

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
    colormap = None

# read sequences
if len(fastafiles) > 0:
    seqs = fasta.readFasta(fastafiles[0])
else:
    seqs = None


#=============================================================================
# read alignment
aligns = []
for filename in alignfiles:
    aligns.append(fasta.readFasta(filename))



#=============================================================================
# create PhyloViewer

vis = phylovis.PhyloViewer(trees, distmats, aligns,
                  
                           # tree config
                           stree=stree,
                           gene2species=gene2species,
                           treeColormap=treeColormap,

                           # distmat config
                           distlabels=distmatLabels, 
                           matrixColormap=colormap,
                           seqs=seqs,

                           # filenames
                           treeNames=treefiles,
                           distmatNames=distmatNames,                  
                           alignNames=alignfiles)
vis.show()

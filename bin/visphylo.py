#!/usr/bin/env python-i
#
# Thu May 17 11:42:05 EDT 2007 
# ultimate tree/alignment/distmat viewer
#


# python libs
import sys

# summon libs
from summon.core import *
from summon import multiwindow, sumtree, matrix
import summon

# rasmus libs
from rasmus import treelib, util
from rasmus.bio import alignlib, blast, fasta, phylip, genomeutil
from rasmus.vis import distmatrixvis, alignvis, treevis
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
    
    
    ["M:", "maxdist=", "maxdist", "<maximum distance>",
        {"single": True,
         "parser": float}]
    ]


conf = util.parseOptions(sys.argv, options)

      

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
    
    # blast hits have their own color map
    mat.colormap = util.ColorMap([[0, (0, 0, 0)],
                                  [100, (0, 0, 1)],
                                  [300, (0, 1, 0)],
                                  [500, (1, 1, 0)],
                                  [1000, (1, 0, 0)]])

    return mat




        
#=============================================================================
#




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


#=============================================================================
# read species information

if "smap" in conf:
    gene2species = genomeutil.readGene2species(conf["smap"])
else:
    gene2species = None

if "stree" in conf:
    stree = treelib.readTree(conf["stree"])
else:
    stree = None


#=============================================================================


class PhyloViewer (object):
    
    def __init__(self, trees=[], distmats=[], aligns=[],
                       stree=None,
                       gene2species=None,
                       
                       distlabels=None, 
                       matrixColormap=None,
                       distlabelsFromAlign=True,
                       seqs=None,
                       
                       treeWinSize=(350, 500),
                       distmatWinSize=None,
                       alignWinSize=(580, 500),
                       
                       treeNames=None, distmatNames=None, alignNames=None):
        """
        Initialize PhyloViewer
        
        NOTE: distlabelsFromAlign (bool) determines whether a distance matrix 
        should ignore its own labels (which are probablu dummy labels) and 
        assume its rows and columns are ordered the same as the given alignment.
        """
        
        # data
        self.trees = trees
        self.distmats = distmats
        self.aligns = aligns
        
        # optional data
        self.stree = stree
        self.gene2species = gene2species
        self.distlabels = distlabels
        self.distlabelsFromAlign = distlabelsFromAlign
        self.seqs = seqs
        
        # names
        self.treeNames = treeNames
        self.distmatNames = distmatNames
        self.alignNames = alignNames
        
        # default name
        if self.treeNames == None and len(self.trees) > 0:
            self.treeNames = ["tree %d" for i in xrange(1, len(self.trees)+1)]
        if self.distmatNames == None and len(self.distmats) > 0:
            self.distmatNames = ["matrix %d" for i in xrange(1, len(self.distmats)+1)]
        if self.alignNames == None and len(self.aligns) > 0:
            self.alignNames = ["alignment %d" for i in xrange(1, len(self.aligns)+1)]
            
        
        self.order = None
        self.alignOrder = None
        
        # currently displayed
        self.currentTree = None
        self.currentMatrix = None
        self.currentAlign = None
        
        # window sizes
        self.treeWinSize = treeWinSize
        self.distmatWinSize = distmatWinSize
        self.alignWinSize = alignWinSize
        
        # window ensemble properties
        self.windows = []
        self.coords = []
        
        # colormap
        if matrixColormap != None:
            self.matrixColormap = matrixColormap
        else:
            # setup default
            self.matrixColormap = util.ColorMap([[-1e-10, (1, .7, .8)],
                                                 [0, util.black],
                                                 [1e-10, util.blue],
                                                 [.5, util.green],
                                                 [1.0, util.yellow],
                                                 [5.0, util.red],
                                                 [10, util.white]])
        
        # initialize data
        self.initTrees()
        self.initAligns()
        self.initDistmats()        
    
    
    def initTrees(self):
        """Initialize trees"""
        
        if len(self.trees) > 0:
            self.currentTree = self.trees[0]
            self.vistree = treevis.TreeViewer(self.currentTree, name=self.treeNames[0],
                                         xscale=100.0,
                                         stree=self.stree,
                                         gene2species=self.gene2species)            
            self.order = self.currentTree.leafNames()
        else:
            self.vistree = None
            self.order = None
    
    
    def initDistmats(self):
        """Initialize distance matrices
        
           Initialization should by done after trees and alignments
        """
    
        if len(self.distmats) > 0:
            self.matrices = []
            
            # setup matrices            
            for i, distmat in enumerate(self.distmats):
                # convert distmatrix to summon Matrix
                if isinstance(distmat, matrix.Matrix):
                    mat = distmat
                else:
                    mat = matrix.Matrix()
                    mat.from2DList(distmat)            

                # set default colormap
                if mat.colormap == None:
                    mat.colormap = self.matrixColormap
                
                # determine labels
                if self.distlabelsFromAlign and self.alignOrder != None:
                    # determine row/col labels from alignment if it exists
                    mat.rowlabels = self.alignOrder
                    mat.collabels = self.alignOrder
                
                elif self.distlabels != None:
                    mat.rowlabels = self.distlabels[i]
                    mat.collabels = self.distlabels[i]
                    
                else:
                    raise Exception("no labels given for matrix")
                
                # reorder according to any given tree
                if self.order != None:
                    lookup = util.list2lookup(mat.rowlabels)
                    mat.rperm = util.mget(lookup, self.order)
                    mat.cperm = util.mget(lookup, self.order)
                
                mat.setup()

                self.matrices.append(mat)
        
            if self.seqs == None:
                seqs = self.currentAlign
            else:
                seqs = self.seqs
            
            # create matrix vis
            self.currentMatrix = self.matrices[0]
            self.visdist = distmatrixvis.DistMatrixViewer(self.currentMatrix, 
                                                          seqs=seqs, 
                                                          bgcolor=(1,1,1))
        else:
            self.visdist = None
    
    
    def initAligns(self):
        """Initialize sequence alignments"""
    
        if len(self.aligns) > 0:
            self.currentAlign = self.aligns[0]
            
            self.alignOrder = self.currentAlign.keys()
            if self.order != None:
                self.currentAlign.names = self.order
            else:
                self.order = self.currentAlign.keys()

            self.visalign = alignvis.AlignViewer(self.currentAlign, 
                                                 size=self.alignWinSize,
                                                 title=self.alignNames[0])
        else:
            self.visalign = None
    
    
    
    def show(self):
        """Make visualization window visible"""
        
        # tree visualization
        if self.vistree != None:
            self.vistree.show()
            self.vistree.win.set_size(*self.treeWinSize)
            
            self.windows.append(self.vistree.win)
            self.coords.append(max(node.y for node in self.currentTree.nodes.itervalues()))
        
        # distance matrix visualization
        if self.visdist != None:            
            self.visdist.show()
            self.visdist.win.set_name(self.distmatNames[0])
            
            self.windows.append(self.visdist.win)
            self.coords.append(0)
            
            # add additional menu options
            self.visdist.bar.addItem(hud.MenuItem("next matrix (n)", self.nextMatrix))
            self.visdist.bar.addItem(hud.MenuItem("next matrix (p)", self.prevMatrix))
            
            # add additional key binding
            self.visdist.win.set_binding(input_key("n"), self.nextMatrix)
            self.visdist.win.set_binding(input_key("p"), self.prevMatrix)


        # alignment visualization
        if self.visalign != None:
            self.visalign.show()        
            self.windows.append(self.visalign.vis.win)
            self.coords.append(-1.5)        
        
        
        # tie all windows by their y-coordinate
        if len(self.windows) > 1:
            self.ensembl = multiwindow.WindowEnsemble(self.windows, stacky=True, 
                                                      sameh=True, tiey=True, 
                                                      coordsy=self.coords)


    #=============================================
    # allow easy switching between matrices
    def nextMatrix(self):
        matindex = self.matrices.index(self.currentMatrix)
        matindex = (matindex + 1) % len(self.matrices)
        self.currentMatrix = self.matrices[matindex]
        self.visdist.setMatrix(self.currentMatrix)
        self.visdist.win.set_name(self.distmatNames[matindex])
        self.visdist.redraw()

    def prevMatrix(self):
        matindex = self.matrices.index(self.currentMatrix)
        matindex = (matindex - 1) % len(self.matrices)
        self.currentMatrix = self.matrices[matindex]
        self.visdist.setMatrix(self.currentMatrix)
        self.visdist.win.set_name(self.distmatNames[matindex])
        self.visdist.redraw()



    
#=============================================================================
# read tree
trees = []
for filename in treefiles:
    trees.append(treelib.readTree(filename))


#=============================================================================
# read alignment
aligns = []
for filename in alignfiles:
    aligns.append(fasta.readFasta(filename))



#=============================================================================
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


#=============================================================================
# create PhyloViewer

vis = PhyloViewer(trees, distmats, aligns,
                  
                  # tree config
                  stree=stree,
                  gene2species=gene2species,
                  
                  # distmat config
                  distlabels=distmatLabels, 
                  matrixColormap=colormap,
                  seqs=seqs,
                  
                  # filenames
                  treeNames=treefiles,
                  distmatNames=distmatNames,                  
                  alignNames=alignfiles)
vis.show()

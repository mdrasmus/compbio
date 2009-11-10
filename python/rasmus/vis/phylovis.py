"""

    Phylogeny Visualization
    
    Visualizes trees, distance matrices, and sequence alignments simultaneously.
    
"""

from rasmus import util
from rasmus.vis import distmatrixvis, alignvis, treevis

import summon
from summon import matrix, hud, multiwindow
from summon.core import *




class PhyloViewer (object):
    """
    Phylogeny Visualization
    
    Visualizes trees, distance matrices, and sequence alignments simultaneously.
    """
    
    def __init__(self, trees=[], distmats=[], aligns=[],
                       stree=None,
                       gene2species=None,
                       treeColormap=lambda x: (0, 0, 0),
                       
                       distlabels=None, 
                       matrixColormap=None,
                       distlabelsFromAlign=True,
                       seqs=None,
                       
                       treeWinSize=(350, 500),
                       distmatWinSize=None,
                       alignWinSize=(580, 500),
                       
                       treeNames=None, distmatNames=None, alignNames=None):
        """
        arguments:
        
        trees    -- a list of trees
        distmats -- a list of distance matrices (2D list or summon.matrix.Matrix)
        aligns   -- a list of sequence alignments 
        
        * Tree optional configuration
            If both stree and gene2species are given, duplications and losses
            will be plotted.
        stree        -- species tree
        gene2species -- gene name to species name mapping function
                       
        * Distance matrix configuration
        distlabels      -- a list of row/col labels (one for each matrix)
        matrixColormap  -- default colormap for matrices
        distlabelsFromAlign -- distlabelsFromAlign (bool) determines whether a
                               distance matrix  should ignore its own labels
                               (which are probably dummy labels) and  assume its
                               rows and columns are ordered the same as the
                               given alignment.
        seqs            -- a dict of sequences for visualizing BLAST matrices
                       
        * Window sizes
        treeWinSize    -- tree visualization window size. default: (350, 500),
        distmatWinSize -- matrix visualization window size. 
        alignWinSize   -- alignment visualization window size. default: (580, 500)
        
        * Window titles (optional)
        treeNames    -- a list of tree names (one for each tree) to be used in 
                        tree window title. default: None
        distmatNames -- a list fo matric names (one for each matrix) to be used
                        in matrix window title. default: None
        alignNames   -- a list of alignment names (one for each alignment) to 
                        be used in alignment window title. default: None
        """
        
        # data
        self.trees = trees
        self.distmats = distmats
        self.aligns = aligns
        
        # optional data
        self.stree = stree
        self.gene2species = gene2species
        self.treeColormap = treeColormap
        self.distlabels = distlabels
        self.distlabelsFromAlign = distlabelsFromAlign
        self.seqs = seqs
        self.matrices = []
        
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
            self.vistree = PhyloTreeViewer(self.currentTree, name=self.treeNames[0],
                                         phyloViewer=self,
                                         xscale=100.0,
                                         stree=self.stree,
                                         gene2species=self.gene2species,
                                         winsize=self.treeWinSize,
                                         colormap=self.treeColormap)
            self.vistree.setTree(self.currentTree)
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
            
            winpos = None
            self.visalign = alignvis.AlignViewer(self.currentAlign, 
                                                 winsize=self.alignWinSize,
                                                 winpos=winpos,
                                                 title=self.alignNames[0])                                                 
        else:
            self.visalign = None
    
    
    
    def show(self):
        """Make visualization window visible"""
        
        # tree visualization
        if self.vistree != None:
            self.vistree.show()
            #self.vistree.win.set_size(*self.treeWinSize)
            
            self.windows.append(self.vistree.win)
            self.coords.append(max(node.y for node in self.currentTree.nodes.itervalues()))
        
            # add additional menu options
            self.vistree.bar.addItem(hud.MenuItem("next tree (n)", self.nextTree))
            self.vistree.bar.addItem(hud.MenuItem("prev tree (p)", self.prevTree))
            
            # add additional key binding
            self.vistree.win.set_binding(input_key("n"), self.nextTree)
            self.vistree.win.set_binding(input_key("p"), self.prevTree)
        
        
        # distance matrix visualization
        if self.visdist != None:            
            self.visdist.show()
            self.visdist.win.set_name(self.distmatNames[0])
            
            self.windows.append(self.visdist.win)
            self.coords.append(0)
            
            # add additional menu options
            self.visdist.bar.addItem(hud.MenuItem("next matrix (n)", self.nextMatrix))
            self.visdist.bar.addItem(hud.MenuItem("prev matrix (p)", self.prevMatrix))
            
            # add additional key binding
            self.visdist.win.set_binding(input_key("n"), self.nextMatrix)
            self.visdist.win.set_binding(input_key("p"), self.prevMatrix)


        # alignment visualization
        if self.visalign != None:
            self.visalign.show()        
            self.windows.append(self.visalign.vis.win)
            self.coords.append(-1.5)        

            # add additional key binding
            self.visalign.win.set_binding(input_key("n"), self.nextAlign)
            self.visalign.win.set_binding(input_key("p"), self.prevAlign)
                
        
        # tie all windows by their y-coordinate
        if len(self.windows) > 1:
            self.ensembl = multiwindow.WindowEnsemble(self.windows, stacky=True, 
                                                      sameh=True, tiey=True, 
                                                      coordsy=self.coords)


    def onReorderLeaves(self):
        leaves = self.currentTree.leafNames()
        
        # reorder matrix
        for mat in self.matrices:
            lookup = util.list2lookup(mat.rowlabels)
            mat.rperm = util.mget(lookup, leaves)
            mat.cperm = util.mget(lookup, leaves)
            mat.setup()
        if self.visdist:
            self.visdist.redraw()
        
        
        # reorder alignment
        for aln in self.aligns:
            aln.names = leaves
        if self.visalign:
            self.visalign.show()


    #=============================================
    # allow easy switching between trees
    def nextTree(self):
        treeindex = self.trees.index(self.currentTree)
        treeindex = (treeindex + 1) % len(self.trees)
        self.currentTree = self.trees[treeindex]        
        self.vistree.setTree(self.currentTree)
        self.vistree.win.set_name(self.treeNames[treeindex])
        self.vistree.show()
        self.onReorderLeaves()

    def prevTree(self):
        treeindex = self.trees.index(self.currentTree)
        treeindex = (treeindex - 1) % len(self.trees)
        self.currentTree = self.trees[treeindex]
        self.vistree.setTree(self.currentTree)
        self.vistree.win.set_name(self.treeNames[treeindex])
        self.vistree.show()
        self.onReorderLeaves()


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

    #=============================================
    # allow easy switching between alignments
    def nextAlign(self):
        alnindex = self.aligns.index(self.currentAlign)
        alnindex = (alnindex + 1) % len(self.aligns)
        self.currentAlign = self.aligns[alnindex]
        self.visalign.setAlign(self.currentAlign)
        self.visalign.win.set_name(self.alignNames[alnindex])
        self.visalign.show()

    def prevAlign(self):
        alnindex = self.aligns.index(self.currentAlign)
        alnindex = (alnindex - 1) % len(self.aligns)
        self.currentAlign = self.aligns[alnindex]
        self.visalign.setAlign(self.currentAlign)
        self.visalign.win.set_name(self.alignNames[alnindex])
        self.visalign.show()

        

class PhyloTreeViewer (treevis.TreeViewer):
    def __init__(self, *args, **options):
        
        if "phyloViewer" in options:
            phyloViewer = options["phyloViewer"]
            del options["phyloViewer"]
        treevis.TreeViewer.__init__(self, *args, **options)
        
        self.phyloViewer = phyloViewer

    def onReorderLeaves(self):
        self.phyloViewer.onReorderLeaves()


    def setTree(self, tree):

        if max(node.dist for node in tree) == 0.0:
            self.xscale = 0.0
        else:
            self.xscale = 100.0
            
        treevis.TreeViewer.setTree(self, tree)

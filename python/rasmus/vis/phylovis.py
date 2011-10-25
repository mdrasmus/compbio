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
                       tree_colormap=lambda x: (0, 0, 0),
                       
                       dist_labels=None, 
                       matrix_colormap=None,
                       dist_labels_from_align=True,
                       seqs=None,
                       
                       tree_win_size=(350, 500),
                       distmat_win_size=None,
                       align_win_size=(580, 500),
                       
                       tree_names=None, distmat_names=None, align_names=None):
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
        dist_labels      -- a list of row/col labels (one for each matrix)
        matrix_colormap  -- default colormap for matrices
        dist_labels_from_align -- distlabelsFromAlign (bool) determines
                               whether a distance matrix  should ignore its
                               own labels (which are probably dummy labels)
                               and  assume its rows and columns are ordered
                               the same as the given alignment.
        seqs            -- a dict of sequences for visualizing BLAST matrices
                       
        * Window sizes
        tree_win_size    -- tree visualization window size. default: (350, 500),
        distmat_win_size -- matrix visualization window size. 
        align_win_size   -- alignment visualization window size. default: (580, 500)
        
        * Window titles (optional)
        tree_names    -- a list of tree names (one for each tree) to be used in 
                         tree window title. default: None
        distmat_names -- a list fo matric names (one for each matrix) to be used
                         in matrix window title. default: None
        align_names   -- a list of alignment names (one for each alignment) to 
                         be used in alignment window title. default: None
        """
        
        # data
        self.trees = trees
        self.distmats = distmats
        self.aligns = aligns
        
        # optional data
        self.stree = stree
        self.gene2species = gene2species
        self.tree_colormap = tree_colormap
        self.distlabels = dist_labels
        self.dist_labels_from_align = dist_labels_from_align
        self.seqs = seqs
        self.matrices = []
        
        # names
        self.tree_names = tree_names
        self.distmat_names = distmat_names
        self.align_names = align_names
        
        # default name
        if self.tree_names == None and len(self.trees) > 0:
            self.tree_names = ["tree %d" for i in xrange(1, len(self.trees)+1)]
        if self.distmat_names == None and len(self.distmats) > 0:
            self.distmat_names = ["matrix %d" for i in xrange(1, len(self.distmats)+1)]
        if self.align_names == None and len(self.aligns) > 0:
            self.align_names = ["alignment %d" for i in xrange(1, len(self.aligns)+1)]
            
        
        self.order = None
        self.align_order = None
        
        # currently displayed
        self.current_tree = None
        self.current_matrix = None
        self.current_align = None
        
        # window sizes
        self.tree_win_size = tree_win_size
        self.distmat_win_size = distmat_win_size
        self.align_win_size = align_win_size
        
        # window ensemble properties
        self.windows = []
        self.coords = []
        
        # colormap
        if matrix_colormap != None:
            self.matrix_colormap = matrix_colormap
        else:
            # setup default
            self.matrix_colormap = util.ColorMap([[-1e-10, (1, .7, .8)],
                                                  [0, util.black],
                                                  [1e-10, util.blue],
                                                  [.5, util.green],
                                                  [1.0, util.yellow],
                                                  [5.0, util.red],
                                                  [10, util.white]])
        
        # initialize data
        self.init_trees()
        self.init_aligns()
        self.init_distmats()        
    
    
    def init_trees(self):
        """Initialize trees"""
        
        if len(self.trees) > 0:
            self.current_tree = self.trees[0]
            self.vistree = PhyloTreeViewer(self.current_tree,
                                           name=self.tree_names[0],
                                           phylo_viewer=self,
                                           xscale=100.0,
                                           stree=self.stree,
                                           gene2species=self.gene2species,
                                           winsize=self.tree_win_size,
                                           colormap=self.tree_colormap)
            self.vistree.set_tree(self.current_tree)
            self.order = self.current_tree.leaf_names()
        else:
            self.vistree = None
            self.order = None
    
    
    def init_distmats(self):
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
                    mat.colormap = self.matrix_colormap
                
                # determine labels
                if self.dist_labels_from_align and self.align_order != None:
                    # determine row/col labels from alignment if it exists
                    mat.rowlabels = self.align_order
                    mat.collabels = self.align_order
                
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
                seqs = self.current_align
            else:
                seqs = self.seqs
            
            # create matrix vis
            self.current_matrix = self.matrices[0]
            self.visdist = distmatrixvis.DistMatrixViewer(self.current_matrix, 
                                                          seqs=seqs, 
                                                          bgcolor=(1,1,1))
        else:
            self.visdist = None
    
    
    def init_aligns(self):
        """Initialize sequence alignments"""
    
        if len(self.aligns) > 0:
            self.current_align = self.aligns[0]
            
            self.align_order = self.current_align.keys()
            if self.order != None:
                self.current_align.names = self.order
            else:
                self.order = self.current_align.keys()
            
            winpos = None
            self.visalign = alignvis.AlignViewer(self.current_align, 
                                                 winsize=self.align_win_size,
                                                 winpos=winpos,
                                                 title=self.align_names[0])                                                 
        else:
            self.visalign = None
    
    
    
    def show(self):
        """Make visualization window visible"""
        
        # tree visualization
        if self.vistree != None:
            self.vistree.show()
            #self.vistree.win.set_size(*self.tree_win_size)
            
            self.windows.append(self.vistree.win)
            self.coords.append(max(node.y for node in self.current_tree.nodes.itervalues()))
        
            # add additional menu options
            self.vistree.bar.add_item(hud.MenuItem("next tree (n)", self.next_tree))
            self.vistree.bar.add_item(hud.MenuItem("prev tree (p)", self.prev_tree))
            
            # add additional key binding
            self.vistree.win.set_binding(input_key("n"), self.next_tree)
            self.vistree.win.set_binding(input_key("p"), self.prev_tree)
        
        
        # distance matrix visualization
        if self.visdist != None:            
            self.visdist.show()
            self.visdist.win.set_name(self.distmat_names[0])
            
            self.windows.append(self.visdist.win)
            self.coords.append(0)
            
            # add additional menu options
            self.visdist.bar.add_item(hud.MenuItem("next matrix (n)", self.next_matrix))
            self.visdist.bar.add_item(hud.MenuItem("prev matrix (p)", self.prev_matrix))
            
            # add additional key binding
            self.visdist.win.set_binding(input_key("n"), self.next_matrix)
            self.visdist.win.set_binding(input_key("p"), self.prev_matrix)


        # alignment visualization
        if self.visalign != None:
            self.visalign.show()        
            self.windows.append(self.visalign.vis.win)
            self.coords.append(-1.5)        

            # add additional key binding
            self.visalign.win.set_binding(input_key("n"), self.next_align)
            self.visalign.win.set_binding(input_key("p"), self.prev_align)
                
        
        # tie all windows by their y-coordinate
        if len(self.windows) > 1:
            self.ensembl = multiwindow.WindowEnsemble(
                self.windows, stacky=True, 
                sameh=True, tiey=True, 
                coordsy=self.coords)


    def on_reorder_leaves(self):
        leaves = self.current_tree.leaf_names()
        
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
    def next_tree(self):
        treeindex = self.trees.index(self.current_tree)
        treeindex = (treeindex + 1) % len(self.trees)
        self.current_tree = self.trees[treeindex]        
        self.vistree.set_tree(self.current_tree)
        self.vistree.win.set_name(self.tree_names[treeindex])
        self.vistree.show()
        self.on_reorder_leaves()

    def prev_tree(self):
        treeindex = self.trees.index(self.current_tree)
        treeindex = (treeindex - 1) % len(self.trees)
        self.current_tree = self.trees[treeindex]
        self.vistree.set_tree(self.current_tree)
        self.vistree.win.set_name(self.tree_names[treeindex])
        self.vistree.show()
        self.on_reorder_leaves()


    #=============================================
    # allow easy switching between matrices
    def next_matrix(self):
        matindex = self.matrices.index(self.current_matrix)
        matindex = (matindex + 1) % len(self.matrices)
        self.current_matrix = self.matrices[matindex]
        self.visdist.set_matrix(self.current_matrix)
        self.visdist.win.set_name(self.distmat_names[matindex])
        self.visdist.redraw()

    def prev_matrix(self):
        matindex = self.matrices.index(self.current_matrix)
        matindex = (matindex - 1) % len(self.matrices)
        self.current_matrix = self.matrices[matindex]
        self.visdist.set_matrix(self.current_matrix)
        self.visdist.win.set_name(self.distmat_names[matindex])
        self.visdist.redraw()

    #=============================================
    # allow easy switching between alignments
    def next_align(self):
        alnindex = self.aligns.index(self.current_align)
        alnindex = (alnindex + 1) % len(self.aligns)
        self.current_align = self.aligns[alnindex]
        self.visalign.set_align(self.current_align)
        self.visalign.win.set_name(self.align_names[alnindex])
        self.visalign.show()

    def prev_align(self):
        alnindex = self.aligns.index(self.current_align)
        alnindex = (alnindex - 1) % len(self.aligns)
        self.current_align = self.aligns[alnindex]
        self.visalign.set_align(self.current_align)
        self.visalign.win.set_name(self.align_names[alnindex])
        self.visalign.show()

        

class PhyloTreeViewer (treevis.TreeViewer):
    def __init__(self, *args, **options):
        
        if "phylo_viewer" in options:
            phylo_viewer = options["phylo_viewer"]
            del options["phylo_viewer"]
        treevis.TreeViewer.__init__(self, *args, **options)
        
        self.phylo_viewer = phylo_viewer

    def on_reorder_leaves(self):
        self.phylo_viewer.on_reorder_leaves()


    def set_tree(self, tree):

        if max(node.dist for node in tree) == 0.0:
            self.xscale = 0.0
        else:
            self.xscale = 100.0
            
        treevis.TreeViewer.set_tree(self, tree)

"""
    
    Visualizations for Distance Matrices of molecular sequences

"""

# rasmus libs
from rasmus import util
from rasmus.vis import genomebrowser

from compbio import muscle

# summon libs
from summon.core import *
import summon
from summon import matrix
from summon import hud




class DistMatrixViewer (matrix.MatrixViewer):
    def __init__(self, distmat, seqs=None, 
                 colormap=None, 
                 rlabels=None, clabels=None, cutoff=-util.INF,
                 rperm=None, cperm=None, rpart=None, cpart=None,
                 style="quads",
                 **options):
        
        # setup matrix
        if isinstance(distmat, list):
            # create sparse matrix from dense matrix
            mat = Matrix()
            mat.from_dmat(data, cutoff=cutoff)
        else:
            # already sparse
            mat = distmat
        
        if rlabels != None: mat.rowlabels = rlabels
        if clabels != None: mat.collabels = clabels
        if rperm != None:   mat.rperm = rperm
        if cperm != None:   mat.cperm = cperm
        if rpart != None:   mat.rpart = rpart
        if cpart != None:   mat.cpart = cpart
        mat.setup()
        
        if colormap != None:
            mat.colormap = colormap
        
        matrix.MatrixViewer.__init__(self, mat, on_click=self.on_click,
                                     style=style, drawzeros=True,
                                     **options)
        
        
        # distmatrix specific initialization
        self.seqs = seqs
        self.selgenes = set()
        self.aln = None


    def show(self):
        matrix.MatrixViewer.show(self)
        
        # build sidebar menu
        self.bar = hud.SideBar(self.win, width=150)
        self.bar.add_item(hud.MenuItem("align gene (a)", self.show_align))
        self.bar.add_item(hud.MenuItem("clear genes (d)", self.clear_selection))
        self.bar.add_item(hud.MenuItem("show genes (s)", self.show_selection))
        self.bar.add_item(hud.MenuItem("toggle labels (l)", self.toggle_label_windows))
        
        # register key bindings
        self.win.set_binding(input_key("a"), self.show_align)
        self.win.set_binding(input_key("d"), self.clear_selection)        
        self.win.set_binding(input_key("s"), self.show_selection)
    
    
    def on_click(self, row, col, i, j, val):
        self.selgenes.add(row)
        self.selgenes.add(col)

        print "%s %s (%d, %d) = %f" % (row, col, i, j, val)

    def clear_selection(self):
        self.selgenes.clear()
        print "clear selection"


    def show_selection(self):
        # sort genes by order in matrix
        genes = list(self.selgenes)
        lookup = util.list2lookup(self.mat.rowlabels)
        genes.sort(key=lambda x: lookup[x])
        
        print
        print "selected genes:"
        for i, gene in enumerate(genes):
            print "%3d %s" % (i+1, gene)
    
    
    def show_align(self):
        if self.seqs == None:
            print "cannot build alignment: no sequences are loaded"
            return
        
        if len(self.selgenes) == 0:
            print "cannot build alignment: no sequences selected"
            return
        
        genes = list(self.selgenes)
        lookup = util.list2lookup(self.mat.rowlabels)
           
        self.aln = muscle.muscle(self.seqs.get(genes))
        self.aln.names.sort(key=lambda x: lookup[x])
        
        self.visaln = genomebrowser.show_align(self.aln)
        self.visaln.win.set_name("alignment")

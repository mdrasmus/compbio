"""
    
    Visualizations for Distance Matrices of molecular sequences

"""

# rasmus libs
from rasmus import util
from rasmus.bio import muscle
from rasmus.vis import genomebrowser as gb


# summon libs
from summon.core import *
import summon
from summon import matrix
from summon import hud




class AlignViewer (object):
    def __init__(self, aln, title="alignviewer",
                 winsize=[580, 500], winpos=None, 
                 showColorBases=True, colorBases=True, showBases=True,
                 **options):
        
        self.aln = aln
        self.title = title
        self.winsize = winsize
        self.winpos = winpos
        
        # config options
        self.showColorBases = showColorBases
        self.colorBases = colorBases
        self.showBases = showBases
        
        view = gb.Region("", "", "", 1, aln.alignlen())
        self.vis = gb.GenomeStackBrowser(view=view, winsize=winsize, 
                                         winpos=winpos,
                                         **options)
        self.vis.addTrack(gb.RulerTrack(bottom=-len(aln)))
        self.alntrack = None
        self.setAlign(aln)
        self.bar = None
        

    def setAlign(self, aln):
        if self.alntrack != None:
            self.vis.removeTrack(self.alntrack)
        
        self.aln = aln        
        self.alntrack = gb.AlignTrack(aln, colorBases=self.colorBases,
                                      showColorBases=self.showColorBases,
                                      showBases=self.showBases)
        self.vis.addTrack(self.alntrack)
        
        
    
    def show(self):
        newwin = (self.vis.win == None)
        
        self.vis.show()
        self.win = self.vis.win # for convenience
        
        if newwin:
            self.win.set_name(self.title)
            self.win.set_size(* self.winsize)
        
        
        if self.bar == None:
            # build sidebar menu
            self.bar = hud.SideBar(self.win, width=150)
            self.bar.addItem(hud.MenuItem("toggle color (c)", self.toggleColorBases))
            self.bar.addItem(hud.MenuItem("toggle leftwin (l)", self.toggleLeftWindow))
            self.bar.addItem(hud.MenuItem("always color", 
                lambda: self.setMinColorBases(0)))
        
        # register key bindings
        self.win.set_binding(input_key("c"), self.toggleColorBases)
        self.win.set_binding(input_key("l"), self.toggleLeftWindow)
    
    def enableColorBases(self, enabled=True):
        self.alntrack.showColorBases = enabled
        self.vis.update()
    
    
    def toggleColorBases(self):
        self.enableColorBases(not self.alntrack.showColorBases)
    
    def toggleLeftWindow(self):
        self.vis.enableSideWindows(left=not self.vis.showLeftWindow)

    def setMinColorBases(self, size):
        #self.alntrack.min_base_width = size
        pass



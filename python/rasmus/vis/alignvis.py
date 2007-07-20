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
    def __init__(self, aln, colorBases=True, title="alignviewer",
                 size=[580, 500], showColorBases=True, showBases=True,
                 **options):
        
        self.aln = aln
        self.title = title
        self.size = size
        
        self.vis = gb.GenomeStackBrowser(**options)
        self.vis.addTrack(gb.RulerTrack(bottom=-len(aln)))
        self.vis.addTrack(gb.AlignTrack(aln, colorBases=colorBases,
                                        showColorBases=showColorBases,
                                        showBases=showBases))
    
    def show(self):
        self.vis.show()
        self.vis.win.set_name(self.title)
        self.vis.win.set_size(* self.size)
        
        #self.set_binding(input_key("c"), self.toggleColorBases)
        
    
    def enableColorBases(self, enabled=True):
        self.showColorBases = self.vis.showColorBases = enabled
        self.vis.update()
    
    
    def toggleColorBases(self):
        self.enableColorBases(not self.showColorBases)
    

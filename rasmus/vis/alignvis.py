"""
    
    Visualizations for alignments of molecular sequences

"""

# rasmus libs
from rasmus import util
from rasmus.vis import genomebrowser as gb

from compbio import muscle

# summon libs
from summon.core import *
import summon
from summon import matrix
from summon import hud




class AlignViewer (object):
    def __init__(self, aln, title="alignviewer",
                 winsize=[580, 500], winpos=None, 
                 show_color_bases=True, color_bases=True, show_bases=True,
                 **options):
        
        self.aln = aln
        self.title = title
        self.winsize = winsize
        self.winpos = winpos
        
        # config options
        self.show_color_bases = show_color_bases
        self.color_bases = color_bases
        self.show_bases = show_bases
        
        view = gb.Region("", "", "", 1, aln.alignlen())
        self.vis = gb.GenomeStackBrowser(view=view, winsize=winsize, 
                                         winpos=winpos,
                                         **options)
        self.vis.add_track(gb.RulerTrack(bottom=-len(aln)))
        self.alntrack = None
        self.set_align(aln)
        self.bar = None
        

    def set_align(self, aln):
        if self.alntrack != None:
            self.vis.remove_track(self.alntrack)
        
        self.aln = aln        
        self.alntrack = gb.AlignTrack(aln, color_bases=self.color_bases,
                                      show_color_bases=self.show_color_bases,
                                      show_bases=self.show_bases)
        self.vis.add_track(self.alntrack)
        
        
    
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
            self.bar.add_item(hud.MenuItem("toggle color (c)",
                                           self.toggle_color_bases))
            self.bar.add_item(hud.MenuItem("toggle leftwin (l)",
                                           self.toggle_left_window))
            self.bar.add_item(hud.MenuItem("always color", 
                lambda: self.set_min_color_bases(0)))
        
        # register key bindings
        self.win.set_binding(input_key("c"), self.toggle_color_bases)
        self.win.set_binding(input_key("l"), self.toggle_left_window)
    
    def enable_color_bases(self, enabled=True):
        self.alntrack.show_color_bases = enabled
        self.vis.update()
    
    
    def toggle_color_bases(self):
        self.enable_color_bases(not self.alntrack.show_color_bases)
    
    def toggle_left_window(self):
        self.vis.enableSideWindows(left=not self.vis.showLeftWindow)

    def set_min_color_bases(self, size):
        #self.alntrack.min_base_width = size
        pass



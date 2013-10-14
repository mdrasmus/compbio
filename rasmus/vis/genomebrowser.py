
import math
import copy

from summon.core import *
import summon
from summon import multiwindow
from summon import select
from summon import hud
from summon.multiscale import Multiscale
from summon import VisObject

from rasmus.vis import visual
from rasmus import util, stats

from compbio import regionlib
from compbio.regionlib import Region
from compbio import gff, fasta, alignlib, seqlib


# TODO: refactor the concept of height


class Browser (VisObject):
    """Base class for containing and managing multiple graphical elements"""
    
    def __init__(self, show_left_window=False, show_top_window=False):
        VisObject.__init__(self)
        
        self.leftwin = None     # left margin window
        self.topwin = None      # top margin window
        self.show_left_window = show_left_window
        self.show_top_window = show_top_window
        self._tracks = []
        
    
    def add_track(self, track, index=None):
        track.set_browser(self)
        if index == None:
            self._tracks.append(track)
        else:
            self._tracks.insert(index, track)
    
    def remove_track(self, track):
        track.set_browser(None)
        self._tracks.remove(track)
    
    def clear_tracks(self):
        for track in self._tracks:
            track.set_browser(None)
        self._tracks = []
    
    def update(self):
    
        # propogate updates to tracks
        for track in self._tracks:
            track.update()
    
    def show(self):    
        pass
    
    
    def get_left_window(self):
        return self.leftwin
    
    def get_top_window(self):
        return self.topwin

    def enable_side_windows(self, left=None, top=None):
        if left != None:
            self.show_left_window = left
        if top != None:
            self.show_top_window = top
        
        self.show()
    


class Track:
    def __init__(self, pos=[0.0, 0.0], size=[0.0, 1.0], 
                       pos_offset=[0.0, 0.0],
                       view=None):
        self.height = size[1]
        self.size = size
        self.pos = pos[:]
        self.pos_offset = pos_offset[:]
        self.view = copy.copy(view)
        self.browser = None
    
    
    def set_pos(self, x, y):
        self.pos = [x, y]
    
    def set_size(self, w, h):
        self.size = [w, h]
    
    def get_size(self):
        return self.size
        
    def set_view(self, species, chrom, start, end):
        self.view = Region(species, chrom, "", start, end)
        self.size[0] = end - start + 1

    def get_view(self):
        return self.view
    
    def update(self):
        """update graphics"""
    
    def draw(self):
        """Returns initial graphics for track"""
        return group()
    
    def draw_left(self):
        """Returns initial graphics for track that should appear in the left window"""
        return group()
    
    def draw_top(self):
        """Returns initial graphics for track that should appear in the top window"""
        return group()
    
    def get_window(self):
        """Gets the main browser window"""
        
        # get window from browser
        return self.browser.get_window()

    def set_browser(self, browser):
        self.browser = browser
    
    def getBrowser(self):
        return self.browser




class GenomeStackBrowser (Browser):
    """A browser where tracks are stacked vertically"""

    def __init__(self, view=None, winsize=(800,400), winpos=None,
                 **options):
        Browser.__init__(self, **options)
        
        self.win = None         # main window
        self.winsize = winsize
        self.winpos = winpos
        self.gid = None
        self.leftgid = None
        
        # create stack track
        self.stack = StackTrack()
        if view != None:
            self.stack.set_view(view.species, view.seqname, view.start, view.end)
        Browser.add_track(self, self.stack)
    
    
    def set_view(self, species, chrom, start, end):
        self.stack.set_view(species, chrom, int(start), int(end))
        
    
    def add_track(self, track, index=None):
        self.stack.add_track(track, index)

    def remove_track(self, track):
        self.stack.remove_track(track)
        
    
    def clear_tracks(self):
        self.stack.clear_tracks()
    
    
    def show(self, species=None, chrom=None, start=None, end=None,
                   width=800, height=400):
        self.enable_updating(False)
        
        view_change = False
        
        # initialize window if needed
        if self.win == None:
            self.win = summon.Window("StackBrowser", 
                                     size=self.winsize, position=self.winpos)
            self.win.set_bgcolor(1, 1, 1)
            
            # TODO: add binding query for zooming
            self.win.set_binding(input_click("right", "down"), "focus")
            self.win.set_binding(input_motion("right", "down"), "zoomx")
            self.win.set_binding(input_click("right", "down", "shift"), "focus")
            self.win.set_binding(input_motion("right", "down", "shift"), "zoomy")
            
            # initialize root group
            self.gid = self.win.add_group(group())
            
            view_change=True
        
        # initialize left window if needed
        if self.show_left_window:
            if self.leftwin == None:
                self.leftwin = summon.Window(" ", size=(150, height))
                self.leftwin.set_bgcolor(1, 1, 1)

                # initialize root group                
                self.leftgid = self.leftwin.add_group(group())
                
                self.left_ensemble = multiwindow.WindowEnsemble(
                                     [self.leftwin, self.win], 
                                      stacky=True, sameh=True,
                                      tiey=True, piny=True,
                                      master=self.win)
                self.leftwin.set_boundary(0, util.INF, -150, -util.INF)
        else:
            if self.leftwin:
                self.leftwin.close()
                self.leftwin = None
        
        
        # initialize view
        if species != None:
            view_change = True
            self.set_view(species, chrom, start, end)
                
        
        # perform drawing 
        self.gid = self.win.replace_group(self.gid, self.stack.draw())
        if self.leftwin:
            self.leftid = self.leftwin.replace_group(self.leftgid, self.stack.draw_left())
        
        
        # center window view
        if view_change:
            stacksize = self.stack.get_size()
            stackpos = self.stack.pos
            
            w, h = self.win.get_size()
            self.win.set_visible(stackpos[0], stackpos[1],
                                 stackpos[0] + stacksize[0], stackpos[1]-stacksize[1])
            self.win.focus(w/2, h/2)
            self.win.zoom(1, stacksize[0] / (5.0 * stacksize[1]))
            self.win.zoom(.9, .9)
        
        self.enable_updating()

    




class GenomeOverview (Browser):
    """Get an overview of elements across the genome"""

    def __init__(self, chroms=None, chrom_step=-3, 
                       winsize=(400,400), winpos=None,
                       **options):
        """Initialize browser
        
           chroms    -- chromosomes to display (Region objects)
           chrom_step -- separation between chromosomes (default: -3)
        """
        
        Browser.__init__(self, **options)
        self.win = None
        self.leftwin = None
        self.metatracks = []
        self.show_ruler = True
        self.winsize = winsize
        self.winpos = winpos
        
        
        if chroms == None:
            self.chroms = []
        else:
            self.chroms = chroms
        
        # layout chromosomes
        self.chrom_step = chrom_step
                
        self.chrom_pos = {}
        y = 0
        x = 0
        for chrom in self.chroms:
            self.chrom_pos[chrom] = [x, y]
            y += self.chrom_step
            
        self.rulers = self.add_track(RulerTrack, 0, top=1.0, bottom=-1.0,
                                     text_align="top")
                

    def add_regions(self, regions, **options):
        """Add regions to browser
        
           options:
            col    -- color (default: color(0,0,1))
            style  -- drawing style (default: "box")
            height -- height with respect to ruler (default: 0.5)
        """
    
        options.setdefault("col", color(0,0,1))
        options.setdefault("style", "box")
        options.setdefault("height", 0.5)
        
        return self.add_track(RegionTrack, -.5, regions, **options)
          
    
    
    def show(self):
        """Display browser"""
        self.enable_updating(False)
        
        # create window
        if self.win == None or not self.win.is_open():
            newwin = True
            self.win = summon.Window(position=self.winpos, size=self.winsize)
                
            self.win.set_bgcolor(1, 1, 1)
            self.win.select.set_callback(self.on_select)
        

            # zooming
            self.win.set_binding(input_click("right", "down"), "focus")
            self.win.set_binding(input_motion("right", "down"), "zoomx")
            self.win.set_binding(input_click("right", "down", "shift"), "focus")
            self.win.set_binding(input_motion("right", "down", "shift"), "zoomy")

        
            # create left window
            pos = self.win.get_position()
            self.leftwin = summon.Window(" ", size=(150, self.win.get_size()[1]),
                      position=(pos[0]-150-self.win.get_decoration()[0], pos[1]))
                
            self.leftwin.set_bgcolor(1, 1, 1)
        
            self.left_ensemble = multiwindow.WindowEnsemble(
                                 [self.leftwin, self.win], 
                                  stacky=True, sameh=True,
                                  tiey=True, piny=True,
                                  master=self.win)
            
            # menu
            self.sidebar = hud.SideBar(self.win, width=150)
            self.sidebar.add_item(hud.MenuItem("toggle ruler (r)",
                                               self.toggle_ruler))
            self.win.set_binding(input_key("r"), self.toggle_ruler)
            
        else:
            newwin = False
            self.win.clear_groups()
        
        # determine largest chromosome size
        maxchrom = max(x.length() for x in self.chroms)
        
        # add line to left window
        maxname = max(len(x.seqname) for x in self.chroms) * 20
        self.leftwin.add_group(
            lines(color(1,1,1), 
                  0, 0, -maxname, 0,
                  0, 0, 0, self.chrom_step * len(self.chroms),
                  0, self.chrom_step * len(self.chroms), 
                  -maxname, self.chrom_step * len(self.chroms)))
        
        
        # draw chromosome labels
        self.win.add_group(self.draw_chrom_labels())
        self.leftwin.add_group(self.draw_chrom_labels())
        
        # draw all tracks
        for track in self._tracks:
            self.win.add_group(track.draw())
        

        if newwin:    
            # setup left window
            self.leftwin.home()
        
            # setup window
            w, h = self.win.get_size()
            bottom = float(min(self.chrom_pos[i][1] for i in self.chroms) +
                           2*self.chrom_step)
        

            self.win.set_visible(0, -2*self.chrom_step, maxchrom, bottom)
            self.win.focus(w/2, h/2)
            self.win.zoom(1, maxchrom / abs(-2*self.chrom_step - bottom))
            self.win.zoom(.9, .9)

        self.enable_updating(True)
    
    
    def draw_chrom_labels(self):
        """draw chromosome labels"""
        
        vis = [color(0, 0, 0)]
        
        # determine largest chromosome size
        maxchrom = max(x.length() for x in self.chroms)
        
        # draw labels
        y = 0
        for chrom in self.chroms:
            x, y = self.chrom_pos[chrom]
            
            # draw chromosome name (main window)
            vis.append(text_clip(chrom.seqname + " ", 
                                 -maxchrom, y-self.chrom_step/2.0, 
                                 0, y+self.chrom_step/2.0,
                                 4, 20, 
                                 "right", "middle"))
        return group(*vis)
       
    
    def add_track(self, track_class, offset, *args, **kargs):
        tracks = []
    
        # make a track for each chromosome
        for chrom in self.chroms:
            x, y = self.chrom_pos[chrom]
            
            track = track_class(*args, **kargs)
            track.set_pos(x, y+offset)
            track.set_view(chrom.species, 
                          chrom.seqname,
                          chrom.start,
                          chrom.end)
            Browser.add_track(self, track)
            tracks.append(track)
        
        metatrack = (track_class, offset, args, kargs, tracks)
        self.metatracks.append(metatrack)
        return metatrack
    
    
    def remove_track(self, metatrack):
        (track_class, offset, args, kargs, tracks) = metatrack
        
        # remove meta track
        self.metatracks.remove(metatrack)
        
        # remove tracks
        for track in tracks:
            Browser.remove_track(self, track)
    
    
    def pos2chrom(self, x, y):
        """Returns a chromosome that intersects the given point (x, y)"""
    
        for chrom in self.chroms:
            pos = self.chrom_pos[chrom]
            
            if pos[0] <= x <= pos[0] + chrom.length() and \
               pos[1] -1 <= y <= pos[1] + 1:
                return chrom
        return None
             
             
    
    def on_select(self, pos1, pos2):
        """Callback for region selections"""
         
        # determine selected chrom
        chrom = self.pos2chrom(*pos1)
        chrom2 = self.pos2chrom(*pos2)

        if chrom == None:
           print "no chromosome selected"
           return

        if chrom != chrom2:
           print "please select only one chromosome"
           return

        x1 = pos1[0]
        x2 = pos2[0]
        if x1 > x2:
           x1, x2 = x2, x1
        chrom_pos = self.chrom_pos[chrom][0]

        newchrom = Region(chrom.species, chrom.seqname, chrom.feature, 
                          int(max(chrom.start + x1 - chrom_pos, 0)),
                          int(min(chrom.start + x2 - chrom_pos, chrom.end)))

        print "displaying region: %s:%s:%s-%s:%s" % \
           (newchrom.species, newchrom.seqname, 
            util.int2pretty(newchrom.start),
            util.int2pretty(newchrom.end),
            ["+", "+", "-"][newchrom.strand])

        # create sub-browser
        subbrowser = GenomeOverview([newchrom],
                                    winpos=self.win.get_position(),
                                    winsize=self.win.get_size())

        for track_class, offset, args, kargs, tracks in self.metatracks:
            subbrowser.add_track(track_class, offset, *args, **kargs)
         
        subbrowser.show()
         

    def toggle_ruler(self):
        self.show_ruler = not self.show_ruler
        
        for ruler in self.rulers[4]:
            ruler.set_visible(self.show_ruler)



class StackTrack (Track):
    """This track manages a stack of child tracks"""

    def __init__(self, *args, **options):
        Track.__init__(self, *args, **options)
        self._tracks = []
    
    def set_view(self, species, chrom, start, end):
        # set personal view
        Track.set_view(self, species, chrom, start, end)
        
        # propogate view change to children tracks
        for track in self._tracks:
            track.set_view(species, chrom, start, end)
    
    def add_track(self, track, index=None):
    
        # propogate browser
        track.set_browser(self.getBrowser())
        if index == None:
            self._tracks.append(track)
        else:
            self._tracks.insert(index, track)
    
    def remove_track(self, track):
        track.set_browser(None)
        self._tracks.remove(track)
    
    def clear_tracks(self):
        for track in self._tracks:
            track.browser = None
        self._tracks = []
        
    def set_browser(self, browser):
        Track.set_browser(self, browser)
        
        # propogate browser to children
        for track in self._tracks:
            track.set_browser(browser)

    def update(self):
    
         # propogate updates to child tracks
        for track in self._tracks:
            track.update()

    def get_size(self):
        self._layout()
        return self.size
        
    
    def _layout(self):
        """Layout the child tracks vertically"""
        top = 2.0
        y = 0.0
        x = 0.0
        maxend = 0.0
        for track in self._tracks:
            track.set_view(self.view.species, self.view.seqname, 
                          self.view.start, self.view.end)
            tracksize = track.get_size()
            y -= tracksize[1]
            track.set_pos(track.pos_offset[0] + x, track.pos_offset[1] + y)
        self.size = [self.view.end - self.view.start + 1, 0 - y]
            
    
    def draw(self):
        self._layout()
        
        vis = group()
        for track in self._tracks:
            vis.append(track.draw())
        return vis
        
    
    def draw_left(self):
        self._layout()
        
        vis = group()
        for track in self._tracks:
            vis.append(track.draw_left())
        return vis


class DividerTrack (Track):
    """Track for visually separating other tracks
    
       This track is especailly useful with the GenomeStackBrowser
    """

    def __init__(self, top_color=color(0,0,0,0),
                       bottom_color=color(0,0,0,0),
                       height=1,
                       fill_color=color(0,0,0,0),
                       **options):
        Track.__init__(self, **options)
        
        self.size[1] = height
        self.top_color = top_color
        self.bottom_color = bottom_color
        self.fill_color = fill_color
    
    
    def draw(self):
        assert self.view != None, "Track view not initialized"
    
        return group(
            self.bottom_color,
            lines(self.pos[0], self.pos[1],
                  self.pos[0] + self.view.length(), self.pos[1]),
            self.top_color,
            lines(self.pos[0], self.pos[1] + self.size[1],
                  self.pos[0] + self.view.length(), self.pos[1] + self.size[1]),
            self.fill_color,
            quads(self.pos[0], self.pos[1],
                  self.pos[0] + self.view.length(), self.pos[1],
                  self.pos[0] + self.view.length(), self.pos[1] + self.size[1],
                  self.pos[0], self.pos[1] + self.size[1]))
    


class RulerTrack (Track):
    """Track for displaying a base-pair ruler along the genome"""

    def __init__(self, top=2.0, bottom=0.0, 
                 minicolor=color(.8,.8,.8), maincolor = color(0,0,0),
                 align=None,
                 text_align="middle",
                 text_color=color(0, 0, 0),
                 show=True,
                 fixed_height=True,
                 **options):
        
        Track.__init__(self, **options)
        self.top = top
        self.bottom = bottom
        self.minicolor = minicolor
        self.maincolor = maincolor
        self.text_align = text_align
        self.text_color = text_color
        self.show      = show
        self.shown     = show
        self.fixed_height = fixed_height
        
        if align != None:
            self.coords = alignlib.CoordConverter(align)
        else:
            self.coords = None
        
        self.multiscale = Multiscale(marginx=.5, marginy=.5,
                                     scalex=10, scaley=10)
       
    
    def draw(self):
        self.multiscale.init(self.get_window())

        self.start = self.view.start-1
        self.end = self.view.end
        self.height = self.top
        self.multiscale.reset()
        
        self.gid = group()
        return group(self.gid)
        
    
    def update(self):
        self.win = self.get_window()
        
        if not self.show:
            if self.shown:
                self.gid = self.win.replace_group(self.gid, group())
                self.shown = False
        else:
            self.shown = True
            
            if not self.multiscale.same_view():
                if self.coords == None:
                    g = self.draw_ruler(self.pos, 
                                       self.start, 
                                       self.end)
                else:
                    g = self.draw_align_ruler(self.pos, self.start, self.end)
                self.gid = self.win.replace_group(self.gid, g)
    
    
    def set_visible(self, visible):
        self.show = visible
    
    
    def draw_ruler(self, pos, start, end):
        worldx1, worldy1, worldx2, worldy2 = self.win.get_visible()
        screenwidth, screenheight = self.win.get_size()
        
        worldwidth = worldx2 - worldx1
        worldx1 -= worldwidth / 2.0
        worldx2 += worldwidth / 2.0

        # find appropriate unit if one is not given
        unit = visual.getRulerAutoSize(screenwidth, worldwidth)
        order = int(math.log10(unit))
        unit2, unitstr = visual.getUnitSuffix(unit)
        
        
        x, y = pos
        vis = []

        # make mini hashes
        if unit >= 10:
            vis.append(self.minicolor)
            i = unit * (max(start, worldx1 - x + start) // unit)
            while x + i - start <= worldx2 and i < end:
                if i >= start:
                    vis.append(lines(x + i - start, y+self.bottom, 
                                     x + i - start, y+self.top))
                i += unit // 10


        # make main hashes
        i = unit * (max(start, worldx1 - x + start) // unit)
        while x + i - start <= worldx2 and i < end:
            if i >= start:
                vis.append(self.maincolor)            
                vis.append(lines(x + i - start, y, x + i - start, y + self.top))
                vis.append(self.text_color)
                vis.append(text(str(int(i//unit2)) + unitstr, 
                                x + i - start, y, x + i -start - unit, y + self.top, 
                                self.text_align, "right"))
            i += unit

        # base line
        vis.append(lines(self.maincolor, x, y, x + end - start, y))

        return group(* vis)
    
    
    def draw_align_ruler(self, pos, start, end):
        worldx1, worldy1, worldx2, worldy2 = self.win.get_visible()
        screenwidth, screenheight = self.win.get_size()
        
        worldwidth = worldx2 - worldx1
        worldx1 -= worldwidth / 2.0
        worldx2 += worldwidth / 2.0

        # find appropriate unit if one is not given
        unit = visual.getRulerAutoSize(screenwidth, worldwidth)
        order = int(math.log10(unit))
        unit2, unitstr = visual.getUnitSuffix(unit)
        
        x, y = pos
        vis = []
        
        # make mini hashes
        vis.append(self.minicolor)
        i = unit * (max(start, worldx1 - x + start) // unit)
        while x + i - start <= worldx2 and i < end:
            if i >= start:
                vis.append(lines(x + i - start, y + self.bottom, 
                                 x + i - start, y + self.height))
            i += max(unit // 10, 1)


        # make main hashes
        
        # find starting local coord
        seqi = unit * ((start + self.coords.align2local(max(0, worldx1 - x), 
                                                        clamp=True)) // unit) \
                                                        - start-1
        # find starting align coord
        i = self.coords.local2align(seqi)
        endseqi = min(self.coords.align2local(end, clamp=True), 
                      self.coords.align2local(worldx2-x, clamp=True))
        
        # draw all hashes in view
        while seqi <= endseqi:
            vis.append(self.maincolor)
            vis.append(lines(x + i+1, y, x + i+1, y + self.height))
            vis.append(self.text_color)
            vis.append(text(str(int((seqi+start+1)//unit2)) + unitstr, 
                            x + i+1, y, x + i+1 - unit, y + self.height, 
                            self.text_align, "right"))
            seqi += unit
            i = self.coords.local2align(seqi, clamp=True)

        # base line
        vis.append(lines(self.maincolor, x, y, x + end - start, y))

        return group(* vis)



class SeqTrack (Track):
    """Track for displaying a single sequence along the genome"""

    def __init__(self, seqs, **options):
        options.setdefault("size", [0.0, 0.5])
        Track.__init__(self, **options)
        
        self.seqs = seqs
        self.shown = False
        self.multiscale = Multiscale(marginx=.5, marginy=.5)

    
    def draw(self):
        self.shown = False
        self.multiscale.init(self.get_window())
        
        # return placeholder group
        self.gid = group()
        return self.gid
    
    
    def update(self):
        if self.view == None:
            raise Exception("Track view not initialized")
        
        win = self.get_window()
        
        view = win.get_visible()
        x, y = self.pos
        
        if self.multiscale.atleast(4, .1, view=view):
            if not self.shown or \
               not self.multiscale.same_view(view):
                self.shown = True
                start = max(int(self.multiscale.worldx1 - x + self.view.start), 
                            int(self.view.start))
                end = min(int(self.multiscale.worldx2 - x + self.view.start), 
                          int(self.view.end))
                
                seq = self.seqs.getseq(self.view.seqname, start, end)
                
                # convert from inclusive to exclusive
                end = len(seq) + start - 1
                
                self.gid = win.replace_group(self.gid, 
                    group(translate(x, y,
                        color(0,0,0),
                        scale(1, self.size[1], 
                        text_scale(seq, start-self.view.start, 
                                        0, end-self.view.start+1, 2, 
                                        "left", "bottom")))))
            
        elif self.shown:
            self.shown = False
            self.gid = win.replace_group(self.gid, group())


class CurveTrack (Track):
    """Track for displaying a curve along the genome"""

    def __init__(self, xdata, ydata, **options):
        Track.__init__(self, **options)
        
        self.xdata = xdata
        self.ydata = ydata
        self.multiscale = Multiscale(marginx=.25, marginy=.25, 
                                            scalex=4.0, scaley=4.0)
        self.shown = False
        
    
    def draw(self):
        self.multiscale.init(self.get_window())
        self.shown = False
        
        # return placeholder group
        self.gid = group()
        return self.gid
        
    
    def update(self):
        assert self.view != None, "Track view not initialized"
        
        win = self.get_window()
    
        if not self.shown or not self.multiscale.same_view():
            self.shown = True
            x, y = self.pos
            
            start2 = self.view.start
            start = int(max(0, start2-1, self.multiscale.worldx1-x+start2))
            end = int(min(len(self.xdata), self.view.end, self.multiscale.worldx2-x+start2))
            step = max(1, (end - start) // 400)
            
            vis = []
            vis2 = []
            for i in xrange(start, end, step):
                dat = self.ydata[i:i+step]

                assert len(dat) > 0, (start, end, step)

                y1 = min(dat)
                y2 = max(dat)
                #y = self.ydata[i]
                y1 = (util.clamp(y1, .33, .66) - .33) / .33
                y2 = (util.clamp(y2, .33, .66) - .33) / .33                
                vis.extend([self.xdata[i], y2 * self.size[1]])
                vis2.extend([self.xdata[i], y1 * self.size[1]])

            # draw curve on middle of base (.5)
            self.gid = win.replace_group(self.gid, 
                group(translate(x - self.view.start + .5, y,
                      color(0,1,0), 
                      line_strip(* vis),
                      line_strip(* vis2),
                      color(0,0,0),
                      lines(self.view.start - 0.5, 0, self.view.end + 0.5, 0,
                            self.view.start - 0.5, self.size[1], 
                            self.view.end + 0.5, self.size[1]))))




class RegionTrack (Track):
    """Track for displaying genomic regions (genes, regulator elements, etc)"""

    def __init__(self, regions, height=0.5, col=color(0,0,1,.5), 
                       text_color=color(1, 1, 1),
                       textSize=12,
                       style="box", on_click=None,
                       **options):
        Track.__init__(self, **options)
        
        self.regions = regions
        self.color = col
        self.text_color = text_color
        self.textSize = textSize
        self.style = style
        self.height = height
        self.on_click = on_click
    
    
    def get_region_pos(self, reg):
        if reg.seqname == self.view.seqname and \
           util.overlap(self.view.start, self.view.end, reg.start, reg.end):
            return (self.pos[0] + reg.start - self.view.start, self.pos[1])
        else:
            return None
    
    
    def draw(self):
        assert self.view != None, "Track view not initialized"
    
        species = self.view.species
        chrom = self.view.seqname
        start = self.view.start
        end = self.view.end
    
        height = self.height
        regions = filter(lambda x: x.seqname == chrom and 
                                   x.species == species, self.regions)
        
        
        def click_region(region):
            return lambda: self.on_click(region)
        
        if self.style == 'box':
            # box style
            vis = []
            vis2 = []
            hotspots = []
            names = []
        
            for reg in regions:
                if util.overlap(start, end, reg.start, reg.end):
                    if reg.strand == 1:
                        # positive strand
                        bot = 0
                        top = height
                    elif reg.strand == -1:
                        # negative strand
                        bot = -height
                        top = 0
                    else:
                        bot = -height
                        top = height
                        
                    vis.extend([reg.start-start, bot,
                                reg.end-start+1, bot,
                                reg.end-start+1, top,
                                reg.start-start, top])
                    vis2.extend([reg.start-start, bot, 
                                 reg.start-start, top])
                    if 'ID' in reg.data:
                        names.append(text_clip(reg.data['ID'], 
                                               reg.start-start, bot,
                                               reg.end-start+1, top,
                                               4, 
                                               self.textSize))
                    if self.on_click:
                        hotspots.append(hotspot("click",
                                                reg.start-start, top,
                                                reg.end-start+1, bot,
                                                click_region(reg)))
            
            return group(translate(
                self.pos[0], self.pos[1] + self.size[1] / 2.0,
                color(0,0,0), lines(0, 0, end-start+1, 0),
                self.color, quads(*vis), lines(*vis2),
                group(*hotspots),
                self.text_color, *names))
                
        elif self.style == 'line':
            # line style
            vis = []
        
            for reg in regions:
                if util.overlap(start, end, reg.start, reg.end):
                   if reg.strand != 0:
                       vis.extend([reg.start-start, 0, 
                                   reg.start-start, reg.strand*height])
                   else:
                       vis.extent([reg.start-start, -.5 * height, 
                                   reg.start-start, .5 * height])
            
            return group(translate(self.pos[0], self.pos[1]  + self.size[1] / 2.0,
                                   color(0,0,0), lines(0, 0, end-start+1, 0),
                                   self.color, lines(* vis)))
        
        else:
            # custom style
            regions2 = (reg for reg in regions 
                        if util.overlap(start, end, reg.start, reg.end))
            
            return group(translate(self.pos[0], self.pos[1] + self.size[1] / 2.0, 
                         self.style(regions2, start)))
        



class GenomeAlignTrack (Track):
    def __init__(self, galign, collapse=None, color_bases=False, 
                 height=1,
                 seqtype="dna", show_color_bases=True, show_bases=True, **options):
        Track.__init__(self, **options)
        
        self.galign = galign
        self.collapse = collapse
        self.color_bases = color_bases
        self.seqtype = seqtype
        self.show_color_bases = show_color_bases
        self.show_bases = show_bases
        
        self.size = [None, None]
        self.aligns = None
        self.tracks = []

    def set_view(self, species, chrom, start, end):
        Track.set_view(self, species, chrom, start, end)
        self.aligns = self.galign.get_aligns(species, chrom, start, end,
                                             collapse=self.collapse)
    
    def get_size(self):
        assert self.view != None
        
        if self.aligns == None:
            self.aligns = self.galign.get_aligns(self.view.species,
                                                 self.view.seqname,
                                                 self.view.start, self.view.end,
                                                 collapse=self.collapse)
        x = 0           
        for aln in self.aligns:
            x += aln.alignlen()
        self.size = [x, len(aln)]
        return self.size
                
    
    def draw(self):
        if self.view == None:
            raise Exception("view must be set before calling draw()")
        
        self.tracks = []
        vis = []
        x = 0
        height = 0
        show_labels = True
        
        for aln in self.aligns:
            self.tracks.append(AlignTrack(aln,
                                          color_bases=self.color_bases, 
                                          seqtype=self.seqtype, 
                                          show_color_bases=self.show_color_bases, 
                                          show_bases=self.show_bases,
                                          show_labels=show_labels))
            self.tracks[-1].set_pos(self.pos[0] + x, self.pos[1])
            self.tracks[-1].browser = self.browser
            x += aln.alignlen()
            height = max(height, len(aln))
            vis.append(self.tracks[-1].draw())
            
            # only the first track needs labels
            show_labels = False
        
        return group(*vis)
    
    
    def draw_left(self):
        if len(self.tracks) > 0:
            return self.tracks[0].draw_left()
        else:
            return group()
    
    
    def update(self):
        for track in self.tracks:
            track.update()
        

class AlignTrack (Track):
    def __init__(self, aln, collapse=None, cols=None, color_bases=False, 
                 seqtype=None, show_color_bases=True, show_bases=True,
                 show_labels=True,
                 rowspacing=None,
                 **options):
        Track.__init__(self, **options)
        self.size = [aln.alignlen(), len(aln)]
        self.multiscale = Multiscale(marginx=.5, marginy=.5)
        self.collapse = collapse
        
        self.show_color_bases = show_color_bases
        self.show_bases = show_bases
        self.show_labels = show_labels
        self.rowspacing = rowspacing
        self.always_color = False
        self.color_bases_vis = None
        
        if seqtype == None:
            self.seqtype = guessAlign(aln)
        else:
            self.seqtype = seqtype
        
        if color_bases == True:
            if self.seqtype == "dna":
                self.color_bases = dna_colors
            elif self.seqtype == "pep":
                self.color_bases = pep_colors
        else:
            self.color_bases = color_bases
        
        if collapse != None:
            assert collapse in aln.keys()
            cols = util.findneq('-', aln[collapse])
            
        if cols != None:
            self.aln = alignlib.subalign(aln, cols)
        else:
            self.aln = aln
        
    
        
    def draw(self):
        self.text_shown = False
        self.multiscale.init(self.get_window())
    
        BASE    = 0
        GAP     = 1
        NOBASE  = 2
        
        if self.seqtype == "dna":
            baseclasses = {'A': BASE, 'C': BASE, 'T': BASE, 'G': BASE,
                           '-': GAP, 'N': NOBASE, '*': NOBASE, 'X': NOBASE}
            
        elif self.seqtype == "pep":
            baseclasses = {'-': GAP, 'N': NOBASE, '*': NOBASE, 'X': NOBASE}
            
            for aa in 'ARNDCEQGHILKMFPSTWYVU':
                baseclasses[aa] = BASE
        
        else:
            raise Exception("unknown seqtype '%s'" % self.seqtype)
        
        # init row spacing
        if self.rowspacing == None:
            self.rowspacing = range(len(self.aln))
        
        def getRegions(selectedClass):
            boxpts = []
            diagpts = []
            diagpts2 = []
            
            for row, (key, val) in zip(self.rowspacing, self.aln.iteritems()):
                lastbase = None
                lastclass = None
                lasti = 0
                for i in xrange(len(val)+1):
                    # this extra is being used to handle the case when
                    # a sequence is all bases
                    if i < len(val):
                        base = val[i]
                    else:
                        base = '-'
                    if base not in baseclasses:
                        baseclass = NOBASE
                    else:
                        baseclass = baseclasses[base]

                    if baseclass == lastclass:
                        continue
                    
                    if lastbase != None and lastclass == selectedClass:
                        boxpts.extend([lasti, -row, lasti, -row-1, i, -row-1, i, -row])
                        diagpts.extend([i, -row, i, -row-1])
                        #diagpts.extend([lasti, -row, i, -row-1])
                        #diagpts2.extend([lasti, -row, lasti, -row-1])

                    lasti = i
                    lastbase = base
                    lastclass = baseclass
            return boxpts, diagpts, diagpts2
        
        base_boxpts, base_diagpts, base_diagpts2 = getRegions(BASE)
        nobase_boxpts, nobase_diagpts, nobase_diagpts2 = getRegions(NOBASE)
        
        
        # build labels
        if self.show_labels:
            labels = []
            for i, key in zip(self.rowspacing, self.aln):
                labels.append(text_clip(key, -util.INF, -i, 0, -i-1,
                                    4, 12, "middle", "right"))
            labelsgroup = group(*labels)
        else:
            labelsgroup = group()
        
        
        # build hotspot
        click = hotspot("click", 0, 0, self.aln.alignlen(), -self.size[1],
                        self.on_click_callback)
        
        self.text_group = group()
        
        return group(translate(self.pos[0], self.pos[1] + self.size[1],
                     color(0, 0, 0),
                     labelsgroup,
                     
                     click,
                     
                     color(.5, .5, .5), 
                     quads(* base_boxpts),
                     lines(* base_diagpts),
                     #lines(* base_diagpts2),
                     
                     color(.7, .2, .2),
                     quads(* nobase_boxpts),
                     lines(* nobase_diagpts),
                     #lines(* nobase_diagpts2),
                     group(self.text_group)))
    
    def draw_left(self):
        labels = []
        maxsize = max(map(len, self.aln.keys()))
        
        for i, key in enumerate(self.aln):
            labels.append(text_clip(key, -util.INF, -i, 0, -i-1,
                                    4, 12, "middle", "right"))
        return group(translate(self.pos[0], self.pos[1] + self.size[1],
                     color(0,0,0), 
                     #lines(0, 0, 0, -len(self.aln),
                     #      -maxsize, 0, 0, 0,
                     #      -maxsize, -len(self.aln), 0, -len(self.aln)),
                     *labels))
        
        
    
    def update(self):
        win = self.get_window()
        view = win.get_visible()
        size = win.get_size()
        x, y = self.pos
        
        mintextSize = 4
        minblockSize = 1
        
        color_bases = self.color_bases
        
        if self.multiscale.atleast(minblockSize, .1, 
                                   view=view, size=size):
            if not self.text_shown or \
               not self.multiscale.same_view(view):
                self.text_shown = True
                
                start = max(int(self.multiscale.worldx1 - x - 1), 
                            0)
                end = max(int(self.multiscale.worldx2 - x), start)
                
                vis = []
                vis2 = []
                for i, row in enumerate(self.aln.itervalues()):
                    seq = row[start:end]
                    seq = seq.replace("-", " ")
                    
                    # color bases
                    if self.show_color_bases and color_bases != False:
                        for j in xrange(len(seq)):
                            base = seq[j].upper()
                            if base in color_bases:
                                vis2.extend([color_bases[base], 
                                             quads(start + j, -i,
                                                   start + j, -i-1,
                                                   start + j+1, -i-1,
                                                   start + j+1, -i)])
                    
                    end2 = start + len(seq)
                    
                    # draw text
                    if self.show_bases and \
                       self.multiscale.atleast(mintextSize, 2, 
                                               view=view, size=size):
                        vis.append(text_scale(seq, 
                                              start, -i+1, 
                                              end2, -i-1, 
                                              "left", "bottom"))
                
                self.text_group = win.replace_group(self.text_group, 
                    group(group(*vis2), color(0,0,0), * vis))
            
        elif self.text_shown:
            self.text_shown = False
            self.text_group = win.replace_group(self.text_group, group())


    def on_click_callback(self):
        # TODO: should not rely on get_mouse_pos("world")
        win = self.get_window()
        x, y = win.get_mouse_pos('world')
        x -= self.pos[0]
        y = self.size[1] - (y - self.pos[1])
        self.on_click(x, y)


    def on_click(self, x, y):
        y = int(y)
    
        if 0 <= y < len(self.aln):
            print self.aln.keys()[y]
        
        

class TrackOverlay (Track):
    """Overlay several tracks on top of one another to create a single track"""
    
    def __init__(self, tracks, **options):
        Track.__init__(self, **options)
        
        self.tracks = tracks
        

    
    
    def set_view(self, species, chroms, start, end):
        for track in self.tracks:
            track.set_view(species, chroms, start, end)
    
    
    def get_size(self):
        maxheight = 0
        maxwidth = 0
        for track in tracks:
            w, h = track.get_size()
            maxwidth = max(maxwidth, h)            
            maxheight = max(maxheight, h)
        self.size = [maxwidth, maxheight]
        return self.size
    
    
    def draw(self):
        """Returns initial display graphics"""
    
        vis = []
        for track in self.tracks:
            track.view = self.view
            vis.append(track.draw())
        return group(translate(self.pos[0], self.pos[1], * vis))
    
    
    def update(self):
        """Updates track periodically if needed"""
    
        for track in self.tracks:
            track.update()



def show_align(* alns):
    """Quick way to visualize an alignment"""
    
    def colorAlign(aln):
        if guessAlign(aln) == "pep":
            return pep_colors
        else:
            return dna_colors
    
    view = Region("", "", "", 1, 1)
    colors = []
    
    height = 0
    for aln in alns:
        view.end = max(view.end, alns[-1].alignlen())
        height += len(alns[-1])
        colors.append(colorAlign(alns[-1]))
    
    browser = GenomeStackBrowser(view=view)
    browser.add_track(RulerTrack(bottom=-height))
    for aln, col in zip(alns, colors):
        browser.add_track(AlignTrack(aln, color_bases=col))
    browser.show()
    
    return browser
    


#=============================================================================
# sequence coloring
#


def prop2color(prop, t=0):
    return {    
    "hydrophobic":          color(1, t, t),
    "weakly hydrophobic":   color(1, .5, t),
    "charged":              color(1, 1, t),
    "polar":                color(t, t, 1),
    "turn":                 color(t, 1, t),
    "met":                  color(t, 1, t),
    "stop":                 color(t, t, .2),
    }[prop]


def make_pep_colors(prop2color=prop2color):
    pep_colors = util.Dict(default=color(.5, .5, .5))

    AA = 'ARNDCEQGHILKMFPSTWYVU*'
    pep_per_prop = util.hist_dict(util.mget(seqlib.AA_PROPERTY, AA))

    prop_counts = util.Dict(default=0)
    for char in AA:
        prop = seqlib.AA_PROPERTY[char]
        tint = prop_counts[prop] / float(pep_per_prop[prop])
        pep_colors[char] = prop2color(prop, tint * .5)
        prop_counts[prop] += 1
    
    return pep_colors


dna_colors = util.Dict({"A": color(1, .5, .5),
                        "T": color(1, 1, .5),
                        "C": color(.5, 1, .5),
                        "G": color(.5, .5, 1)},
                       default=color(.5, .5, .5))

pep_colors = make_pep_colors(prop2color=prop2color)


def guess_seq(seq):
    """Guesses whether a sequence is 'dna' or 'pep'"""
    dna = "ACTG-N"
    
    chars = util.unique(seq.upper())
    
    for char in chars:
        if char not in dna:
            return "pep"
    return "dna"
guessSeq = guess_seq

def guess_align(aln):
    """Guesses whether an alignment is 'dna' or 'pep'"""
    
    if "pep" in [guess_seq(seq) for seq in aln.itervalues()]:
        return "pep"
    else:
        return "dna"
guessAlign = guess_align

'''

def LabelsTrack (Track):
    def __init__(self, labels, color=color(0,0,0), **options):
        Track.__init__(self, **options)
        self.labels = labels
        self.color = color
        self.size = [0, len(labels)]
    
        
    def draw(self):
        # build labels
        labels = []
        for i, key in enumerate(self.aln):
            labels.append(text_clip(key, -1000 * 10, -i, 0, -i-1,
                                    4, 12, "middle", "right"))
                
        return group(translate(self.pos[0], self.pos[1] + self.size[1],
                     self.color,
                     group(*labels)))

'''

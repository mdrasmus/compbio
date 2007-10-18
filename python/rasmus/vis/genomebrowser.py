
import math
import copy

from summon.core import *
import summon
from summon import multiwindow
from summon import select
from summon import hud

from rasmus.vis import visual
from rasmus import util, stats, regionlib
from rasmus.regionlib import Region
from rasmus.bio import gff, fasta, alignlib, seqlib






class Browser (visual.VisObject):
    """Base class for containing and managing multiple graphical elements"""
    
    def __init__(self, showLeftWindow=False, showTopWindow=False):
        visual.VisObject.__init__(self)
        
        self.leftwin = None     # left margin window
        self.topwin = None      # top margin window
        self.showLeftWindow = showLeftWindow
        self.showTopWindow = showTopWindow
        self._tracks = []
        
    
    def addTrack(self, track):
        track.browser = self
        self._tracks.append(track)
    
    
    def removeTrack(self, elm):
        track.browser = None
        self._tracks.remove(track)
    
    def clearTracks(self):
        for track in self._tracks:
            track.browser = None
        self._tracks = []
    
    
    def update(self):
        for track in self._tracks:
            track.update()
    
    def show(self):    
        pass
    
    
    def getLeftWindow(self):
        return self.leftwin
    
    
    def getTopWindow(self):
        return self.topwin

    def enableSideWindows(self, left=None, top=None):
        if left != None:
            self.showLeftWindow = left
        if top != None:
            self.showTopWindow = top
        
        self.show()

# NOTE: I probably don't need to subclass from VisObject
class Track (visual.VisObject):
    def __init__(self, pos=[0.0, 0.0], size=[0.0, 1.0], 
                       posOffset=[0.0, 0.0],
                       view=None):
        visual.VisObject.__init__(self)
        self.height = size[1]
        self.size = size
        self.pos = pos[:]
        self.posOffset = posOffset[:]
        self.view = copy.copy(view)
        self.browser = None
    
    
    def setPos(self, x, y):
        self.pos = [x, y]
    
    
    def setSize(self, w, h):
        self.size = [w, h]
    
        
    def setView(self, species, chrom, start, end):
        self.view = Region(species, chrom, "", start, end)
        
    def draw(self):
        """Returns initial graphics for track"""
        return group()
    
    def drawLeft(self):
        """Returns initial graphics for track that should appear in the left window"""
        return group()
    
    def drawTop(self):
        """Returns initial graphics for track that should appear in the top window"""
        return group()
    
    def get_window(self):
        """Gets the main browser window"""
        
        # get window from browser
        return self.browser.get_window()


    '''
    def show(self, species, chrom, start, end):
        """Draw the track for sequence from the given
           'species' and 'chrom' from base 'start' to 'end'.
        """
        return group()
    '''

class GenomeStackBrowser (Browser):
    """A browser where tracks are stacked vertically"""

    def __init__(self, view=None, **options):
        Browser.__init__(self, **options)
        
        self.tracks = []
        self.win = None         # main window
        self.view = view        # region of genome in view
        self.gid = None
        self.leftgid = None
        
    
    def addTrack(self, track):
        self.tracks.append(track)
        if track.view == None:
            track.view = self.view
        
        # also add track to master list
        Browser.addTrack(self, track)

    
    def show(self, species=None, chrom=None, start=None, end=None,
                   width=800, height=400):
        summon.stop_updating()
        
        viewChange = False
        
        # initialize window if needed
        if self.win == None:
            self.win = summon.Window()
            self.win.set_size(width, height) 
            
            # zooming
            self.win.set_binding(input_click("right", "down"), "focus")
            self.win.set_binding(input_motion("right", "down"), "zoomx")
            self.win.set_binding(input_click("right", "down", "shift"), "focus")
            self.win.set_binding(input_motion("right", "down", "shift"), "zoomy")
            
            # background color
            self.win.set_bgcolor(1, 1, 1)
            
            viewChange=True
        
        # initialize left window if needed
        if self.showLeftWindow:
            if self.leftwin == None:
                self.leftwin = summon.Window(" ")
                self.leftwin.set_bgcolor(1, 1, 1)
                self.leftwin.set_size(150, height)
                
                self.leftEnsemble = multiwindow.WindowEnsemble(
                                     [self.leftwin, self.win], 
                                      stacky=True, sameh=True,
                                      tiey=True, piny=True,
                                      master=self.win)
        else:
            if self.leftwin:
                self.leftwin.close()
                self.leftwin = None
        
        
        # initialize view
        if species != None:
            viewChange = True
            view.species = species
            view.seqname = chrom
            view.start = start
            view.end = end
        
        
        # initialize root groups
        if self.gid == None:
            self.gid = self.win.add_group(group())
        else:
            self.gid = self.win.replace_group(self.gid, group())
        
        if self.leftwin:
            if self.leftgid == None:
                self.leftgid = self.leftwin.add_group(group())
            else:
                self.leftgid = self.leftwin.replace_group(self.leftgid, group())
        else:
            self.leftgid = None
        
        
        
        # draw tracks
        top = 2.0
        y = 0.0
        x = 0.0
        maxend = 0.0
        for track in self.tracks:
            y -= track.size[1]
            track.setPos(track.posOffset[0] + x, track.posOffset[1] + y)
            track.setView(self.view.species, self.view.seqname, 
                          self.view.start, self.view.end)
            
            if track.view != None:
                maxend = max(maxend, track.pos[0] + track.view.end - track.view.start)
            else:
                maxend = max(maxend, track.pos[0] + track.size[0])
            
            
            # perform drawing 
            self.win.insert_group(self.gid, track.draw())
            if self.leftwin:
                self.leftwin.insert_group(self.leftgid, track.drawLeft())
        
        # setup left window
        if viewChange and self.leftwin:
            self.leftwin.home()
        
        # setup window
        if viewChange:
            w, h = self.win.get_size()
            self.win.set_visible(0, y, maxend, top)
            self.win.focus(w/2, h/2)
            self.win.zoom(1, maxend / (5.0 * (top - y)))
            self.win.zoom(.9, .9)
        
        self.setVisible()
        summon.begin_updating()

    




class GenomeOverview (Browser):
    """Get an overview of elements across the genome"""

    def __init__(self, chroms=None, chromStep=-3, **options):
        """Initialize browser
        
           chroms    -- chromosomes to display (Region objects)
           chromStep -- separation between chromosomes (default: -3)
        """
        
        Browser.__init__(self, **options)
        self.win = None
        self.leftwin = None
        self.metatracks = []
        self.showRuler = True
        
        
        if chroms == None:
            self.chroms = []
        else:
            self.chroms = chroms
        
        # layout chromosomes
        self.chromStep = chromStep
                
        self.chromPos = {}
        y = 0
        x = 0
        for chrom in self.chroms:
            self.chromPos[chrom] = [x, y]
            y += self.chromStep
            
        self.rulers = self.addTrack(RulerTrack, 0, top=1.0, bottom=-1.0, textAlign="top")
                

    def addRegions(self, regions, **options):
        """Add regions to browser
        
           options:
            col    -- color (default: color(0,0,1))
            style  -- drawing style (default: "box")
            height -- height with respect to ruler (default: 0.5)
        """
    
        options.setdefault("col", color(0,0,1))
        options.setdefault("style", "box")
        options.setdefault("height", 0.5)
        
        return self.addTrack(RegionTrack, -.5, regions, **options)
          
    
    
    def show(self, winsize=(400,400), winpos=(-1, -1)):
        """Display browser"""
        self.enableUpdating(False)
        
        # create window
        if self.win == None or not self.win.is_open():
            newwin = True
            try:
                raise
                self.win = summon.Window(position=winpos, size=winsize)
            except:
                # this is for back compatiability to summon-1.8.2
                self.win = summon.Window()
                self.win.set_position(*winpos)
                self.win.set_size(*winsize)
                
            self.win.set_bgcolor(1, 1, 1)
            self.select = select.Select(self.win, self.onSelect, fillColor=(0,0,1,.3))
        

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
        
            self.leftEnsemble = multiwindow.WindowEnsemble(
                                 [self.leftwin, self.win], 
                                  stacky=True, sameh=True,
                                  tiey=True, piny=True,
                                  master=self.win)
            
            # menu
            self.sidebar = hud.SideBar(self.win, width=150)
            self.sidebar.addItem(hud.MenuItem("toggle ruler (r)", self.toggleRuler))
            self.win.set_binding(input_key("r"), self.toggleRuler)
            
        else:
            newwin = False
            self.win.clear_groups()
        
        # determine largest chromosome size
        maxchrom = max(x.length() for x in self.chroms)
        
        # add line to left window
        maxname = max(len(x.seqname) for x in self.chroms) * 20
        self.leftwin.add_group(lines(color(1,1,1), 
                                     0, 0, -maxname, 0,
                                     0, 0, 0, self.chromStep * len(self.chroms),
                                     0, self.chromStep * len(self.chroms), 
                                     -maxname, self.chromStep * len(self.chroms)))
        
        
        # draw chromosome labels
        self.win.add_group(self.drawChromLabels())
        self.leftwin.add_group(self.drawChromLabels())
        
        # draw all tracks
        for track in self._tracks:
            self.win.add_group(track.draw())
        

        if newwin:    
            # setup left window
            self.leftwin.home()
        
            # setup window
            w, h = self.win.get_size()
            bottom = float(min(self.chromPos[i][1] for i in self.chroms) +
                           2*self.chromStep)
        

            self.win.set_visible(0, -2*self.chromStep, maxchrom, bottom)
            self.win.focus(w/2, h/2)
            self.win.zoom(1, maxchrom / abs(-2*self.chromStep - bottom))
            self.win.zoom(.9, .9)

        self.enableUpdating(True)
        summon.begin_updating()
    
    
    def drawChromLabels(self):
        """draw chromosome labels"""
        
        vis = [color(0, 0, 0)]
        
        # determine largest chromosome size
        maxchrom = max(x.length() for x in self.chroms)
        
        # draw labels
        y = 0
        for chrom in self.chroms:
            x, y = self.chromPos[chrom]
            
            # draw chromosome name (main window)
            vis.append(text_clip(chrom.seqname + " ", 
                                 -maxchrom, y-self.chromStep/2.0, 
                                 0, y+self.chromStep/2.0,
                                 4, 20, 
                                 "right", "middle"))
        return group(*vis)
       
    
    def addTrack(self, trackClass, offset, *args, **kargs):
        tracks = []
    
        for chrom in self.chroms:
            x, y = self.chromPos[chrom]
                                    
            # build ruler for this chromosome
            track = trackClass(*args, **kargs)
            track.setPos(x, y+offset)
            track.setView(chrom.species, 
                          chrom.seqname,
                          chrom.start,
                          chrom.end)
            Browser.addTrack(self, track)
            tracks.append(track)
        
        metatrack = (trackClass, offset, args, kargs, tracks)
        self.metatracks.append(metatrack)
        return metatrack
    
    
    def removeTrack(self, metatrack):
        (trackClass, offset, args, kargs, tracks) = metatrack
        
        # remove meta track
        self.metatracks.remove(metatrack)
        
        # remove tracks
        for track in tracks:
            Browser.removeTrack(self, track)
    
    
    def pos2chrom(self, x, y):
        """Returns a chromosome that intersects the given point (x, y)"""
    
        for chrom in self.chroms:
            pos = self.chromPos[chrom]
            
            if pos[0] <= x <= pos[0] + chrom.length() and \
               pos[1] -1 <= y <= pos[1] + 1:
                return chrom
        return None
             
             
    
    def onSelect(self, pos1, pos2):
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
        chromPos = self.chromPos[chrom][0]

        newchrom = Region(chrom.species, chrom.seqname, chrom.feature, 
                          int(max(chrom.start + x1 - chromPos, 0)),
                          int(min(chrom.start + x2 - chromPos, chrom.end)))

        print "displaying region: %s:%s:%s-%s:%s" % \
           (newchrom.species, newchrom.seqname, 
            util.int2pretty(newchrom.start),
            util.int2pretty(newchrom.end),
            ["+", "+", "-"][newchrom.strand])

        # create sub-browser
        subbrowser = GenomeOverview([newchrom])

        for trackClass, offset, args, kargs, tracks in self.metatracks:
            subbrowser.addTrack(trackClass, offset, *args, **kargs)
         
        subbrowser.show(winpos=self.win.get_position(),
                        winsize=self.win.get_size())
         

    def toggleRuler(self):
        self.showRuler = not self.showRuler
        
        for ruler in self.rulers[4]:
            ruler.setVisible(self.showRuler)


class DividerTrack (Track):
    """Track for visually separating other tracks
    
       This track is especailly useful with the GenomeStackBrowser
    """

    def __init__(self, topColor=color(0,0,0,0),
                       bottomColor=color(0,0,0,0),
                       height=1,
                       fillColor=color(0,0,0,0),
                       **options):
        Track.__init__(self, **options)
        
        self.size[1] = height
        self.topColor = topColor
        self.bottomColor = bottomColor
        self.fillColor = fillColor
    
    
    def draw(self):
        assert self.view != None, "Track view not initialized"
    
        return group(
            self.bottomColor,
            lines(self.pos[0], self.pos[1],
                  self.pos[0] + self.view.length(), self.pos[1]),
            self.topColor,
            lines(self.pos[0], self.pos[1] + self.size[1],
                  self.pos[0] + self.view.length(), self.pos[1] + self.size[1]),
            self.fillColor,
            quads(self.pos[0], self.pos[1],
                  self.pos[0] + self.view.length(), self.pos[1],
                  self.pos[0] + self.view.length(), self.pos[1] + self.size[1],
                  self.pos[0], self.pos[1] + self.size[1]))
    


class RulerTrack (Track):
    """Track for displaying a base-pair ruler along the genome"""

    def __init__(self, top=2.0, bottom=0.0, 
                 minicolor=color(.8,.8,.8), maincolor = color(0,0,0),
                 align=None,
                 textAlign="middle",
                 textColor=color(0, 0, 0),
                 show=True,
                 **options):
        
        Track.__init__(self, **options)
        self.top = top
        self.bottom = bottom
        self.minicolor = minicolor
        self.maincolor = maincolor
        self.textAlign = textAlign
        self.textColor = textColor
        self.show      = show
        self.shown     = show
        
        if align != None:
            self.coords = alignlib.CoordConverter(align)
        else:
            self.coords = None
        
        self.multiscale = visual.Multiscale(marginx=.5, marginy=.5, scalex=10, scaley=10)
       
    
    def draw(self):
        self.multiscale.init(self.get_window())

        self.start = self.view.start-1
        self.end = self.view.end
        self.height = self.top
        self.bottom = self.bottom
        
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
            
            if not self.multiscale.sameScale() or not self.shown:
                if self.coords == None:
                    g = self.drawRuler(self.pos, 
                                       self.start, 
                                       self.end)
                else:
                    g = self.drawAlignRuler(self.pos, self.start, self.end)
                self.gid = self.win.replace_group(self.gid, g)
    
    
    def setVisible(self, visible):
        self.show = visible
    
    
    def drawRuler(self, pos, start, end):
        worldx1, worldy1, worldx2, worldy2 = self.win.get_visible()
        screenwidth, screenheight = self.win.get_size()
        height = self.height
        bottom = self.bottom
        
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
                    vis.append(lines(x + i - start, y+bottom, x + i - start, y+height))
                i += unit // 10


        # make main hashes
        i = unit * (max(start, worldx1 - x + start) // unit)
        while x + i - start <= worldx2 and i < end:
            if i >= start:
                vis.append(self.maincolor)            
                vis.append(lines(x + i - start, y, x + i - start, y + height))
                vis.append(self.textColor)
                vis.append(text(str(int(i//unit2)) + unitstr, 
                                x + i - start, y, x + i -start - unit, y + height, 
                                self.textAlign, "right"))
            i += unit

        # base line
        vis.append(lines(self.maincolor, x, y, x + end - start, y))

        return group(* vis)
    
    
    def drawAlignRuler(self, pos, start, end):
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
            vis.append(self.textColor)
            vis.append(text(str(int((seqi+start+1)//unit2)) + unitstr, 
                            x + i+1, y, x + i+1 - unit, y + self.height, 
                            self.textAlign, "right"))
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
        self.multiscale = visual.Multiscale(marginx=.5, marginy=.5)

    
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
               not self.multiscale.sameScale(view):
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
        self.multiscale = visual.Multiscale(marginx=.25, marginy=.25, 
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
    
        if not self.shown or not self.multiscale.sameScale():
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
                       textColor=color(1, 1, 1),
                       textSize=12,
                       style="box", onClick=None,
                       **options):
        Track.__init__(self, **options)
        
        self.regions = regions
        self.color = col
        self.textColor = textColor
        self.textSize = textSize
        self.style = style
        self.height = height
        self.onClick = onClick
        
    
    def draw(self):
        assert self.view != None, "Track view not initialized"
    
        species = self.view.species
        chrom = self.view.seqname
        start = self.view.start
        end = self.view.end
    
        height = self.height
        regions = filter(lambda x: x.seqname == chrom and 
                                   x.species == species, self.regions)
        
        
        def clickRegion(region):
            return lambda: self.onClick(region)
        
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
                        vis.extend([reg.start-start, 0,
                                    reg.end-start+1, 0,
                                    reg.end-start+1, height,
                                    reg.start-start, height])
                        vis2.extend([reg.start-start, 0, 
                                     reg.start-start, height])
                        if 'ID' in reg.data:
                            names.append(text_clip(reg.data['ID'], 
                                                   reg.start-start, 0,
                                                   reg.end-start+1, height,
                                                   4, 
                                                   self.textSize))
                        if self.onClick:
                            hotspots.append(hotspot("click",
                                                    reg.start-start, height,
                                                    reg.end-start+1, 0,
                                                    clickRegion(reg)))
                    else:
                        # negative strand
                        vis.extend([reg.start-start, -height,
                                    reg.end-start+1, -height,
                                    reg.end-start+1, 0,
                                    reg.start-start, 0])
                        vis2.extend([reg.start-start, 0, 
                                     reg.start-start, -height])
                        if 'ID' in reg.data:                    
                            names.append(text_clip(reg.data['ID'], 
                                                   reg.start-start, -height,
                                                   reg.end-start+1, 0,
                                                   4, 
                                                   self.textSize))
                    
                        if self.onClick:
                            hotspots.append(hotspot("click",
                                                    reg.start-start, -height,
                                                    reg.end-start+1, 0,
                                                    clickRegion(reg)))
            
            return group(translate(self.pos[0], self.pos[1] + self.size[1] / 2.0,
                     color(0,0,0), lines(0, 0, end-start+1, 0),
                     self.color, quads(*vis), lines(*vis2),
                     group(*hotspots),
                     self.textColor, *names))
                
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
    def __init__(self, galign, collapse=None, colorBases=False, 
                 height=1,
                 seqtype=None, showColorBases=True, showBases=True, **options):
        Track.__init__(self, **options)
        
        self.galign = galign
        self.collapse = collapse
        self.colorBases = colorBases
        self.seqtype = seqtype
        self.showColorBases = showColorBases
        self.showBases = showBases
        
        self.size = [0, height]
        self.tracks = []
                
    def draw(self):
        self.aligns = self.galign.getAligns(self.view.species, self.view.seqname,
                                            self.view.start, self.view.end,
                                            collapse=self.collapse)
        
        self.tracks = []
        vis = []
        x = 0
        height = 0
        showLabels = True
        
        for aln in self.aligns:
            self.tracks.append(AlignTrack(aln,
                                          colorBases=self.colorBases, 
                                          seqtype=self.seqtype, 
                                          showColorBases=self.showColorBases, 
                                          showBases=self.showBases,
                                          showLabels=showLabels))
            self.tracks[-1].setPos(self.pos[0] + x, self.pos[1])
            self.tracks[-1].browser = self.browser
            x += aln.alignlen()
            height = max(height, len(aln))
            vis.append(self.tracks[-1].draw())
            
            # only the first track needs labels
            showLabels = False
        
        self.size = [x, height]
        
        return group(*vis)
    
    
    def drawLeft(self):
        if len(self.tracks) > 0:
            return self.tracks[0].drawLeft()
        else:
            return group()
    
    
    def update(self):
        for track in self.tracks:
            track.update()
        

class AlignTrack (Track):
    def __init__(self, aln, collapse=None, cols=None, colorBases=False, 
                 seqtype=None, showColorBases=True, showBases=True,
                 showLabels=True,
                 rowspacing=None,
                 **options):
        Track.__init__(self, **options)
        self.size = [aln.alignlen(), len(aln)]
        self.multiscale = visual.Multiscale(marginx=.5, marginy=.5)
        self.collapse = collapse
        
        self.showColorBases = showColorBases
        self.showBases = showBases
        self.showLabels = showLabels
        self.rowspacing = rowspacing
        
        if seqtype == None:
            self.seqtype = guessAlign(aln)
        else:
            self.seqtype = seqtype
        
        if colorBases == True:
            if self.seqtype == "dna":
                self.colorBases = dna_colors
            elif self.seqtype == "pep":
                self.colorBases = pep_colors
        else:
            self.colorBases = colorBases
        
        if collapse != None:
            assert collapse in aln.keys()
            cols = util.findneq('-', aln[collapse])
            
        if cols != None:
            self.aln = alignlib.subalign(aln, cols)
        else:
            self.aln = aln
        
    
        
    def draw(self):
        self.textShown = False
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
            diagpts3 = []
            
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
                        diagpts.extend([lasti, -row, i, -row-1])
                        diagpts2.extend([lasti, -row, lasti, -row-1])

                    lasti = i
                    lastbase = base
                    lastclass = baseclass
            return boxpts, diagpts, diagpts2
        
        base_boxpts, base_diagpts, base_diagpts2 = getRegions(BASE)
        nobase_boxpts, nobase_diagpts, nobase_diagpts2 = getRegions(NOBASE)
        
        
        # build labels
        if self.showLabels:
            labels = []
            for i, key in zip(self.rowspacing, self.aln):
                labels.append(text_clip(key, -util.INF, -i, 0, -i-1,
                                    4, 12, "middle", "right"))
            labelsgroup = group(*labels)
        else:
            labelsgroup = group()
        
        
        # build hotspot
        click = hotspot("click", 0, 0, self.aln.alignlen(), -self.size[1],
                        self.onClickCallback)
        
        self.textGroup = group()
        
        return group(translate(self.pos[0], self.pos[1] + self.size[1],
                     color(0, 0, 0),
                     labelsgroup,
                     
                     click,
                     
                     color(.5, .5, .5), 
                     quads(* base_boxpts),
                     lines(* base_diagpts),
                     lines(* base_diagpts2),
                     
                     color(.7, .2, .2),
                     quads(* nobase_boxpts),
                     lines(* nobase_diagpts),
                     lines(* nobase_diagpts2),
                     group(self.textGroup)))
    
    def drawLeft(self):
        labels = []
        maxsize = max(map(len, self.aln.keys()))
        
        for i, key in enumerate(self.aln):
            labels.append(text_clip(key, -util.INF, -i, 0, -i-1,
                                    4, 12, "middle", "right"))
        return group(translate(self.pos[0], self.pos[1] + self.size[1],
                     color(0,0,0), 
                     lines(0, 0, 0, -len(self.aln),
                           -maxsize, 0, 0, 0,
                           -maxsize, -len(self.aln), 0, -len(self.aln)),
                     *labels))
        
        
    
    def update(self):
        win = self.get_window()
        view = win.get_visible()
        size = win.get_size()
        x, y = self.pos
        
        mintextSize = 4
        minblockSize = 1
        
        colorBases = self.colorBases
        
        if self.multiscale.atleast(minblockSize, .1, view=view, size=size):
            if not self.textShown or \
               not self.multiscale.sameScale(view):
                self.textShown = True
                
                start = max(int(self.multiscale.worldx1 - x - 1), 
                            0)
                end = max(int(self.multiscale.worldx2 - x), start)
                
                vis = []
                vis2 = []
                for i, row in enumerate(self.aln.itervalues()):
                    seq = row[start:end]
                    seq = seq.replace("-", " ")
                    
                    # color bases
                    if self.showColorBases and colorBases != False:
                        for j in xrange(len(seq)):
                            base = seq[j].upper()
                            if base in colorBases:
                                vis2.extend([colorBases[base], 
                                             quads(start + j, -i,
                                                   start + j, -i-1,
                                                   start + j+1, -i-1,
                                                   start + j+1, -i)])
                    
                    end2 = start + len(seq)
                    
                    # draw text
                    if self.showBases and \
                       self.multiscale.atleast(mintextSize, 2, 
                                               view=view, size=size):
                        vis.append(text_scale(seq, 
                                              start, -i+1, 
                                              end2, -i-1, 
                                              "left", "bottom"))
                
                self.textGroup = win.replace_group(self.textGroup, 
                    group(group(*vis2), color(0,0,0), * vis))
            
        elif self.textShown:
            self.textShown = False
            self.textGroup = win.replace_group(self.textGroup, group())


    def onClickCallback(self):
        # TODO: should not rely on get_mouse_pos("world")
        win = self.get_window()
        x, y = win.get_mouse_pos('world')
        x -= self.pos[0]
        y = self.size[1] - (y - self.pos[1])
        self.onClick(x, y)


    def onClick(self, x, y):
        y = int(y)
    
        if 0 <= y < len(self.aln):
            print self.aln.keys()[y]
        
        

class TrackOverlay (Track):
    """Overlay several tracks on top of one another to create a single track"""
    
    def __init__(self, tracks, **options):
        Track.__init__(self, **options)
        
        self.tracks = tracks
        
        # process tracks
        maxheight = 0
        for track in tracks:
            track.view = self.view
            maxheight = max(maxheight, track.size[1])
        self.size = [0, maxheight]
    
    
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



def showAlign(* alns):
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
    browser.addTrack(RulerTrack(bottom=-height))
    for aln, col in zip(alns, colors):
        browser.addTrack(AlignTrack(aln, colorBases=col))
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
    pep_per_prop = util.histDict(util.mget(seqlib.AA_PROPERTY, AA))

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


def guessSeq(seq):
    """Guesses whether a sequence is 'dna' or 'pep'"""
    dna = "ACTG-N"
    
    chars = util.unique(seq.upper())
    
    for char in chars:
        if char not in dna:
            return "pep"
    return "dna"


def guessAlign(aln):
    """Guesses whether an alignment is 'dna' or 'pep'"""
    
    if "pep" in [guessSeq(seq) for seq in aln.itervalues()]:
        return "pep"
    else:
        return "dna"


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

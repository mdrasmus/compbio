import math
import summon
from summon.core import *
from summon.plot import ScatterPlot

from rasmus import util

from summon import VisObject


# TODO: make multiscale also track y-axis

class Multiscale (object):
    """Manage detecting when the zoom and scroll of the visualization is 
       sufficently different to justify a redraw"""

    def __init__(self, marginx=.5, marginy=.5, scalex=4, scaley=4):
        self.worldx1 = None
        self.worldx2 = None
        self.worldy1 = None
        self.worldy2 = None
        self.marginx = marginx
        self.marginy = marginy
        self.scalex = scalex
        self.scaley = scaley      
        self.win    = None
        
    
    def init(self, win=None, view=None):
        if win != None:
            self.win = win
        elif self.win == None:
            # backward compatibility
            self.win = summon.get_summon_window()
        
        if view == None:
            view = win.get_visible()
        self.worldx1, self.worldy1, self.worldx2, self.worldy2 = view
        
        # define outer bound with margins
        self.worldwidth = self.worldx2 - self.worldx1
        self.worldheight = self.worldy2 - self.worldy1
        marginx = self.worldwidth * self.marginx
        marginy = self.worldheight * self.marginx
        
        self.worldx1 -= marginx
        self.worldx2 += marginx
        self.worldy1 -= marginy
        self.worldy2 += marginy
        
    
    def sameScale(self, view=None):
        if view == None:
            view = self.win.get_visible()
        
        worldx1, worldy1, worldx2, worldy2 = view
        
        # test for scrolling
        if worldx1 < self.worldx1 or \
           worldx2 > self.worldx2 or \
           worldy1 < self.worldy1 or \
           worldy2 > self.worldy2:
            self.init(view=view)
            return False
        
        worldwidth = worldx2 - worldx1
        worldheight = worldy2 - worldy1
        
        if self.worldwidth == 0 or \
           self.worldheight == 0:
            return True
        
        # test for zooming
        if abs(math.log10(worldwidth / self.worldwidth)) > 1./self.scalex or \
           abs(math.log10(worldheight / self.worldheight)) > 1./self.scaley:
            self.init(view=view)
            return False
        
        return True


    def atleast(self, xminres, yminres, view=None, size=None):
        if view == None:
            view = self.win.get_visible()
        if size == None:
            size = self.win.get_size()
        
        worldx1, worldy1, worldx2, worldy2 = view
        screenwidth, screenheight = size
        worldwidth = worldx2 - worldx1
        worldheight = worldy2 - worldy1
        
        return (worldwidth == 0 or
                screenwidth / worldwidth > xminres) and \
               (worldheight == 0 or 
                screenheight / worldheight > yminres)




class Ruler (summon.VisObject):
    """ Ruler visualization object """
    
    def __init__(self, win, start, end, height=20, bottom=0, unitstr="", 
                 minicolor=color(.8,.8,.8), maincolor=color(0,0,0),
                 pos=[0.0, 0.0]):
        summon.VisObject.__init__(self)
         
        self.win = win
        self.start = start
        self.end = end
        self.height = height
        self.bottom = bottom
        self.pos = pos[:]
        self.minicolor = minicolor
        self.maincolor = maincolor
        self.multiscale = Multiscale(marginx=.5, marginy=.5, scalex=10, scaley=10)
    
    
    def draw(self):
        self.multiscale.init(self.get_window())
        
        g = group()
        self.gid = get_group_id(g)
        return g
        
    
    def update(self):
        if not self.multiscale.sameScale():
            g = drawRuler(self.pos, 
                          self.start, 
                          self.end, 
                          height=self.height, 
                          bottom=self.bottom,
                          minicolor=self.minicolor,
                          maincolor=self.maincolor,
                          win=self.win)
            self.gid = self.win.replace_group(self.gid, g)


    def drawRuler(self, pos, start, end, height=20, bottom=0, unit=None, 
                  unitstr="", 
                  minicolor=color(.8,.8,.8), 
                  maincolor=color(0,0,0)):
        win = self.win
        
        worldx1, worldy1, worldx2, worldy2 = win.get_visible()
        screenwidth, screenheight = win.get_size()

        worldwidth = worldx2 - worldx1
        worldx1 -= worldwidth / 2.0
        worldx2 += worldwidth / 2.0

        # find appropriate unit if one is not given
        unit = getRulerAutoSize(screenwidth, worldwidth)
        order = int(math.log10(unit))
        unit2, unitstr = getUnitSuffix(unit)


        x, y = pos
        vis = []

        # make mini hashes
        if unit >= 10:
            vis.append(minicolor)
            i = unit * (max(start, worldx1 - x + start) // unit)
            while x + i - start <= worldx2 and i < end:
                if i >= start:
                    vis.append(lines(x + i - start, y+bottom, x + i - start, y+height))
                i += unit // 10


        # make main hashes
        vis.append(maincolor)
        i = unit * (max(start, worldx1 - x + start) // unit)
        while x + i - start <= worldx2 and i < end:
            if i >= start:
                vis.append(lines(x + i - start, y, x + i - start, y + height))
                vis.append(text(str(int(i//unit2)) + unitstr, 
                                x + i - start, y, x + i -start - unit, y + height, "middle", "right"))
            i += unit

        # base line
        vis.append(lines(color(0,0,0), x, y, x + end - start, y))

        return group(* vis)


def getRulerAutoSize(screenwidth, worldwidth):
    """get most appropriate unit for zoom level"""
    unit = 1
    order = 0
    
    if worldwidth == 0.0:
        return unit
    
    while True:
        # find pixels per unit
        pixelsize = screenwidth / (worldwidth / float(unit))

        if pixelsize < 50 and order < 20:
            unit *= 10
            order += 1
        else:
            break
    
    return unit


def getUnitSuffix(unit):
    """get the sufffix for a unit"""
    
    order = int(math.log10(max(unit, 1)))
    if order < 3:
        unitstr = ""
        unit2 = 1
    elif 3 <= order < 6:
        unitstr = "K"
        unit2 = 1000
    elif 6 <= order < 9:
        unitstr = "M"
        unit2 = 1e6
    elif 9 <= order < 12:
        unitstr = "G"
        unit2 = 1e9
    elif 12 <= order < 15:
        unitstr = "T"
        unit2 = 1e12
    elif 15 <= order:
        unitstr = "e" + str(order)
        unit2 = unit
    
    return unit2, unitstr
    

def drawRuler(pos, start, end, height=20, bottom=0, unit=None, unitstr="", 
              minicolor=color(.8,.8,.8), maincolor=color(0,0,0),
              win=None):
    if win == None:
        # backwards compat.
        win = summon.get_summon_window()
    
    worldx1, worldy1, worldx2, worldy2 = win.get_visible()
    screenwidth, screenheight = win.get_size()
    
    worldwidth = worldx2 - worldx1
    worldx1 -= worldwidth / 2.0
    worldx2 += worldwidth / 2.0
    
    # find appropriate unit if one is not given
    unit = getRulerAutoSize(screenwidth, worldwidth)
    order = int(math.log10(unit))
    if order < 3:
        unitstr = ""
        unit2 = 1
    elif 3 <= order < 6:
        unitstr = "K"
        unit2 = 1000
    elif 6 <= order < 9:
        unitstr = "M"
        unit2 = 1e6
    elif 9 <= order < 12:
        unitstr = "G"
        unit2 = 1e9
    elif 12 <= order < 15:
        unitstr = "T"
        unit2 = 1e12
    elif 15 <= order:
        unitstr = "e" + str(order)
        unit2 = unit
    
    
    x, y = pos
    vis = []
    
    # make mini hashes
    if unit >= 10:
        vis.append(minicolor)
        i = unit * (max(start, worldx1 - x + start) // unit)
        while x + i - start <= worldx2 and i < end:
            if i >= start:
                vis.append(lines(x + i - start, y+bottom, x + i - start, y+height))
            i += unit // 10

    
    # make main hashes
    vis.append(maincolor)
    i = unit * (max(start, worldx1 - x + start) // unit)
    while x + i - start <= worldx2 and i < end:
        if i >= start:
            vis.append(lines(x + i - start, y, x + i - start, y + height))
            vis.append(text(str(int(i//unit2)) + unitstr, 
                            x + i - start, y, x + i -start - unit, y + height, "middle", "right"))
        i += unit
    
    # base line
    vis.append(lines(color(0,0,0), x, y, x + end - start, y))
    
    return group(* vis)

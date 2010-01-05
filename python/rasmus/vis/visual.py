import math
import summon
from summon.core import *
from summon.plot import ScatterPlot

from rasmus import util

from summon import VisObject
from summon.multiscale import Multiscale



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
        if not self.multiscale.same_view():
            g = draw_ruler(self.win,
                           self.pos, 
                           self.start, 
                           self.end, 
                           height=self.height, 
                           bottom=self.bottom,
                           minicolor=self.minicolor,
                           maincolor=self.maincolor,
                           win=self.win)
            self.gid = self.win.replace_group(self.gid, g)


    def draw_ruler(self, pos, start, end, height=20, bottom=0, unit=None, 
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
    drawRuler = draw_ruler


def get_ruler_auto_size(screenwidth, worldwidth):
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
getRulerAutoSize = get_ruler_auto_size

def get_unit_suffix(unit):
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
getUnitSuffix = get_unit_suffix
    

def draw_ruler(win, pos, start, end, height=20, bottom=0, unit=None,
               unitstr="", 
               minicolor=color(.8,.8,.8), maincolor=color(0,0,0)):

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
drawRuler = draw_ruler


import math
from summon import *
import summonlib

class VisObject:
    def __init__(self):
        self.vis = group()
    
    def update(self):
        pass
    
    def setVisible(self, visible=True):
        if visible:    
            if not summonlib.is_update_func(self.update):
                summonlib.add_update_func(self.update)
        else:
            if summonlib.is_update_func(self.update):
                summonlib.remove_update_func(self.update)
        
    

class Ruler (VisObject):
    def __init__(self, gid, start, end, height=20, bottom=0, unitstr="", 
                 minicolor=color(.8,.8,.8), maincolor = color(0,0,0)):
        VisObject.__init__(self)
        
        self.gid = insert_group(gid, drawRuler(start, end, height=height, 
                                               bottom=bottom,
                                               unit=None, 
                                               unitstr=unitstr, 
                                               minicolor=minicolor, 
                                               maincolor=maincolor))

        worldx1, worldy1, worldx2, worldy2 = get_visible()
        screenwidth, screenheight = get_window_size()
        worldwidth = worldx2 - worldx1
        
        self.start = start
        self.end = end
        self.height = height
        self.bottom = bottom
        self.unit = getRulerAutoSize(screenwidth, worldwidth)
        self.worldx1 = worldx1 - worldwidth / 2.0
        self.worldx2 = worldx2 + worldwidth / 2.0
        
        self.setVisible()
    
    
    def update(self):
        worldx1, worldy1, worldx2, worldy2 = get_visible()
        screenwidth, screenheight = get_window_size()    
        unit = getRulerAutoSize(screenwidth, worldx2 - worldx1)
        
        if unit != self.unit or \
           worldx1 < self.worldx1 or \
           worldx2 > self.worldx2:
            g = drawRuler(self.start, self.end, height=self.height, 
                          bottom=self.bottom)
            self.gid = replace_group(self.gid, g)
            
        self.unit = unit
        self.worldx1 = worldx1
        self.worldx2 = worldx2


def getRulerAutoSize(screenwidth, worldwidth):
    unit = 1
    order = 0
    
    if worldwidth == 0.0:
        return unit
    
    while True:
        # find pixels per unit
        pixelsize = screenwidth / (worldwidth / float(unit))

        if pixelsize < 50:
            unit *= 10
            order += 1
        else:
            break
    
    return unit

    

def drawRuler(start, end, height=20, bottom=0, unit=None, unitstr="", 
              minicolor=color(.8,.8,.8), maincolor = color(0,0,0)):
    
    worldx1, worldy1, worldx2, worldy2 = get_visible()
    screenwidth, screenheight = get_window_size()
    
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

    
    vis = []
    
    # make mini hashes
    if unit >= 10:
        vis.append(minicolor)
        i = unit * (max(0, worldx1) // unit)
        while i <= worldx2 and i < end:
            vis.append(lines(i, bottom, i, height))
            i += unit // 10

    
    # make main hashes
    vis.append(maincolor)
    i = unit * (max(0, worldx1) // unit)
    while i <= worldx2 and i < end:
        vis.append(lines(i, 0, i, height))
        vis.append(text(str(int(i//unit2)) + unitstr, 
                        i, 0, i+unit, height, "middle", "left"))
        i += unit
    
    # base line
    vis.append(lines(color(0,0,0), start, 0, end, 0))
    
    return group(* vis)

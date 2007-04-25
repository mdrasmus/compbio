import math
from summon.core import *
import summon
from rasmus import util



class VisObject (object):
    """Base class of visualization objects"""
    def __init__(self):
        self.win = None
    
    def __del__(self):
        if summon != None:
            self.setVisible(False)
    
    def update(self):
        pass
    
    def show(self):
        pass
    
    def enableUpdating(self, visible=True):
        if visible:    
            if not summon.is_update_func(self.update):
                assert self.win != None, "must set window"
                summon.add_update_func(self.update, self.win)
        else:
            if summon.is_update_func(self.update):
                summon.remove_update_func(self.update)

    setVisible = enableUpdating
    

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
        
    
    def init(self, view=None):
        if view == None:
            view = get_visible()
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
            view = get_visible()
        
        worldx1, worldy1, worldx2, worldy2 = view
        
        # test for scrolling
        if worldx1 < self.worldx1 or \
           worldx2 > self.worldx2 or \
           worldy1 < self.worldy1 or \
           worldy2 > self.worldy2:
            self.init(view)
            return False
        
        worldwidth = worldx2 - worldx1
        worldheight = worldy2 - worldy1
        
        if self.worldwidth == 0 or \
           self.worldheight == 0:
            return True
        
        # test for zooming
        if abs(math.log10(worldwidth / self.worldwidth)) > 1./self.scalex or \
           abs(math.log10(worldheight / self.worldheight)) > 1./self.scaley:
            self.init(view)
            return False
        
        return True


    def atleast(self, xminres, yminres, view=None):
        if view == None:
            view = get_visible()
        
        worldx1, worldy1, worldx2, worldy2 = view
        screenwidth, screenheight = get_window_size()
        worldwidth = worldx2 - worldx1
        worldheight = worldy2 - worldy1
        
        return (worldwidth == 0 or
                screenwidth / worldwidth > xminres) and \
               (worldheight == 0 or 
                screenheight / worldheight > yminres)



class ScatterPlot (object):
    def __init__(self, scale=[1.0, 1.0], onClick=None):
        self.data = []
        self.scale = scale
        self.selgroup = None
        
        if onClick != None:
            self.onClick = onClick
        
    
    def plot(self, x, y, names, col=[0,0,0], style="points", size=1.0):
        self.data.append(util.Bundle(x=x, y=y, 
                                     names=names, 
                                     color=col,
                                     style=style,
                                     size=size))
    
    
    def show(self):
        self.win = summon.Window()
        self.win.set_bgcolor(1, 1, 1)
        
        self.win.add_group(group(scale(self.scale[0], self.scale[1],
                                       self.draw())))
        
        self.win.reset_binding(input_key('d'), self.clearSelection)
        
        minx = util.INF
        maxx = -util.INF
        miny = util.INF
        maxy = -util.INF
            
        for dat in self.data:
            minx = min(minx, min(dat.x))
            maxx = max(maxx, max(dat.x))
            miny = min(miny, min(dat.y))
            maxy = max(maxy, max(dat.y))
        self.win.set_visible(minx * self.scale[0], miny * self.scale[1], 
                             maxx * self.scale[0], maxy * self.scale[1])
    
        
    def draw(self):
        vis = []
        
        minx = util.INF
        maxx = -util.INF
        miny = util.INF
        maxy = -util.INF        
        
        # draw data
        for dat in self.data:
            minx = min(minx, min(dat.x))
            maxx = max(maxx, max(dat.x))
            miny = min(miny, min(dat.y))
            maxy = max(maxy, max(dat.y))
        
        
            if dat.style == "points":
                vis2 = [color(* dat.color)]

                for x, y in zip(dat.x, dat.y):
                    vis2.extend([x, y])
                vis.append(points(* vis2))
            
            elif dat.style == "diamonds":
                vis2 = [color(* dat.color)]
                vis3 = []
                
                for i in xrange(len(dat.x)):
                    x = dat.x[i]
                    y = dat.y[i]
                    size = dat.size
                
                    vis2.extend([x + size, y,
                                 x       , y + size,
                                 x - size, y,
                                 x       , y - size])
                    vis3.extend([x, y + size, x, y - size])
                vis.append(quads(*vis2))
                vis.append(lines(*vis3))
        
        # add hotspot
        width = maxx * self.scale[0] - minx * self.scale[0]
        height = maxy * self.scale[1] - miny * self.scale[1]
        vis.append(hotspot("click", minx * self.scale[0] - width*.1, 
                                    miny * self.scale[1] - height*.1, 
                                    maxx * self.scale[1] + height*.1, 
                                    maxy * self.scale[1] + height*.1, 
                                    self.clickCallback))
    
        selgroup = group()
        self.selgroup = get_group_id(selgroup)
        vis.append(group(selgroup))
            
        return group(*vis)
    
    
    def clickCallback(self):
        x, y = self.win.get_mouse_pos('world')
        x /= float(self.scale[0])
        y /= float(self.scale[1])
        
        closest_dist = util.INF
        closest = None
        
        for dat in self.data:
            x2 = dat.x
            y2 = dat.y
        
            for i in xrange(len(dat.x)):
                dist = math.sqrt((x2[i] - x)**2 + (y2[i] - y)**2)
                
                if dist < closest_dist:
                    closest_dist = dist
                    closest = (dat, i)
        
        if closest != None:
            dat, i = closest
            x, y, size = dat.x[i], dat.y[i], dat.size
            self.onClick(x, y, dat.names[i])
            
            self.win.insert_group(self.selgroup, 
                                  self.drawSelect(x, y, size, col=[1, 0, 0]))

    def drawSelect(self, x, y, size, col):
        return group(color(* col), 
                     line_strip(x - size, y - size,
                                x - size, y + size,
                                x + size, y + size,
                                x + size, y - size,
                                x - size, y - size))

    def onClick(self, x, y, name):
        print name, x, y
        

    def clearSelection(self):
        self.selgroup = self.win.replace_group(self.selgroup, group())

    
    def select(self, func, col=[0, 0, 1]):
        vis = []
        sel = []
    
        for dat in self.data:
            x = dat.x
            y = dat.y
            names = dat.names
            size = dat.size
            
            for i in xrange(len(dat.x)):
                if func(x[i], y[i], names[i]):
                    sel.append([x[i], y[i], names[i]])
                    vis.append(self.drawSelect(x[i], y[i], size, col=col))
        
        if len(sel) > 0:
            self.win.insert_group(self.selgroup, group(*vis))
        
        return sel
            
            


# TODO: convert to use multiscale
class Ruler (VisObject):
    """ Ruler visualization object """
    
    def __init__(self, gid, start, end, height=20, bottom=0, unitstr="", 
                 minicolor=color(.8,.8,.8), maincolor=color(0,0,0),
                 pos=[0.0, 0.0]):
        VisObject.__init__(self)
        
        self.gid = gid 
        
        worldx1, worldy1, worldx2, worldy2 = get_visible()
        screenwidth, screenheight = get_window_size()
        worldwidth = worldx2 - worldx1
        
        self.start = start
        self.end = end
        self.height = height
        self.bottom = bottom
        self.pos = pos[:]
        self.minicolor = minicolor
        self.maincolor = maincolor        
        self.unit = getRulerAutoSize(screenwidth, worldwidth)
        self.worldx1 = worldx1 - worldwidth / 2.0
        self.worldx2 = worldx2 + worldwidth / 2.0
    
    
    def update(self):
        worldx1, worldy1, worldx2, worldy2 = get_visible()
        screenwidth, screenheight = get_window_size()    
        unit = getRulerAutoSize(screenwidth, worldx2 - worldx1)
        
        worldwidth = worldx2 - worldx1
        xmargin = worldwidth / 2.0
        
        if unit != self.unit or \
           worldx1 < self.worldx1 or \
           worldx2 > self.worldx2:
            g = drawRuler(self.pos, 
                          self.start, 
                          self.end, 
                          height=self.height, 
                          bottom=self.bottom,
                          minicolor=self.minicolor,
                          maincolor=self.maincolor)
            self.gid = replace_group(self.gid, g)
            
            self.worldx1 = worldx1 - xmargin
            self.worldx2 = worldx2 + xmargin            
            self.unit = unit
        
        


def getRulerAutoSize(screenwidth, worldwidth):
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


def drawRuler(pos, start, end, height=20, bottom=0, unit=None, unitstr="", 
              minicolor=color(.8,.8,.8), maincolor=color(0,0,0)):
    
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

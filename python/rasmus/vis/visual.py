from summon import *

def drawRuler(start, end, unit=None, unitstr="", 
              minicolor=color(.8,.8,.8), maincolor = color(0,0,0)):
    
    worldx1, worldy1, worldx2, worldy2 = map(int, get_visible())
    screenwidth, screenheight = get_window_size()
    
    # find appropriate unit if one is not given
    if unit == None:
        unit = 1
        order = 0
        
        while True:
            # find pixels per unit
            pixelsize = screenwidth / ((worldx2 - worldx1) / float(unit))
            
            if pixelsize < 20:
                unit *= 10
                order += 1
            else:
                break
        
        if order > 0:
            unitstr = "e" + str(order)
    
    vis = []
    
    

    # make mini hashes
    if unit >= 10:
        vis.append(minicolor)
        i = 0
        while i < start: i += unit/10
        while i < end:
            if worldx1 <= i <= worldx2:
                vis.append(lines(i, 0, i, 1))
            i += unit/10

    
    # make main hashes
    vis.append(maincolor)
    i = 0
    while i < start: i += unit
    while i < end:
        if worldx1 <= i <= worldx2:
            vis.append(lines(i, 0, i, 1))
            vis.append(text(str(i/unit) + unitstr, 
                            i, 0, i+unit, 1, "middle", "left"))
        i += unit
    
    # base line
    vis.append(lines(color(0,0,0), start, 0, end, 0))
    
    return group(* vis)

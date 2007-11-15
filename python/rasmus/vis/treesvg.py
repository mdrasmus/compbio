
from rasmus import svg
from rasmus import util
from rasmus import treelib
from rasmus.bio import phylo

import sys
import math
import os



def drawTree(tree, labels={}, xscale=100, yscale=20, canvas=None,
             labelOffset=None, fontSize=10, labelSize=None,
             minlen=1, maxlen=util.INF, filename=sys.stdout,
             rmargin=100, lmargin=10, tmargin=0, bmargin=None,
             colormap=None,
             stree=None,
             gene2species=None,
             lossColor=(0, 0, 1),
             dupColor=(1, 0, 0),
             eventSize=4,
             legendScale=False, autoclose=None):
    
    # set defaults
    fontRatio = 8. / 11.
    
    if labelSize == None:
        labelSize = .7 * fontSize
    
    if labelOffset == None:
        labelOffset = -1
    
    if bmargin == None:
        bmargin = yscale
    
    if sum(x.dist for x in tree.nodes.values()) == 0:
        legendScale = False
        minlen = yscale * 2
    
    if colormap == None:
        for node in tree:
            node.color = (0, 0, 0)
    else:
        colormap(tree)
    
    if stree and gene2species:
        recon = phylo.reconcile(tree, stree, gene2species)
        events = phylo.labelEvents(tree, recon)
        losses = phylo.findLoss(tree, stree, recon)
    else:
        events = None
        losses = None
    
    # layout tree
    coords = treelib.layoutTree(tree, xscale, yscale, minlen, maxlen)
    
    xcoords, ycoords = zip(* coords.values())
    maxwidth = max(xcoords)
    maxheight = max(ycoords) + labelOffset
    
    
    # initialize canvas
    if canvas == None:
        canvas = svg.Svg(util.openStream(filename, "w"))
        width = int(rmargin + maxwidth + lmargin)
        height = int(tmargin + maxheight + bmargin)
        
        canvas.beginSvg(width, height)
        
        if autoclose == None:
            autoclose = True
    else:
        if autoclose == None:
            autoclose = False
    
    
    # draw tree
    def walk(node):
        x, y = coords[node]
        if node.parent:
            parentx = coords[node.parent][0]
        else:
            parentx = 0
        
        # draw branch
        canvas.line(parentx, y, x, y, color=node.color)
        if node.name in labels:
            branchlen = x - parentx
            lines = str(labels[node.name]).split("\n")
            labelwidth = max(map(len, lines))
            labellen = min(labelwidth * fontRatio * fontSize, 
                           max(int(branchlen-1), 0))
            
            for i, line in enumerate(lines):
                canvas.text(line,
                            parentx + (branchlen - labellen)/2., 
                            y + labelOffset 
                            +(-len(lines)+1+i)*(labelSize+1),
                            labelSize)
        
        if node.isLeaf():
            canvas.text(str(node.name), 
                        x + fontSize, y+fontSize/2., fontSize,
                        fillColor=node.color)
        else:
            top = coords[node.children[0]][1]
            bot = coords[node.children[-1]][1]
            
            # draw children
            canvas.line(x, top, x, bot)
            
            for child in node.children:
                walk(child)
    
    canvas.beginTransform(("translate", lmargin, tmargin))
    walk(tree.root)
        
    if stree and gene2species:
        drawEvents(canvas, tree, coords, events, losses,
                   lossColor=lossColor,
                   dupColor=dupColor,
                   size=eventSize)
    canvas.endTransform()
    
    # draw legend
    if legendScale:
        if legendScale == True:
            # automatically choose a scale
            length = maxwidth / float(xscale)
            order = math.floor(math.log10(length))
            length = 10 ** order
    
        drawScale(lmargin, tmargin + maxheight + bmargin - fontSize, 
                  length, xscale, fontSize, canvas=canvas)
    
    if autoclose:
        canvas.endSvg()
    
    return canvas


def drawEvents(canvas, tree, coords, events, losses,
               lossColor=(0, 0, 1),
               dupColor=(1, 0, 0),
               size=4):

    # draw duplications
    for node in tree:
        x, y = coords[node]
        if events[node] == "dup":
            canvas.rect(x - size/2.0, y - size/2.0,
                        size, size,  fillColor=dupColor, strokeColor=(0,0,0,0))

    # draw losses
    losses_per_branch = util.histDict([node for node, schild in losses])

    for node, nlosses in losses_per_branch.iteritems():
        if node.parent == None:
            continue

        x1 = coords[node.parent][0]
        x2, y1 = coords[node]
        step = (x2 - x1) / float(nlosses + 1)

        for x in util.frange(x1 + step, x2-(step/2.0), step):
            canvas.line(x, y1 - size/2.0, x, y1 + size/2.0, color=lossColor)
            

def drawScale(x, y, length, xscale, fontSize, canvas=None):
    assert canvas != None
    
    canvas.line(x, y, x + length * xscale, y)
    canvas.line(x, y+1, x, y-1)
    canvas.line(x + length * xscale, y+1, x + length * xscale, y-1)
    canvas.text("%.3f" % length, x, y-1, fontSize)
    


def drawDistRuler(names, dists, scale=500,
                  padding=10, textsize=12, notchsize=2,
                  labelpadding=5, distsize=9,
                  filename=sys.stdout):
    """Produce a ruler of pairwise distances"""

    nameswidth = textsize * max(map(len, names))
    
    out = svg.Svg(util.openStream(filename, "w"))
    out.beginSvg(scale * max(dists) + 2*padding,
                 2*padding+nameswidth + 5*distsize)

    
    
    out.beginTransform(("translate", padding, nameswidth+padding))

    # draw ruler
    out.line(0, 0, scale*max(dists), 0)

    for name, dist in zip(names, dists):
        x = scale*dist
        out.text(name, x + textsize/2.0, - labelpadding, textsize, angle=-90)
        out.line(x, notchsize, x, - notchsize)
        out.text("%.3f" % dist, x + textsize/2.0, labelpadding + distsize*3.5,
                 distsize, angle=-90)
    
    out.endTransform()
    out.endSvg()


def getPairwiseDists(tree, mainsp):
    """get pairwise distances between the main taxon and the rest of the taxa"""
    lst = []

    for sp in tree.leafNames():
        lst.append((sp, findDist(tree, mainsp, sp)))

    lst.sort(key=lambda x: x[1])
    names, dists = zip(* lst)
    return names, dists

    
#=============================================================================
# additional wrapper functions

def drawTreeLens(tree, *args, **kargs):
    labels = {}
    for node in tree.nodes.values():
        if node == tree.root:
            continue
        labels[node.name] = "%.3f" % node.dist
    
    drawTree(tree, labels, *args, **kargs)
    

def showTree(tree, *args, **kargs):
    if not os.fork():
        kargs.setdefault('filename', os.popen("display", "w"))
        kargs.setdefault('legendScale', True)
        drawTree(tree, *args, **kargs)
        sys.exit(0)

    

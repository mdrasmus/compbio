
from rasmus import svg
from rasmus import util

import sys
import math
import os


def layoutTree(tree, xscale, yscale, minlen, maxlen):
    sizes = {}
    nodept = {}
    maxwidth = [0]
    maxheight = len(tree.leaves()) * yscale
    
    def walk(node, y):
        dist = min(maxlen, node.dist * xscale)
        dist = max(minlen, dist)
        y = y + dist
        
        if node.isLeaf():
            sizes[node] = 1
            nodept[node] = yscale - 1
            maxwidth[0] = max(maxwidth[0], y)
        else:
            sizes[node] = 0
        for child in node.children:
            sizes[node] += walk(child, y)
        if not node.isLeaf():
            top = nodept[node.children[0]]
            bot = (sizes[node] - sizes[node.children[-1]])*yscale + \
                  nodept[node.children[-1]]
            nodept[node] = (top + bot) / 2
        return sizes[node]
    walk(tree.root, 0)
    maxwidth = maxwidth[0]
    
    return sizes, nodept, maxwidth, maxheight


def drawTree(tree, labels={}, xscale=100, yscale=20, canvas=None,
             labelOffset=None, fontSize=10, labelSize=None,
             minlen=1, maxlen=10000, filename=sys.stdout,
             rmargin=100, lmargin=10, tmargin=0, bmargin=None,
             legendScale=False, autoclose=None):

    
    if labelSize == None:
        labelSize = .7 * fontSize
    
    if labelOffset == None:
        labelOffset = labelSize - 1
    
    if bmargin == None:
        bmargin = yscale
    
    if sum(x.dist for x in tree.nodes.values()) == 0:
        legendScale = False
        minlen = yscale * 2
    
    
    # determine node sizes
    sizes, nodept, maxwidth, maxheight = layoutTree(tree, xscale, yscale, minlen, maxlen)
    maxheight += labelOffset
    
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
    
    
    def walk(node, x, y):
        # calc coords
        xchildren = x+min(max(node.dist*xscale, minlen), maxlen)
        
        # draw branch
        canvas.line(x, y+nodept[node], xchildren, y+nodept[node])
        if node.name in labels:
            branchlen = xchildren - x
            lines = str(labels[node.name]).split("\n")
            labelwidth = max(map(len, lines))
            
            labellen = min(labelwidth, 
                           max(int(branchlen-1),0))
            canvas.text(labels[node.name],
                        x + (branchlen - labellen)/2., 
                        y+nodept[node] + labelOffset - labelSize,
                        labelSize)
        
        if node.isLeaf():
            canvas.text(str(node.name), xchildren + fontSize, y+yscale+fontSize/2., fontSize)
        else:
            top = y + nodept[node.children[0]]
            bot = y + (sizes[node]-sizes[node.children[-1]]) * yscale + \
                      nodept[node.children[-1]]
        
            # draw children
            canvas.line(xchildren, top, xchildren, bot)
            
            ychild = y
            for child in node.children:
                walk(child, xchildren, ychild)
                ychild += sizes[child] * yscale

            
            #canvas.set(xchildren, y+nodept[node], '+')
            #canvas.set(xchildren, top, '/')
            #canvas.set(xchildren, bot, '\\')
        #canvas.set(x, y+nodept[node], '+')
    walk(tree.root, lmargin, tmargin)
    
    
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


def drawScale(x, y, length, xscale, fontSize, canvas=None):
    assert canvas != None
    
    canvas.line(x, y, x + length * xscale, y)
    canvas.line(x, y+1, x, y-1)
    canvas.line(x + length * xscale, y+1, x + length * xscale, y-1)
    canvas.text("%.3f" % length, x, y-1, fontSize)
    
    


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

    

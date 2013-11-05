
from rasmus import svg
from rasmus import util
from rasmus import treelib
from compbio import phylo

import sys
import math
import os





def draw_tree(tree, labels={}, xscale=100, yscale=20, canvas=None,
              leafPadding=10, leafFunc=lambda x: str(x.name),
              labelOffset=None, fontSize=10, labelSize=None,
              minlen=1, maxlen=util.INF, filename=sys.stdout,
              rmargin=150, lmargin=10, tmargin=0, bmargin=None,
              colormap=None,
              stree=None,
              layout=None,
              gene2species=None,
              lossColor=(0, 0, 1),
              dupColor=(1, 0, 0),
              eventSize=4,
              legendScale=False, autoclose=None,
              extendRoot=True, labelLeaves=True, drawHoriz=True, nodeSize=0):
    
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
        minlen = xscale
    
    if colormap == None:
        for node in tree:
            node.color = (0, 0, 0)
    else:
        colormap(tree)
    
    if stree and gene2species:
        recon = phylo.reconcile(tree, stree, gene2species)
        events = phylo.label_events(tree, recon)
        losses = phylo.find_loss(tree, stree, recon)
    else:
        events = None
        losses = None

    if len(labels) > 0 or (stree and gene2species):
        drawHoriz = True
    
    # layout tree
    if layout is None:
        coords = treelib.layout_tree(tree, xscale, yscale, minlen, maxlen)
    else:
        coords = layout
    
    xcoords, ycoords = zip(* coords.values())
    maxwidth = max(xcoords)
    maxheight = max(ycoords) + labelOffset
    
    
    # initialize canvas
    if canvas == None:
        canvas = svg.Svg(util.open_stream(filename, "w"))
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
            parentx, parenty = coords[node.parent]
        else:
            if extendRoot:
                parentx, parenty = 0, y
            else:
                parentx, parenty = x, y     # e.g. no branch
        
        # draw branch
        if drawHoriz:
            canvas.line(parentx, y, x, y, color=node.color)
        else:
            canvas.line(parentx, parenty, x, y, color=node.color)

        # draw branch labels
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

        # draw nodes
        if nodeSize > 0:
            canvas.circle(x, y, nodeSize, strokeColor=svg.null, fillColor=node.color)

        # draw leaf labels or recur
        if node.is_leaf():
            if labelLeaves:
                canvas.text(leafFunc(node), 
                            x + leafPadding, y+fontSize/2., fontSize,
                            fillColor=node.color)
        else:
            if drawHoriz:
                # draw vertical part of branch
                top = coords[node.children[0]][1]
                bot = coords[node.children[-1]][1]
                canvas.line(x, top, x, bot, color=node.color)
                
            # draw children
            for child in node.children:
                walk(child)
    
    canvas.beginTransform(("translate", lmargin, tmargin))
    walk(tree.root)
        
    if stree and gene2species:
        draw_events(canvas, tree, coords, events, losses,
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


def draw_events(canvas, tree, coords, events, losses,
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
    losses_per_branch = util.hist_dict([node for node, schild in losses])

    for node, nlosses in losses_per_branch.iteritems():
        if node.parent == None:
            continue

        x1 = coords[node.parent][0]
        x2, y1 = coords[node]
        step = (x2 - x1) / float(nlosses + 1)

        for x in util.frange(x1 + step, x2-(step/2.0), step):
            canvas.line(x, y1 - size, x, y1 + size, color=lossColor)

            

def drawScale(x, y, length, xscale, fontSize, canvas=None):
    assert canvas != None

    color = (0,0,0)
    
    canvas.line(x, y, x + length * xscale, y, color=color)
    canvas.line(x, y+1, x, y-1, color=color)
    canvas.line(x + length * xscale, y+1, x + length * xscale, y-1, color=color)
    canvas.text("%.3f" % length, x, y-1, fontSize)
    


def drawDistRuler(names, dists, scale=500,
                  padding=10, textsize=12, notchsize=2,
                  labelpadding=5, distsize=9,
                  filename=sys.stdout):
    """Produce a ruler of pairwise distances"""

    nameswidth = textsize * max(map(len, names))
    
    out = svg.Svg(util.open_stream(filename, "w"))
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


def get_pairwise_dists(tree, mainsp):
    """get pairwise distances between the main taxon and the rest of the taxa"""
    lst = []

    for sp in tree.leaf_names():
        lst.append((sp, treelib.find_dist(tree, mainsp, sp)))

    lst.sort(key=lambda x: x[1])
    names, dists = zip(* lst)
    return names, dists

    
#=============================================================================
# additional wrapper functions

def draw_events_tree(tree, stree, gene2species, **options):
    phylo.init_dup_loss_tree(stree)
    phylo.count_dup_loss_tree(tree, stree, gene2species)
    phylo.count_ancestral_genes(stree)
    
    labels = options.get("labels", {})
    
    for node in stree:
        labels[node.name] = ""

        if node.data['dup'] > 0:
            labels[node.name] += " +%d" % node.data['dup']
        else:
            labels[node.name] += "   "

        if node.data['loss'] > 0:
            labels[node.name] += " -%d" %  node.data['loss']
        else:
            labels[node.name] += "   "
        
        if labels[node.name] != "":
            labels[node.name] += ": "
        labels[node.name] += "%d" % node.data['genes']
    
    draw_tree(stree, labels=labels, **options)


def draw_tree_lens(tree, *args, **kargs):
    labels = {}
    for node in tree.nodes.values():
        if node == tree.root:
            continue
        labels[node.name] = "%.3f" % node.dist
    
    draw_tree(tree, labels, *args, **kargs)


def show_tree(tree, *args, **kargs):
    if not os.fork():
        kargs.setdefault('filename', os.popen("display", "w"))
        kargs.setdefault('legendScale', True)
        draw_tree(tree, *args, **kargs)
        sys.exit(0)

    

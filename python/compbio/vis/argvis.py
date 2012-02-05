"""

    Visualization of ARGs (Ancestral Recombination Graphs)

"""

# python imports
from itertools import chain

# rasmus imports
from rasmus import util, treelib, stats

# compbio imports
from compbio import arglib

# summon imports
import summon
from summon import sumtree
from summon.shapes import box
from summon.core import *




#=============================================================================
# visualization

def layout_arg(arg, leaves=None, yfunc=lambda x: x):

    layout = {}

    if leaves is None:
        leaves = sorted((i for i in arg.leaves()), key=lambda x: x.name)

    # layout leaves
    leafx = util.list2lookup(leaves)
    
    for node in arg.postorder():
        if node.is_leaf():
            layout[node] = [leafx[node], yfunc(node.age)]
        else:
            layout[node] = [
                stats.mean(layout[child][0] for child in node.children),
                yfunc(node.age)]

    return layout


def show_arg(arg, layout=None, leaves=None, mut=None,
             recomb_width=.4, recomb_width_expand=0):

    win = summon.Window()
    
    if layout is None:
        layout = layout_arg(arg, leaves)

    def branch_hotspot(node, parent, x, y, y2):
        def func():
            print node.name, parent.name
        return hotspot("click", x-.5, y, x+.5, y2, func)

    for node in layout:
        recomb_width2 = recomb_width + node.age * recomb_width_expand
        
        if not node.is_leaf():
            x, y = layout[node]
            for i, child in enumerate(node.children):
                x2, y2 = layout[child]
                step = 0.0
                
                if child.event == "recomb":
                    if (len(child.parents) == 2 and
                        child.parents[0] == child.parents[1]):
                        step = recomb_width2 * [-1, 1][i]
                    else:
                        step = recomb_width2 * [-1, 1][
                            child.parents.index(node)]
                    win.add_group(line_strip(x, y,
                                             x2+step, y,
                                             x2+step, y2,
                                             x2, y2))                    
                else:
                    win.add_group(line_strip(x, y, x2, y, x2, y2))

                win.add_group(
                    branch_hotspot(child, node, x2+step, y, y2))

            # draw mutation
            if node.event == "recomb":
                win.add_group(zoom_clamp(
                    color(1, 0, 0),
                    box(x-.5, y-.5, x+.5, y+.5, fill=True),
                    color(1,1,1),
                    origin=(x, y),
                    minx=4.0, miny=4.0, maxx=20.0, maxy=20.0,
                    link=True))


    # draw mutations
    if mut:
        for node, parent, pos, t in mut:
            x, y = layout[parent]
            x2, y2 = layout[node]
            recomb_width2 = recomb_width + node.age * recomb_width_expand
            
            if node.event == "recomb":
                if (len(node.parents) == 2 and
                    node.parents[0] == node.parents[1]):
                    step = recomb_width2 * [-1, 1][i]
                else:
                    step = recomb_width2 * [-1, 1][node.parents.index(parent)]
            else:
                step = 0.0

            mx = x2+step
            my = t
            
            win.add_group(zoom_clamp(
                    color(0, 0, 1),
                    box(mx-.5, my-.5, mx+.5, my+.5, fill=True),
                    color(1,1,1),
                    origin=(mx, my),
                    minx=4.0, miny=4.0, maxx=20.0, maxy=20.0,
                    link=True))



    return win
    


def show_marginal_trees2(arg, mut=None):

    def minlog(x, low=10):
        return log(max(x, low))
    #ymap = lambda x: x
    ymap = minlog

    win = summon.Window()
    x = 0
    step = 2
    treewidth = len(list(arg.leaves())) + step

    def trans_camera(win, x, y):
        v = win.get_visible()
        win.set_visible(v[0]+x, v[1]+y, v[2]+x, v[3]+y, "exact")

    win.set_binding(input_key("]"), lambda : trans_camera(win, treewidth, 0))
    win.set_binding(input_key("["), lambda : trans_camera(win, -treewidth, 0))

    blocks = arglib.iter_recomb_blocks(arg)

    for tree, block in izip(arglib.iter_marginal_trees(arg), blocks):
        pos = block[0]
        print pos
        
        leaves = sorted((x for x in tree.leaves()), key=lambda x: x.name)
        layout = layout_arg(tree, leaves, yfunc=ymap)
        win.add_group(translate(x, 0, color(1,1,1),
                                draw_tree(tree, layout)))

        # mark responsible recomb node
        for node in tree:
            if pos != 0.0 and node.pos == pos:
                nx, ny = layout[node]
                win.add_group(draw_mark(x + nx, ny))

        # draw mut
        if mut:
            for node, parent, mpos, t in mut:
                if (node.name in tree and node.name != tree.root.name and
                    block[0] < mpos < block[1]):
                    nx, ny = layout[tree[node.name]]
                    win.add_group(draw_mark(x + nx, ymap(t), col=(0,0,1)))
                if node.name in tree and tree[node.name].parents:
                    nx, ny = layout[tree[node.name]]
                    py = layout[tree[node.name].parents[0]][1]
                    start = arg[node.name].data["ancestral"][0][0]
                    win.add_group(lines(color(0,1,0), 
                                        x+nx, ny, x+nx, py,
                                        color(1,1,1)))
            
                
        x += treewidth

    win.set_visible(* win.get_root().get_bounding() + ("exact",))

    return win


def show_tree_track(tree_track, mut=None, show_labels=False,
                    use_blocks=False):
    """
    tree_track = [((start, end), tree), ...]
    """

    def draw_labels(tree, layout):
        return group(*
                [text_clip(leaf.name, layout[leaf][0], layout[leaf][1],
                          1, layout[leaf][1] + 1e4, 4, 20, "middle", "left")
                 for leaf in tree.leaves()])
    


    def minlog(x, low=10):
        return log(max(x, low))
    #ymap = lambda x: x
    ymap = minlog

    tree_track = iter(tree_track)
    block, tree = tree_track.next()

    win = summon.Window()
    x = 0
    step = 2
    treewidth = len(list(tree.leaves())) + step

    def trans_camera(win, x, y):
        v = win.get_visible()
        win.set_visible(v[0]+x, v[1]+y, v[2]+x, v[3]+y, "exact")

    win.set_binding(input_key("]"), lambda : trans_camera(win, treewidth, 0))
    win.set_binding(input_key("["), lambda : trans_camera(win, -treewidth, 0))

    for block, tree in chain([(block, tree)], tree_track):
        pos = block[0]
        print pos

        #if use_blocks:
        #    treewidth = block[1] - block[0]
        
        layout = treelib.layout_tree(tree, xscale=1, yscale=1)
        treelib.layout_tree_vertical(layout, leaves=0)
        win.add_group(
            translate(x, 0, color(1,1,1),
                      sumtree.draw_tree(tree, layout, 
                                        vertical=True),
                      (draw_labels(tree, layout) if show_labels else group()),
                      text_clip(
                    "%d-%d" % (block[0], block[1]),
                    treewidth*.05, 0, 
                    treewidth*.95, -max(l[1] for l in layout.values()),
                    4, 20, 
                    "center", "top")))
        
        # TODO: update with tree width
        # draw mut
        if mut:
            for node, parent, mpos, t in mut:
                if (node.name in tree and node.name != tree.root.name and
                    block[0] < mpos < block[1]):
                    nx, ny = layout[tree[node.name]]
                    win.add_group(draw_mark(x + nx, ymap(t), col=(0,0,1)))
                if node.name in tree and tree[node.name].parents:
                    nx, ny = layout[tree[node.name]]
                    py = layout[tree[node.name].parents[0]][1]
                    start = arg[node.name].data["ancestral"][0][0]
                    win.add_group(lines(color(0,1,0), 
                                        x+nx, ny, x+nx, py,
                                        color(1,1,1)))
            
                
        x += treewidth

    #win.set_visible(* win.get_root().get_bounding() + ("exact",))
    win.home("exact")

    return win



def show_coal_track(tree_track):
    
    win = summon.Window()
    
    
    bgcolor = (1, 1, 1, .1)
    cmap = util.rainbow_color_map(low=0.5, high=1.0)

    maxage = 0
    for (start, end), tree in tree_track:
        print start
        l = []
        times = treelib.get_tree_timestamps(tree)
        nleaves = len(tree.leaves())
        maxage2 = 0
        for node in tree:            
            if len(node.children) > 1:
                age = times[node]
                sizes = [len(x.leaves()) for x in node.children]
                bias = max(sizes) / float(sum(sizes))
                l.extend([color(*cmap.get(bias)), start, age, end, age])
                if age > maxage2:
                    maxage2 = age
        win.add_group(group(lines(*l), color(*bgcolor),
                      box(start, 0, end, maxage2, fill=True)))
        if maxage2 > maxage:
            maxage = maxage2

    def func():
        x, y = win.get_mouse_pos()
        print "pos=%s age=%f" % (util.int2pretty(int(x)), y)
    win.add_group(hotspot("click", 0, 0, end, maxage,
                          func))
    
    win.home("exact")


    return win


def show_coal_track2(tree_track):
    
    win = summon.Window()
    
    
    bgcolor = (1, 1, 1, .1)
    cmap = util.rainbow_color_map(low=0.0, high=1.0)

    maxage = 0
    for (start, end), tree in tree_track:
        print start
        l = []
        times = treelib.get_tree_timestamps(tree)
        nleaves = len(tree.leaves())
        maxage2 = 0
        for node in tree:            
            if len(node.children) > 1:
                age = times[node]
                freq = len(node.leaves()) / float(nleaves)
                bias = max(sizes) / float(sum(sizes))
                l.extend([color(*cmap.get(freq)), start, age, end, age])
                if age > maxage2:
                    maxage2 = age
        win.add_group(group(lines(*l), color(*bgcolor),
                      box(start, 0, end, maxage2, fill=True)))
        if maxage2 > maxage:
            maxage = maxage2

    def func():
        x, y = win.get_mouse_pos()
        print "pos=%s age=%f" % (util.int2pretty(int(x)), y)
    win.add_group(hotspot("click", 0, 0, end, maxage,
                          func))
    
    win.home("exact")


    return win



'''
def mark_tree(tree, layout, name, y=None, time=None,
              col=(1, 0, 0), ymap=lambda y: y):
    nx, ny = layout[tree[name]]
    if y is not None:
        y += ny
    else:
        y = time
    return draw_mark(nx, ymap(y), col=col)
'''    


def draw_tree(tree, layout, orient="vertical"):
    
    vis = group()
    bends = {}

    for node in tree.postorder():
        # get node coordinates
        nx, ny = layout[node]
        px, py = layout[node.parents[0]] if node.parents else (nx, ny)

        # determine bend point
        if orient == "vertical":
            bends[node] = (nx, py)
        else:
            bends[node] = (px, ny)
        
        # draw branch
        vis.append(lines(nx, ny, bends[node][0], bends[node][1]))

        # draw cross bar
        if len(node.children) > 0:
            a = bends[node.children[-1]]
            b = bends[node.children[0]]
            vis.append(lines(a[0], a[1], b[0], b[1]))

    return vis


def draw_mark(x, y, col=(1,0,0), size=.5, func=None):

    if func:
        h = hotspot("click", x-size, y-size, x+size, y+size, func)
    else:
        h = group()
    
    return zoom_clamp(
        color(*col),
        box(x-size, y-size, x+size, y+size, fill=True),
        h,
        color(1,1,1),
        origin=(x, y),
        minx=10.0, miny=10.0, maxx=20.0, maxy=20.0,
        link=True)



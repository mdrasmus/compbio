#!/usr/bin/env summon


# python libs
import sys, math

# rasmus libs
from rasmus import genomeio, util, fasta, clustalw, muscle, phylip
from rasmus import algorithms, phyloutil

# graphics libs
from summon import *
import summonlib.tree


sumtree = summonlib.tree

options = sumtree.options + [
 ["d:", "desc=", "desc", "AUTO<descriptions>"],
 ["o:", "orth=", "orths", "AUTO<orths>"],
 ["c:", "orthcomp=", "orthcomps", "AUTO<orth comps>"],
 ["f:", "fasta=", "fasta", "AUTO<fasta>"],
 ["v:", "visual=", "visual", "AUTO<outputfile>"],
 ["w:", "winsize=", "winsize", "AUTO<window width>x<window height>"],
 ["r:", "reconroot=", "reconroot", "AUTO<species tree>"]
]


try:
    param, rest = util.parseArgs(sys.argv[1:], options)
except:
    sys.exit(1)


class VisTree (sumtree.SumTree):
    def __init__(self):
        sumtree.SumTree.__init__(self)
        
        self.colorOrths = color(0,0,0)
        self.clickMode = "desc"
        self.alignfunc = muscle.muscle

        self.desc = {}
        self.orths = []
        self.orthcomps = []
        self.seqs = {}


    # override node clicks
    def nodeClick(self, node):
        print "-------------"

        if self.clickMode == "desc":
            self.printNode(node)
            print "%d genes total" % len(self.tree.leaves(node))
        if self.clickMode == "gene":
            self.printNode(node, False)
            print "%d genes total" % len(self.tree.leaves(node))
        elif self.clickMode == "align":
            self.alignNode(node)
        elif self.clickMode == "refine":
            self.refine(node)


    def printNode(self, node, showDesc = True):
        if node.isLeaf():
            if showDesc and node.name in self.desc:
                print "%s\t%s" % (node.name, self.desc[node.name])
            else:
                print node.name
        else:
            for child in node.children:
                self.printNode(child, showDesc)


    def alignNode(self, node):
        # get all protein sequences of sub tree
        seqs2 = util.subdict(self.seqs, self.tree.leaveNames(node))

        self.selalign = self.alignfunc(seqs2, verbose = True)
        clustalw.printAlign(self.selalign, order=tree.leaveNames(node))


    def refine(self, node):
        #seqs2 = util.subdict(seqs, tree.leaveNames(node))
        #aln = alignfunc(seqs2, verbose = True)

        phylip.refineNode(self.tree, node, self.seqs, 
                          phylip.protpars, muscle.muscle)
        clear_groups()

        self.setupTree(self.param, self.tree)
        self.drawTree(self.param, self.tree, self.orths)


    def drawOrthologs(self, param, tree, orths):      
        vis = []

        for orth in orths:
            if orth[0] in tree.nodes and orth[1] in tree.nodes:
                x1 = tree.nodes[orth[0]].x
                y1 = tree.nodes[orth[0]].y
                x2 = tree.nodes[orth[1]].x
                y2 = tree.nodes[orth[1]].y

                length = math.sqrt((x1 - x2)**2 + (y1 - y2)**2)

                alpha = 1 / length
                alpha = max(alpha, .1)

                vis.append(lines(color(0,0,0,alpha), vertices(x1, y1, x2, y2)))

        return list2group(vis)


    def display(self, param, tree):
        sumtree.SumTree.display(self, param, tree)
        add_group(self.drawOrthologs(param, tree, self.orths))



def find(gene):
    return vis.find(gene)

def mark(boxColor = color(1,0,0)):
    return vis.mark(boxColor)

def flag(flagColor = color(1,0,0)):
    return vis.flag(flagColor)

def clear_marks():
    return vis.clearMarks()


vis = VisTree()

# read tree
tree = vis.readTree(param)

if "reconroot" in param:
    stree = algorithms.Tree()
    stree.readNewick(param["reconroot"][-1])
    tree = phyloutil.reconRoot(tree, stree)
    vis.setupTree(param, tree)


# read descriptions
if "desc" in param:
    util.tic("read descriptions")
    for filename in param["desc"]:
        util.log("reading '%s'..." % filename)
        vis.desc.update(genomeio.readGeneDesc(filename))
    util.log("loaded %d desciptions" % len(vis.desc))
    util.toc()

# read orthologs
if "orths" in param:
    util.tic("read orthologs")
    for filename in param["orths"]:
        util.log("reading '%s'..." % filename)
        vis.orths.extend(util.readDelim(filename))
    util.log("loaded %d orthologs" % len(vis.orths))
    util.toc()

if "orthcomps" in param:
    util.tic("read ortholog components")
    for filename in param["orthcomps"]:
        util.log("reading '%s'..." % filename)
        vis.orthcomps.extend(util.readDelim(filename))
    for comp in vis.orthcomps:
        for i in xrange(len(comp)):
            for j in xrange(i+1, len(comp)):
                vis.orths.append([comp[i], comp[j]])
    util.log("loaded %d ortholog components" % len(vis.orthcomps))
    util.toc()

# read sequence
if "fasta" in param:
    util.tic("read sequence")
    for filename in param["fasta"]:
        util.tic("read '%s'" % filename)
        vis.seqs.update(fasta.readFasta(filename))
        util.toc()
    util.log("loaded %d sequences" % len(vis.seqs))
    util.toc()

if "winsize" in param:
    w, h = map(int, param["winsize"][-1].split("x"))
    set_window_size(w, h)



# draw
vis.display(param, tree)
home()

if "visual" in param:
    import summonlib
    x1, y1, x2, y2 = get_visible()
    margin = (x2-x1) * .1
    set_visible(x1-margin, y1-margin, x2+margin, y2+margin)
    from summonlib import svg
    svg.printScreen(param["visual"][-1])
    sys.exit(0)
    


# add key bindings
def press(mode):
    def func():
        print "mode is '"+ mode + "'"
        vis.clickMode = mode
    return func

set_binding(input_key("a"), press("align"))
set_binding(input_key("d"), press("desc"))
set_binding(input_key("g"), press("gene"))
set_binding(input_key("r"), press("refine"))


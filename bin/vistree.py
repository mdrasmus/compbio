#!/usr/bin/python -i


# python libs
import sys, math

# rasmus libs
from rasmus import treelib
from rasmus import util
from rasmus.bio import fasta, clustalw, muscle, phylip, alignlib

# graphics libs
from summon.core import *
from summon import sumtree
from summon import svg
from summon import hud
import summon


options = [
# ["p:", "ptree=", "ptree", "<ptree file>",
#    {"single": True}],
# ["l:", "labels=", "labels", "<leaf labels file>",
#    {"single": True}],
 ["t:", "usedist=", "usedist", "<distance factor>",
    {"single": True,
     "default": 1.0,
     "parser": float}],
 ["L", "noshowlabels", "noshowlabels", "",
    {"single": True}],
 ["V", "vertical", "vertical", "",
    {"single": True}],

 ["o:", "orth=", "orths", "<orths>"],
 ["c:", "orthcomp=", "orthcomps", "<orth comps>"],
 
 ["f:", "fasta=", "fasta", "<fasta>"],
 ["v:", "visual=", "visual", "<outputfile>",
    {"single": True}],
 ["w:", "winsize=", "winsize", "<window width>x<window height>",
    {"single": True}],
 ["r:", "reroot=", "reroot", "<root node>",
    {"single": True}]
]


conf = util.parseOptions(sys.argv, options, resthelp="<newick tree> ...")

selnode = None

class VisTree (sumtree.SumTree):
    def __init__(self, conf, trees, filenames, **options):
        sumtree.SumTree.__init__(self, trees[0], **options)
        
        self.trees = trees
        self.filenames = filenames
        
        self.conf = conf
        self.colorOrths = color(0,0,0)
        self.clickMode = "desc"
        self.alignfunc = muscle.muscle

        self.desc = {}
        self.orths = []
        self.orthcomps = []
        self.seqs = {}
        self.selnode = None
        self.curtree = 0
    


    def switchTree(self, step=1):
        
        self.curtree = (self.curtree + step) % len(self.trees)
        self.tree = self.trees[self.curtree]
        self.show()
        

    # override node clicks
    def nodeClick(self, node):
        self.selnode = node
        
        # set global selnode variable for interactive users
        global selnode        
        selnode = node
        
        
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
        if node.is_leaf():
            if showDesc and node.name in self.desc:
                print "%s\t%s" % (node.name, self.desc[node.name])
            else:
                print node.name
        else:
            for child in node.children:
                self.printNode(child, showDesc)


    def alignNode(self, node):
        # get all protein sequences of sub tree
        seqs2 = self.seqs.get(self.tree.leaf_names(node))

        self.selalign = self.alignfunc(seqs2, verbose = True)
        self.selalign.names = tree.leaf_names(node)
        alignlib.printAlign(self.selalign)


    """
    def refine(self, node):
        #seqs2 = util.subdict(seqs, tree.leaf_names(node))
        #aln = alignfunc(seqs2, verbose = True)

        phylip.refineNode(self.tree, node, self.seqs, 
                          phylip.protpars, muscle.muscle)
        self.show()
    """
    

    def drawOrthologs(self, conf, tree, orths):      
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

                vis.append(lines(color(0,0,0,alpha), x1, y1, x2, y2))

        return group(*vis)


    def show(self):
        sumtree.SumTree.show(self)
            
        # draw ortholog links
        self.win.add_group(self.drawOrthologs(self.conf, self.tree, self.orths))
        
        self.win.set_name(self.filenames[self.curtree])
        
        # add key bindings
        def keypress(mode):
            def func():
                print "mode is '"+ mode + "'"
                self.clickMode = mode
            return func

        
        # build sidebar menu
        self.bar = hud.SideBar(self.win, width=150)
        self.bar.addItem(hud.MenuItem("align mode (a)", keypress("align")))
        self.bar.addItem(hud.MenuItem("gene mode (d)", keypress("gene")))
        #self.bar.addItem(hud.MenuItem("refine mode (r)", keypress("refine")))
        self.bar.addItem(hud.MenuItem("next tree (n)", lambda: self.switchTree(step=1)))
        self.bar.addItem(hud.MenuItem("prev tree (p)", lambda: self.switchTree(step=-1)))
        
        # register key bindings
        self.win.set_binding(input_key("a"), keypress("align"))
        #self.win.set_binding(input_key("d"), keypress("desc"))
        self.win.set_binding(input_key("g"), keypress("gene"))
        #self.win.set_binding(input_key("r"), keypress("refine"))
        self.win.set_binding(input_key("n"), lambda: self.switchTree(step=1))
        self.win.set_binding(input_key("p"), lambda: self.switchTree(step=-1))


        




# read tree
treefilenames = conf["REST"]
trees = []
for filename in treefilenames:
    tree = treelib.Tree()
    tree.readNewick(filename)
    trees.append(tree)
    
    print "%s: %d nodes, %d leaves\n" % \
        (filename, len(tree.nodes), len(tree.leaves()))


# init visualization
vis = VisTree(conf, trees, treefilenames, 
                    xscale=conf["usedist"], 
                    showLabels=not conf["noshowlabels"],
                    vertical=conf["vertical"])

find = vis.find
mark = vis.mark
flag = vis.flag
clearmarks = vis.clearMarks


if "reroot" in conf:
    if conf["reroot"].isdigit():
        conf["reroot"] = int(conf["reroot"])
    tree = treelib.reroot(tree, conf["reroot"], newCopy=False)
    vis.setupTree(conf, tree)



# read orthologs
if "orths" in conf:
    util.tic("read orthologs")
    for filename in conf["orths"]:
        util.log("reading '%s'..." % filename)
        vis.orths.extend(util.readDelim(filename))
    util.log("loaded %d orthologs" % len(vis.orths))
    util.toc()

if "orthcomps" in conf:
    util.tic("read ortholog components")
    for filename in conf["orthcomps"]:
        util.log("reading '%s'..." % filename)
        vis.orthcomps.extend(util.readDelim(filename))
    for comp in vis.orthcomps:
        for i in xrange(len(comp)):
            for j in xrange(i+1, len(comp)):
                vis.orths.append([comp[i], comp[j]])
    util.log("loaded %d ortholog components" % len(vis.orthcomps))
    util.toc()

# read sequence
if "fasta" in conf:
    util.tic("read sequence")
    for filename in conf["fasta"]:
        util.tic("read '%s'" % filename)
        vis.seqs.update(fasta.readFasta(filename))
        util.toc()
    util.log("loaded %d sequences" % len(vis.seqs))
    util.toc()

if "winsize" in conf:
    w, h = map(int, conf["winsize"].split("x"))
    vis.win.set_window_size(w, h)



# draw
vis.show()


if "visual" in conf:
    x1, y1, x2, y2 = vis.win.get_visible()
    margin = (x2-x1) * .1
    vis.win.set_visible(x1-margin, y1-margin, x2+margin, y2+margin)
    svg.printScreen(vis.win, conf["visual"])
    sys.exit(0)
    


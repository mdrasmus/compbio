#!/usr/bin/env python

import sys
from rasmus import util, env, treelib
from rasmus.bio import phylo, genomeutil
from rasmus.vis import treesvg
from rasmus import tablelib


options = [
  ["t:", "tree=", "tree", "<newick file>",
    {"default": []}],
  ["l:", "scale=", "scale", "<scaling>",
    {"default": 20,
     "single": True,
     "parser": float}],
  ["m:", "minlen=", "minlen", "<minimum branch length>",
    {"default": 1,
     "single": True,
     "parser": int}],  
  ["M:", "maxlen=", "maxlen", "<maximum branch length>",
    {"default": 10000,
     "single": True,
     "parser": int}],
  ["i", "hist", "hist", "",  
    {"single": True,
     "help": "output histogram of tree topologies"}],
  ["", "hashes", "hashes", "",
    {"single": True}], 
  ["n", "names", "names", "",
    {"single": True,
     "help": "display internal node names"}],
  ["N", "newick", "newick", "",
    {"single": True,
     "help": "write newick format"}],
  ["", "len", "len", "",
    {"single": True,
     "help": "display branch lengths"}],
  ["r:", "reroot=", "reroot", "<branch to root tree>",
    {"single": True}],
  ["c:", "colormap=", "colormap", "<color map file>",
    {"single": True}],
  ["", "rootby=", "rootby", "dup|loss|duploss",
    {"single": True,
     "default": "duploss"}],
  ["d", "dump", "dump", "",
    {"single": True,
     "help": "covert to easy to parse format"}],
  ["H", "headings", "headings", "",
    {"single": True,
     "help": "show heading information above each tree"}],
  ["g:", "graphical=", "graphical", "<filename>|-",
    {"single": True}],
  ["G", "default-graphical", "default-graphical", "",
    {"single": True}],
  ["e", "events", "events", "",
    {"single": True}],
  ["", "trees=", "trees", "{trees}",
   {"single": True,
    "parser": util.shellparser,
    "default": []}]
] + genomeutil.options


# parse options
conf = util.parseOptions(sys.argv, options, quit=True, resthelp="<trees> ...")
genomeutil.readOptions(conf)
gene2species = conf["gene2species"]
if "stree" in conf:
    stree = conf["stree"]
else:
    stree = None

if "colormap" in conf:
    colormap = treelib.readTreeColorMap(conf["colormap"])
else:
    colormap = None


hashes = []

if "reroot" in conf:
    if conf["reroot"].isdigit():
        conf["reroot"] = int(conf["reroot"])

def iterTrees(treefile):
    ntrees = 0
    infile = util.openStream(treefile)
    
    while True:
        try:
            tree = treelib.readTree(infile)
            ntrees += 1
            yield tree
        except Exception, e:
            if ntrees < 1:
                print >>sys.stderr, e
            break

def processTree(tree):
    global gene2species

    if "stree" in conf and \
       "smap" in conf and "rootby" in conf:
        phylo.reconRoot(tree, stree, gene2species, 
                            rootby=conf["rootby"],
                            newCopy=False)
    
    elif "reroot" in conf:
        tree = treelib.reroot(tree, conf["reroot"])
        
    
    if conf["hist"] or conf["hashes"]:
        # only count the tree in histogram mode
    
        thash = phylo.hashTree(tree, gene2species)
        hashes.append(thash)
    
    
    elif conf["dump"]:
        # dump mode
        
        names = util.sort(tree.nodes.keys())
        
        print "# node  path-to-root"
        for name in names:
            path = []
            ptr = tree.nodes[name]
            while ptr != tree.root:
                ptr = ptr.parent
                path.append(ptr.name)
            
            print "%s\t%s" % (name, ",".join(map(str, path)))
        print 
        
        print "# node branch-length"
        for name in names:
            print "%s\t%f" % (name, tree.nodes[name].dist)
    
    else:
        # default mode: display tree
        
        if conf["headings"]:
            print
            print "------------------------------------------------"
            print "filename: %s" % treefile
            print "treelen:  %f" % sum(x.dist for x in tree.nodes.values())
            print
        
        
        # set default graphical settings
        if conf["default-graphical"]:
            conf.setdefault("graphical", "-")
            conf['scale'] = 500.0


        labels = {}
        
        for node in tree.nodes.values():
            labels[node.name] = ""
        
        if conf["events"]:
            assert stree != None and gene2species != None
            phylo.initDupLossTree(stree)
            phylo.countDupLossTree(tree, stree, gene2species)
            phylo.countAncestralGenes(stree)
            
            for node in stree:
                labels[node.name] = "%d" % node.data['genes']
                
                if node.data['dup'] > 0:
                    labels[node.name] += " +%d" % node.data['dup']
                    
                if node.data['loss'] > 0:
                    labels[node.name] += " -%d" %  node.data['loss']
            tree = stree
            stree = None
            gene2species = None
        
        # create branch labels
        for node in tree.nodes.values():

            # label distances
            if conf["len"]:
                labels[node.name] += "%f" % node.dist
            
            # label bootstraps
            if "boot" in node.data and node.data["boot"] != 0:
                if isinstance(node.data["boot"], int):
                    labels[node.name] = "(%d) %s" % (node.data["boot"], 
                                                     labels[node.name])
                else:
                    labels[node.name] = "(%.2f) %s" % (node.data["boot"], 
                                                       labels[node.name])

            # label node names
            if conf["names"] and not node.isLeaf():
                labels[node.name] = "[%s] %s" % (node.name, 
                                                 labels[node.name])
        
        if "graphical" in conf:            
            if conf["graphical"] == "-":
                treesvg.showTree(tree, labels=labels,
                                       xscale=conf["scale"],
                                       minlen=conf["minlen"],
                                       maxlen=conf["maxlen"],
                                       legendScale=True,
                                       colormap=colormap,
                                       stree=stree,
                                       gene2species=gene2species)
            else:
                treesvg.drawTree(tree, labels=labels,
                                       xscale=conf["scale"],
                                       minlen=conf["minlen"],
                                       maxlen=conf["maxlen"],
                                       filename=conf["graphical"],
                                       legendScale=True,
                                       colormap=colormap,
                                       stree=stree,
                                       gene2species=gene2species)
        elif conf["newick"]:
            tree.write()
        else:
            treelib.drawTree(tree, labels=labels,
                                scale=conf["scale"],
                                minlen=conf["minlen"],
                                maxlen=conf["maxlen"])




for treefile in (conf["REST"] + conf["tree"] + conf["trees"]):
    #try:
    #    tree = treelib.readTree(treefile)        
    #except Exception, e:
    #    print >>sys.stderr, "error reading '%s': %s" % (treefile, e)
    #    continue
    
    for tree in iterTrees(treefile):
        processTree(tree)
    
# display topology histogram
if conf["hist"]:
    histogram = tablelib.histTable(hashes)
    histogram.write(sys.stdout)


if conf["hashes"]:
    for thash in hashes:
        print thash


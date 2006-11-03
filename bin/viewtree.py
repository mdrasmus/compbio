#!/usr/bin/env python

import sys
from rasmus import util, algorithms, phyloutil, genomeutil, env, treelib
from rasmus.vis import treevis


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
  ["n", "names", "names", "",
    {"single": True,
     "help": "display internal node names"}],
  ["", "nolen", "nolen", "",
    {"single": True,
     "help": "do not display branch length"}],
  ["r:", "reroot=", "reroot", "<branch to root tree>",
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
  ["", "trees=", "trees", "{trees}",
   {"single": True,
    "parser": util.shellparser,
    "default": ""}]
] + genomeutil.options


# parse options
conf = util.parseOptions(sys.argv, options, quit=True, resthelp="<trees> ...")
genomeutil.readOptions(conf)
gene2species = conf["gene2species"]
if "stree" in conf:
    stree = conf["stree"]



counts = util.Dict(1, 0)


if "reroot" in conf:
    if conf["reroot"].isdigit():
        conf["reroot"] = int(conf["reroot"])


for treefile in (conf["REST"] + conf["tree"] + conf["trees"].split()):
    try:
        tree = algorithms.readTree(treefile)
    except Exception, e:
        print >>sys.stderr, "error reading '%s': %s" % (treefile, e)
        continue
    
    if "stree" in conf and \
       "smap" in conf:
        phyloutil.reconRoot(tree, stree, gene2species, 
                            rootby=conf["rootby"],
                            newCopy=False)
    
    
    if "reroot" in conf:
        tree = treelib.reroot(tree, conf["reroot"])
        
    
    if conf["hist"]:
        # only count the tree in histogram mode
    
        thash = phyloutil.hashTree(tree, gene2species)
        counts[thash] += 1
    
    
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
            
        
        labels = {}
        
        # create branch labels
        for node in tree.nodes.values():
            
            # label distances
            if conf["nolen"]:
                labels[node.name] = ""
            else:
                labels[node.name] = "%f" % node.dist
            
            # label bootstraps
            if "boot" in node.data and node.data["boot"] != 0:
                labels[node.name] = "(%d) %s" % (node.data["boot"], 
                                                 labels[node.name])
        
            # label node names
            if conf["names"] and not node.isLeaf():
                labels[node.name] = "[%s] %s" % (node.name, 
                                                 labels[node.name])
        
        if "graphical" in conf:
            if conf["graphical"] == "-":
                treevis.showTree(tree, labels=labels,
                                       xscale=conf["scale"],
                                       minlen=conf["minlen"],
                                       maxlen=conf["maxlen"])
            else:
                treevis.drawTree(tree, labels=labels,
                                       xscale=conf["scale"],
                                       minlen=conf["minlen"],
                                       maxlen=conf["maxlen"],
                                       filename=conf["graphical"],
                                       legendScale=True)
        else:
            algorithms.drawTree(tree, labels=labels,
                                scale=conf["scale"],
                                minlen=conf["minlen"],
                                maxlen=conf["maxlen"])


# display topology histogram
if conf["hist"]:
    mat = map(list, counts.items())
    mat.sort(key=lambda x: x[1], reverse=True)
    tot = float(sum(util.cget(mat, 1)))
    
    for row in mat:
        row.append(row[1] / tot)
    
    util.printcols(mat)

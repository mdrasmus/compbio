#!/usr/bin/env python

import sys
from rasmus import util, algorithms, phyloutil, genomeutil, env


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
  ["d", "dump", "dump", "",
    {"single": True,
     "help": "covert to easy to parse format"}],
  ["H", "headings", "headings", "",
    {"single": True,
     "help": "show heading information above each tree"}]
] + genomeutil.options


# parse options
conf = util.parseOptions(sys.argv, options, quit=True)
genomeutil.readOptions(conf)
gene2species = conf["gene2species"]
if "stree" in conf:
    stree = conf["stree"]



counts = util.Dict(1, 0)


for treefile in (conf[""] + conf["tree"]):
    try:
        tree = algorithms.readTree(treefile)
    except Exception, e:
        print >>sys.stderr, "error reading '%s': %s" % (treefile, e)
        continue
    
    if "stree" in conf and \
       "smap" in conf and \
       "tree" in conf:
        tree = phyloutil.reconRoot(tree, stree, gene2species)
    
    
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
            print treefile
            print
            
        
        labels = {}
        
        # create branch labels
        for node in tree.nodes.values():
            
            # label distances
            labels[node.name] = "%f" % node.dist
            
            # label bootstraps
            if "boot" in node.data:
                labels[node.name] = "(%d) %s" % (node.data["boot"], 
                                                 labels[node.name])
        
            # label node names
            if conf["names"] and not node.isLeaf():
                labels[node.name] = "[%s] %s" % (node.name, 
                                                 labels[node.name])
            
        algorithms.drawTree(tree, labels=labels,
                            scale=conf["scale"],
                            minlen=conf["minlen"],
                            maxlen=conf["maxlen"])


# display topology histogram
if conf["hist"]:
    util.printDictByValues(counts, compare=util.invcmp)

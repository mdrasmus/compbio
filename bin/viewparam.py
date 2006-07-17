#!/usr/bin/env python

import sys
from rasmus import util, env, treelib, sindirlib, stats


options = [
  ["s:", "stree=", "stree", "<species tree>",
    {"single": True}],
  ["p:", "params=", "params", "<sindir params file>",
    {"single": True}],
  ["l:", "scale=", "scale", "<scaling>",
    {"default": 20,
     "single": True,
     "parser": float}],
  ["c:", "curves=", "curves", "<curve height>",
    {"single": True,
     "parser": int,
     "default": 0}],
  ["o", "output", "output", "", {"single": True}]
]


# parse options
conf = util.parseOptions(sys.argv, options, quit=True)


# read data
env.addEnvPaths("DATAPATH")
stree = treelib.readTree(env.findFile(conf["stree"]))
params = sindirlib.readParams(env.findFile(conf["params"]))


labels = {}

# set branch lengths to means
for name in stree.nodes:
    stree.nodes[name].dist = params[name][0]
    
    nchars = int(params[name][1] * conf["scale"])
    blen   = int(params[name][0] * conf["scale"]) - 2
    labels[name] = (blen - nchars) * "-" + nchars * "=" + " "
    
    # label internal nodes
    if not stree.nodes[name].isLeaf():
        strnum = str(name)
        labels[name] = labels[name][:-len(strnum)-2] + strnum + "= "

    if conf["curves"] > 0:
        label = ""
        maxy = stats.normalPdf(params[name][0], params[name])
        
        for y in range(conf["curves"]-1 , 0, -1):
            for x in range(blen):
                if stats.normalPdf(x / conf["scale"], params[name]) / \
                    maxy * (conf["curves"]) <= y:
                    
                    label += " "
                else:
                    label += "#"
            label += " \n"
        
        labels[name] = label + labels[name]


# simply output species trees with mean nt sub/site
if conf["output"]:
    stree.write()
    sys.exit(0)
    

# print params
keys = util.sort(util.remove(params.keys(), "baserate", 1)) + ["baserate"]
mat = [["Node", "Mean", "Sdev", "Mean/Sdev"]] + \
      [(key, params[key][0], params[key][1], params[key][0] / params[key][1])
       for key in keys]


util.printcols(mat, spacing=3)
print "avg. baserate:", params["baserate"][0] / params["baserate"][1]
print


# print tree
if conf["curves"] > 0:
    treelib.drawTree(stree, labels=labels,
                     labelOffset=-conf["curves"]+1,
                     spacing=conf["curves"],
                     scale=conf["scale"])

else:
    treelib.drawTree(stree, labels=labels,
                     labelOffset=0,
                     scale=conf["scale"])



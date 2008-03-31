#!/usr/bin/env python

import sys, os

from rasmus import treelib, util
import Spidir


options = [
    ["s:", "stree=", "stree", "", 
     {"single": True}],
    ["p:", "params=", "params", "", 
     {"single": True}],
    ["o:", "out=", "out", "",
     {"single": True}]
]

conf = util.parseOptions(sys.argv, options)


stree = treelib.readTree(conf["stree"])
params = Spidir.readParams(conf["params"])

generate = params["baserate"][0] / params["baserate"][1]

for name, (mu, sigma) in params.iteritems():
    if name in stree.nodes:
        stree.nodes[name].dist = generate * mu

stree.write(conf["out"])

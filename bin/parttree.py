#!/usr/bin/env python

from rasmus import util, cluster

import sys


options = [
    ["p:", "part=", "part", "AUTO<part file>"],
    ["c:", "cols=", "cols", "AUTO<genecol1,genecol2,scorecol>"],
    ["n:", "minscore=", "minscore", "AUTO<min score>"]
]

try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    sys.exit(1)

# parse args
parts = []
for f in param["part"]:
    parts.extend(util.readDelim(f))

matchfiles = rest

cols = map(int, param["cols"][-1].split(","))

minscore = 500
if "minscore" in param:
    minscore = float(param["minscore"][-1])

tree = cluster.partTreeFiles(parts, matchfiles, cols[0], cols[1], cols[2],
    minscore)

for node in tree.nodes.itervalues():
    node.dist = 1
    
tree.writeNewick(sys.stdout)

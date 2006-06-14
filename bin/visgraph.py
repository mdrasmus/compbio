#!/usr/bin/env python

from rasmus import graphviz, util
import sys


options = [
    ["e:", "edges=", "edges", "AUTO<edge file>"],
    ["o:", "out=", "out", "AUTO<output file>"]
    
]



# check args
try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    sys.exit(1)


mat = util.Dict(2, 0)

for line in file(param["edges"][-1]):
    vert1, vert2 = line.rstrip().split()[:2]
    mat[vert1][vert2] = 1


graphviz.visualize(mat, param["out"][-1])

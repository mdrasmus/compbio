#!/usr/bin/env python

import sys
from rasmus import util, graph

options = [
    ["p:", "part=", "part", "<part file>"],
    ["m:", "match=", "match", "<match file>"]
    ]


conf = util.parseOptions(sys.argv, options, quit=True)


vertices = {}


restparts = filter(lambda x: not x.endswith(".match"), conf["REST"])
restmatches = filter(lambda x: x.endswith(".match"), conf["REST"])


i = 0
if "part" in conf:
    for f in conf["part"] + restparts:
        for line in file(f):
            for word in line.split():
                vertices.setdefault(word, {})[i] = 1
                vertices.setdefault(i, {})[word] = 1
            i += 1

if "match" in conf:
    for f in conf["match"] + restmatches:
        for line in file(f):
            gene1, gene2 = line.rstrip().split()[:2]
            vertices.setdefault(gene1, {})[gene2] = 1
            vertices.setdefault(gene2, {})[gene1] = 1


def getNeighbors(vertex):
    return vertices[vertex].keys()

comps = graph.connectedComponents(vertices, getNeighbors)

for comp in comps:
    for word in comp:
        if type(word) == str:
            print word,
    print
    

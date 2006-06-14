#!/usr/bin/env python

import util
import sys
import graph

options = [
    ("g:", "go=", "go", "[-g <go file>] [--go <go file>]"),
    ("r", "rank", "rank", "")
    ]

try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    sys.exit(1)


banned = {
 "GO:0005634":1, "GO:0016021":1, "GO:0006355":1, "GO:0005515":1, "GO:0003677":1,
 "GO:0000004":1, "GO:0005554":1, "GO:0008372":1, "GO:0016020":1, "GO:0016740":1 
}

allowed = {
 "GO:0004984":1, "GO:0008270":1
}

vertices = {}

i = 0
for f in param["go"]:
    for line in file(f):
        for word in line.split():
            #if word not in allowed: #banned:
            #    continue
            vertices.setdefault(word, {})[i] = 1
            vertices.setdefault(i, {})[word] = 1
        i += 1

#keys = vertices.keys()
#for key in keys:
#    if len(vertices[key]) > 20:
#        vertices[key] = {}


def getNeighbors(vertex):
    return vertices[vertex].keys()

if "rank" in param:
    keys = vertices.keys()
    keys.sort(lambda a,b: cmp(len(vertices[b]), len(vertices[a])))

    for key in keys:
        print key, len(vertices[key])

    sys.exit(0)

comps = graph.connectedComponents(vertices, getNeighbors)


for comp in comps:
    comp = filter(lambda x: type(x) == str and x.startswith("ENS"), comp)
    if len(comp) > 1:
        for word in comp:
            print word,
        print
    

#!/usr/bin/python

import util
import getopt
import sys
import graph


if len(sys.argv) < 2:
    print "usage: genepart.py [-p partitions] [-s synteny] [-m matches]"
    sys.exit(1)


params, rest = util.getopt(sys.argv[1:], "m:p:s:")


# read partition
vertices = util.Dict(2)
parts = util.readDelim(params["-p"][0])
for i in xrange(len(parts)):
    for gene in parts[i]:
        pass
        #vertices[gene][i] = 1
        #vertices[i][gene] = 1

# read synteny
for f in params["-s"]:
    for line in file(f):
        genes = line.rstrip().split()
        for i in xrange(len(genes)):
            for j in xrange(i+1, len(genes)):
                vertices[genes[i]][genes[j]] = 1
                vertices[genes[j]][genes[i]] = 1

# merge partitions that share syntenic genes
#def getNeighbors(vertex):
#    return vertices[vertex].keys()
#comps = graph.connectedComponents(vertices.keys(), getNeighbors)

#print len(comps)
#print map(len, comps)

sys.exit(0)

# determine gene seeds
seeds = {}
i = 0
for part in parts:
    for gene in part:
        seeds[gene] = i
    i += 1

keep = {}

# pull in matches
for matchfile in params["-m"]:
    for line in file(matchfile):
        (gene1, gene2, score) = line.split()[:3]
        score = float(score)
        
        if not gene1 in seeds and gene2 in seeds:
            tmp = gene1; gene1 = gene2; gene2 = tmp
        
        if gene1 in seeds and not gene2 in seeds:
            if (not gene2 in keep) or (score > keep[gene2][1]):
                keep[gene2] = (gene1, score)


# add kept genes to parts
for gene in keep:
    parts[seeds[keep[gene][0]]].append(gene)

# output
for part in parts:
    for gene in part:
        print gene,
    print


#!/usr/bin/python

import util
import getopt
import sys
import graph
import algorithms

       


if len(sys.argv) < 2:
    print "usage: genepart.py [-p partitions] [-s synteny] [-m matches]"
    print "       [-maxpart=<size>]"
    sys.exit(1)


params, rest = util.getopt(sys.argv[1:], "m:p:s:", ["maxpart="])

maxpart = int(params["--maxpart"][-1])


vertices = util.Dict(2)

# read partition
sets = {}
parts = util.readDelim(params["-p"][0])
for i in xrange(len(parts)):
    set = algorithms.Set()
    for gene in parts[i]:
        #vertices[gene][i] = 1
        #vertices[i][gene] = 1
        sets[gene] = set
        set.add(gene)


# read synteny
if "-s" in params:
    for f in params["-s"]:
        for line in file(f):
            genes = line.rstrip().split()
            for i in xrange(len(genes)):
                for j in xrange(i+1, len(genes)):
                    vertices[genes[i]][genes[j]] = 1
                    vertices[genes[j]][genes[i]] = 1


# pull in best unidirectional matches
for matchfile in params["-m"]:
    print >>sys.stderr, matchfile
    for line in file(matchfile):
        (gene1, gene2, score) = line.split()[:3]
        score = float(score)
        
        if gene1 == gene2:
            continue
        
        if gene1 not in sets:
            set1 = algorithms.Set()
            set1.add(gene1)
            sets[gene1] = set1
        else:
            set1 = sets[gene1]

        if gene2 not in sets:
            set2 = algorithms.Set()
            set2.add(gene2)
            sets[gene2] = set2
        else:
            set2 = sets[gene2]
        
        if set1.same(set2) or \
           set1.size() + set2.size() > maxpart:
            continue
        
        #size1 = set1.size()
        #size2 = set2.size()
        set1.union(set2)
        #assert set1.size() == size1 + size2
        #assert set1.size() <= maxpart, (size1, size2)


# find unique sets
lst = {}
for set in sets.values():
    lst[set.root()] = 1

for set in lst:
    for gene in set.members():
        print gene,
    print


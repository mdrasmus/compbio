#!/usr/bin/env python

import util
import sys


options = [
    ["p:", "part=", "part", "AUTO<part file>"],
    ["m:", "match=", "match", "AUTO<match file>"],
    ]

try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    print sys.exc_type
    sys.exit(1)




# read in partition and make a gene->part lookup
parts = util.readDelim(param["part"][-1])
lookup = {}
for i in range(len(parts)):
    for gene in parts[i]:
        lookup[gene] = i

outgroups = [None] * len(parts)


# stream in matches
for f in param["match"]:
    print >>sys.stderr, f
    
    for line in file(f):
        gene1, gene2, score = line.split()[:3]
        score = float(score)
        
        # don't use new dog genes (their names are numbers)
        # that don't have sequence
        if gene1[0].isdigit() or gene2[0].isdigit():
            continue
        
        if gene1 in lookup and gene2 in lookup:
            i = lookup[gene1]
            j = lookup[gene2]
            
            # do not use matches within same partition
            if i == j:
                continue
            
            if not outgroups[i] or outgroups[i][1] < score:
                outgroups[i] = [gene2, score]
            if not outgroups[j] or outgroups[j][1] < score:
                outgroups[j] = [gene2, score]
            
for i in outgroups:
    try:
        print i[0]
    except:
        print "ERROR with", i


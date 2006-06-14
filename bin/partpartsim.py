#!/usr/bin/env python

from rasmus import util
import sys


options = [
    ["p:", "part=", "part", "AUTO<part file>"],
    ["m:", "match=", "match", "AUTO<match file>"],
    ["n:", "min=", "min", "AUTO<min score>"],
    ["c:", "cols=", "cols", "AUTO<col1>,<col2>,<scorecol>"]
    ]

try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    print sys.exc_type
    sys.exit(1)



mat = util.Dict(2, [0, 0])

parts = util.readDelim(param["part"][-1])
lookup = {}
for i in range(len(parts)):
    for gene in parts[i]:
        lookup[gene] = i

minscore = float(param["min"][-1])
col1, col2, scorecol = map(int, param["cols"][-1].split(","))

for f in param["match"]:
    print >>sys.stderr, f
    
    for line in file(f):
        tokens = line.split()
        gene1 = tokens[col1]
        gene2 = tokens[col2]
        score = float(tokens[scorecol])
        
        if score < minscore:
            continue
        
        if gene1 in lookup and gene2 in lookup:
            i = lookup[gene1]
            j = lookup[gene2]
            
            if i > j:
                tmp = i; i = j; j = tmp
            
            avg = mat[i][j]
            avg[0] += 1
            avg[1] += score
            

for key1 in mat:
    for key2 in mat[key1]:
        avg = mat[key1][key2]
        print key1, key2, (avg[1] / avg[0])


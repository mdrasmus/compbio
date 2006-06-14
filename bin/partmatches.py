#!/usr/bin/env python

from rasmus import util, ensembl, genomeio, descriptions, genomeutil, fasta
import sys



options = [
    ["p:", "part=", "part", "AUTO<part file>"],
    ["c:", "cols=", "cols", "AUTO<gene1 col>,<gene2 col>,<score col>"],
    ["o:", "outdir=", "outdir", "AUTO<output directory>"]
]

try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    print sys.exc_type
    sys.exit(1)


parts = util.readDelim(param["part"][-1])
cols = map(int, param["cols"][-1].split(","))

lookup = {}
for i in xrange(len(parts)):
    for gene in parts[i]:
        lookup[gene] = i

matches = []
for i in xrange(len(parts)):
    matches.append([])

# read in matches
for fn in rest:
    util.log(fn)
    for line in file(fn):
        tokens = line.rstrip().split()
        gene1 = tokens[cols[0]]
        gene2 = tokens[cols[1]]
        
        try:
            if lookup[gene1] == lookup[gene2]:
                matches[lookup[gene1]].append(line)
                
        except KeyError:
            pass


# write out matches
for i in xrange(len(parts)):
    out = file(param["outdir"][-1] + "/%d.match" % i, "w")
    
    for match in matches[i]:
        print >>out, match,
    out.close()




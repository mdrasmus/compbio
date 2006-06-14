#!/usr/bin/env python

from rasmus import genomeio, util
import sys

options = [
 ["e:", "ensmart=", "ensmart", "AUTO<ensmart file>"],
 ["o:", "output=", "output", "AUTO<output file>"]
]

try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    sys.exit(1)


infilename = param["ensmart"][-1]
outprefix = param["output"][-1]


infile = file(infilename)

# parse header
header = infile.next().rstrip()
fields = header.split("\t")
chromcol = fields.index("Chromosome")
chroms = {}

for line in infile:
    fields = line.split("\t")
    chrom = fields[chromcol]
    
    #if chrom.find("NT_") != -1:
    #    continue
    
    if chrom not in chroms:
        chroms[chrom] = file(outprefix + chrom + ".struct", "w")
        print >>chroms[chrom], header
    chroms[chrom].write(line)





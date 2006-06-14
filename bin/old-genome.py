#!/usr/bin/env python

import sys

name = sys.argv[1] + "_"

print "Ensembl Gene ID	Chromosome	Start Position (bp)	End Position (bp)	Strand"

for line in sys.stdin:
    tokens = line.rstrip().split("\t")
    
    print name+tokens[0] + "\t" + "\t".join(tokens[1:4]),
    if tokens[4] == "+":
        print "\t1"
    else:
        print "\t-1"
    

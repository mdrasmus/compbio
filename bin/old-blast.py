#!/usr/bin/env python

import sys


name1 = sys.argv[1] + "_"
name2 = sys.argv[2] + "_"

sys.stdin.next()

for line in sys.stdin:
    tokens = line.rstrip().split("\t")
    
    if "_" in tokens[0]:
        gene1 = name1 + tokens[0].split("_")[1]
    else:
        gene1 = tokens[0]
    
    if "_" in tokens[1]:
        gene2 = name2 + tokens[1].split("_")[1]
    else:
        gene2 = tokens[1]
    score = tokens[8]
    
    print "\t".join([gene1, gene2, score])
    

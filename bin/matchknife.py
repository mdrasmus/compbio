#!/usr/bin/python

import sys
from rasmus import util


class Param:
    def __init__(self):
        self.matches = []
        self.sel1 = {}
        self.sel2 = {}
        self.minscore = 0
        self.method = "and"
    
param = Param()

def addSel(sel, genes):
    for gene in genes:
        sel[gene] = 1


i = 1
while (i < len(sys.argv)):   
    if sys.argv[i] == "-sel1":
        addSel(param.sel1, util.readStrings(sys.argv[i+1]))
        i += 2
    elif sys.argv[i] == "-sel2":
        addSel(param.sel2, util.readStrings(sys.argv[i+1]))
        i += 2
    elif sys.argv[i] == "-or":
        param.method = "or"
        i += 1
    elif sys.argv[i] == "-min":
        param.minscore = float(sys.argv[i+1])
        i += 2
    elif sys.argv[i] == "-m":
        param.matches.append(sys.argv[i+1])
        i += 2
    else:
        param.matches.append(sys.argv[i])
        i += 1

for match in param.matches:
    print >>sys.stderr, match
    
    if param.method == "and":
        for line in file(match):
            (gene1, gene2, score) = line.split()[:3]
            score = float(score)
            
            if score < param.minscore:
                continue
            
            if not gene1.startswith("ENS") or \
               not gene2.startswith("ENS"):
                continue
            
            if (len(param.sel1) == 0 or gene1 in param.sel1) and \
               (len(param.sel2) == 0 or gene2 in param.sel2):
                print line,
    else:
        for line in file(match):
            (gene1, gene2, score) = line.split()[:3]
            score = float(score)
            
            if score < param.minscore:
                continue
                            
            if gene1 in param.sel1 or gene2 in param.sel2:
                print line,



#!/usr/bin/env python

from rasmus import util
import sys

options = [
  ["n:", "minscore=", "minscore", "AUTO<minscore>"],
]

try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    sys.exit(1)




def blastp2match(infile, outfile, minscore):
    query = ""
    inhits = False
    
    mat = {}
    
    for line in infile:
        if line.startswith("Query="):
            query = line.split("=")[1].strip()
        elif line.startswith("Sequences"):
            inhits = True
        elif inhits:
            if line.startswith(">"):
                inhits = False
            else:
                tokens = line.split()
                if len(tokens) > 0:
                    subject, score = tokens[:2]
                    if float(score) > minscore and \
                     query != subject and \
                     query + subject not in mat:
                        mat[query + subject] = 1
                        print >>outfile, "\t".join([query, subject, score])


minscore = 0
if "minscore" in param:
    minscore = float(param["minscore"][-1])


for fin in rest:
    fout = fin.replace(".blastp", ".match")
    blastp2match(file(fin), file(fout, "w"), minscore)


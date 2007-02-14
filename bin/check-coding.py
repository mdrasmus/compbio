#!/usr/bin/env python

from rasmus import fasta, util, alignlib
import sys


options = [
 ["d:", "dna=", "dna", "<dna fasta>"],
 ["p:", "prot=", "prot", "<protein fasta>"] 
]


param = util.parseOptions(sys.argv, options)

dna  = {}
prot = {}

for f in param["dna"]:
    dna.update(fasta.readFasta(f))

for f in param["prot"]:
    prot.update(fasta.readFasta(f))

good = True

for key in prot:
    if key not in dna:
        print "'%s' not in coding" % key
        good = False
    
    try:
        aa = alignlib.translate(dna[key])
    
        if prot[key] != aa:
            good = False
            print "mismatch", key
            #print aa
            #print prot[key]
            #print 
    except:
        print "error", key
        good = False

if good:
    print "all sequences pass!"
else:
    print "errors occurred"

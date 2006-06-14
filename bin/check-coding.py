#!/usr/bin/env python

from rasmus import fasta, util, genomeutil
import sys


options = [
 ["d:", "dna=", "dna", "AUTO<dna fasta>"],
 ["p:", "prot=", "prot", "AUTO<protein fasta>"] 
]

try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    sys.exit(1)

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
    aa = genomeutil.dna2aa(dna[key])
    
    if prot[key] != aa:
        print "mismatch '%s'" % key
        print aa
        print prot[key]
        print 
    else:
        sys.stdout.write(".")

if good:
    print "all sequences pass!"
else:
    print "errors occurred"

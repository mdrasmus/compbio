#!/usr/bin/env python

import sys
from rasmus import util, fasta, alignlib


options = [
    ["f:", "fasta=", "fasta", "AUTO<fasta file>"],
    ["w:", "width=", "width", "AUTO<width of alignment>"],
    ]


param = util.parseOptions(sys.argv, options, quit=True)


aln = fasta.readFasta(param["fasta"][-1])


width = 59
if "width" in param:
    width = int(param["width"][-1])




percid = alignlib.calcConservation(aln)

print "length:  %d" % len(aln.values()[0])
print "perc id: %f" % (util.countge(.99, percid) / float(len(percid)))
print 

alignlib.printAlign(aln, seqwidth=width)

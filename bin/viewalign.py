#!/usr/bin/env python

import sys
from rasmus import util
from rasmus.bio fasta, alignlib


options = [
    ["f:", "fasta=", "fasta", "<fasta file>",
        {"default": []}],
    ["w:", "width=", "width", "<width of alignment>",
        {"single": True,
         "default": 59,
         "parser": int}],
    ["c:", "check=", "check", "<overlap percent>", 
        {"single": True,
         "default": None,
         "parser": float}],
    ]


conf = util.parseOptions(sys.argv, options, quit=True)

alnfiles = conf["fasta"] + conf["REST"]

def printAlign(conf, alnfile):       
    aln = fasta.readFasta(alnfile)

    percid = alignlib.calcConservation(aln)
    
    print "----------------------------------"
    print "file:    %s" % alnfile
    print "length:  %d" % len(aln.values()[0])
    print "perc id: %f" % (util.countge(.99, percid) / float(len(percid)))
    print 

    alignlib.printAlign(aln, seqwidth=conf["width"])



# main
for alnfile in alnfiles:
    if conf['check'] != None:
        aln = fasta.readFasta(alnfile)
        alignlib.checkAlignOverlap(aln, conf['check'])
    else:
        printAlign(conf, alnfile)
    
 

#!/usr/bin/env python

import sys
from rasmus import util
from rasmus.bio import fasta, alignlib


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
    ["", "codon", "codon", "",
     {"single": True,
      "help": "use codon alignment view"}],
    ]


conf = util.parseOptions(sys.argv, options, quit=True)

alnfiles = conf["fasta"] + conf["REST"]


def codonAlign(aln):

    aln2 = fasta.FastaDict()

    assert aln.alignlen() % 3 == 0

    for key, seq in aln.iteritems():
        seq2 = []

        for i in xrange(0, len(seq), 3):
            seq2.append(seq[i:i+3])
        
        aln2[key] = " ".join(seq2)

    return aln2
        

def printAlign(conf, alnfile):       
    aln = fasta.readFasta(alnfile)

    percid = alignlib.calcConservation(aln)
    
    print "----------------------------------"
    print "file:    %s" % alnfile
    print "length:  %d" % len(aln.values()[0])
    print "perc id: %f" % (util.countge(.99, percid) / float(len(percid)))
    print 

    if conf["codon"]:
        aln = codonAlign(aln)

    alignlib.printAlign(aln, seqwidth=conf["width"])



# main
for alnfile in alnfiles:
    if conf['check'] != None:
        aln = fasta.readFasta(alnfile)
        alignlib.checkAlignOverlap(aln, conf['check'])
    else:
        printAlign(conf, alnfile)
    
 

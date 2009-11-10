#!/usr/bin/env python

import sys, optparse

from rasmus import util
from rasmus.bio import fasta, alignlib, seqlib


#=============================================================================
o = optparse.OptionParser()
o.add_option("-d", "--dna", dest="dna", action="append", default=[])
o.add_option("-p", "--pep", dest="pep", action="append", default=[])
o.add_option("--candida", dest="candida", action="store_true",
             help="use candida coding for Lucine and Serine")
o.add_option("--nols", dest="nols", action="store_true",
             help="ignore Lucine and Serine")


conf, args = o.parse_args()


#=============================================================================
table = seqlib.CODON_TABLE

if conf.candida:
    table = seqlib.CANDIDA_CODON_TABLE


dna  = {}
prot = {}

for f in conf.dna:
    dna.update(fasta.read_fasta(f))

for f in conf.pep:
    prot.update(fasta.read_fasta(f))

good = True

for key in prot:
    if key not in dna:
        print "%s\tnot in coding" % key
        good = False

    if len(dna[key]) != len(prot[key]) * 3:
        print "%s\twrong length\t%d != %d" % (key, len(dna[key]),
                                              len(prot[key])*3)
        print "%s\twrong length pep\t%s" % (key, prot[key])
        print "%s\twrong length trans dna\t%s" % \
              (key,seqlib.translate(dna[key] + "N" * (3 - len(dna[key]) % 3)))
    
    try:
        aa = alignlib.translate(dna[key], table)
        aa2 = prot[key]

        if conf.nols:
            aa.replace("S", "L")
            aa2.replace("S", "L")
        
        if aa2 != aa:
            good = False
            print "%s\tmismatch" % key
    except Exception, e:
        print "%s\terror\t%s" % (key, str(e))
        good = False


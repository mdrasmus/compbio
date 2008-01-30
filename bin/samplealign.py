#!/usr/bin/env python

import os, sys, random

from rasmus import env, util
from rasmus.bio import fasta, genomeutil, alignlib


options = [
    ["l:", "len=", "len", "<length of new alignment>",
     {"req": True,
      "parser": int,
      "single": True}],
]

conf = util.parseOptions(sys.argv, options, 
                         resthelp="<alignment>")


if len(conf["REST"]) == 0:
    sys.exit(1)

aln = fasta.readFasta(conf["REST"][0])
newlen = conf["len"]
oldlen = aln.alignlen()

# sample columns
cols = [random.randint(0, oldlen)
        for i in xrange(newlen)]
    
aln2 = alignlib.subalign(aln, cols)
aln2.write()

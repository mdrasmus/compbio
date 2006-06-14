#!/usr/bin/python

from rasmus import fasta, util
import random, sys


if len(sys.argv) < 3:
    print "usage: fasta-sample.py <fasta> <sample size>"
    sys.exit(1)

filename, size = sys.argv[1:1+2]
size = int(size)

seqs = fasta.readFasta(filename)

keys = seqs.keys()
sample = random.sample(keys, size)

seqs2 = util.getkeys(seqs, sample)

fasta.writeFasta(sys.stdout, seqs2)

    

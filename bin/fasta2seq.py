#!/usr/bin/env python

import sys


if len(sys.argv) < 3:
    print "usage: fasta2seq <fasta> <seq>"
    sys.exit(1)

(infilename, outfilename) = sys.argv[1:3]

infile = file(infilename)

# read header
line = infile.next()
assert line.startswith(">")

out = file(outfilename, "w")
for line in infile:
    if line[0] == ">":
        print "ERROR: multiple sequences!!!"
        sys.exit(1)
    out.write(line.rstrip())

out.close()

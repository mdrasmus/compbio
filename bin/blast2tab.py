#!/usr/bin/env python

from rasmus import util
from rasmus.bio import blast
import sys

options = [
    ["b:", "blast=", "blast", "AUTO<blast file>"],
    ["o:", "out=", "out", "AUTO<output file>"],
    ["c", "cut", "cut", "[-c] [--cut]"]
]


try:
    param, rest = util.parseOptions(sys.argv, options)
except:
    sys.exit(1)


# execute
infile = param["blast"][-1]
outfile = param["out"][-1]


if "cut" in param:
    util.tic("cutting file")
    outs = blast.blastCut(file(infile), util.tempfile(".", "blast2tab", ".tmp"))
    util.toc()
else:
    outs = [infile]


out = file(outfile, "w")
for infile in outs:
    util.tic("converting '%s' to tab file" % infile)
    blast.blast2tab(infile, out)
    util.toc()
out.close()

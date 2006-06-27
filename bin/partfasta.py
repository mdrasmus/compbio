#!/usr/bin/env python

import os
import sys

from rasmus import env
from rasmus import fasta
from rasmus import genomeutil
from rasmus import util


SCRIPT = sys.argv[0]

options = [
    ["p:", "part=", "part", "<part file>"],
    ["f:", "fasta=", "fasta", "<fasta file>"],    
    ["o:", "outdir=", "outdir", "<output directory>",
        {"single": True}],
    ["P:", "paths=", "paths", "<data path1>:<data path2>:...", 
        {"default": ".",
         "single": True}],
    ["d", "dir", "dir", "",
        {"single": True,
         "help": "output each fasta in its own directory"}],
]


conf = util.parseOptions(sys.argv, options, quit=True)


def main(conf):
    # setup paths
    env.addEnvPaths("DATAPATH")
    env.addPaths(conf["paths"])
    
    seqs = fasta.FastaDict()
    nseqs = 0
    util.tic("read sequences")
    for f in conf["fasta"]:
        util.tic("reading '%s'" % f)
        seqs.read(env.findFile(f), useIndex=False)
        util.log("read %d sequences" % (len(seqs.keys()) - nseqs))
        nseqs = len(seqs.keys())
        util.toc()
    util.toc()

    parts = util.readDelim(env.findFile(conf["part"][-1]))
    
    # actually split fasta
    for i in xrange(len(parts)):
        if conf["dir"]:
            if not os.path.exists(os.path.join(conf["outdir"], str(i))):
                os.mkdir(os.path.join(conf["outdir"], str(i)))
            seqfile = os.path.join(conf["outdir"], str(i), "%d.fasta" % i)
        else:
            seqfile = os.path.join(conf["outdir"], "%d.fasta" % i)
        util.log(seqfile)
        seqs2 = seqs.get(parts[i])
        seqs2.write(seqfile)


main(conf)

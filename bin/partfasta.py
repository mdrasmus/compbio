#!/usr/bin/env python

import os
import sys
from itertools import izip

from rasmus import env
from rasmus import fasta
from rasmus import genomeutil
from rasmus import util
from rasmus import genecluster
from rasmus import tablelib



options = [
    ["p:", "part=", "part", "<part file>"],
    ["t:", "famtab=", "famtab", "<family table"],
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
    
    # read partition
    if "part" in conf:
        parts = util.readDelim(env.findFile(conf["part"][-1]))
        famids = map(str, range(len(parts)))
    
    if "famtab" in conf:
        famtab = tablelib.readTable(conf["famtab"][-1])
        famids = famtab.cget("famid")
        parts = genecluster.famtab2parts(famtab)
    
    
    # actually split fasta
    for famid, part in izip(famids, parts): #xrange(len(parts)):
    
        if conf["dir"]:
            if not os.path.exists(os.path.join(conf["outdir"],famid)):
                os.mkdir(os.path.join(conf["outdir"], famid))
            seqfile = os.path.join(conf["outdir"], famid, "%s.fasta" % famid)
        else:
            seqfile = os.path.join(conf["outdir"], "%s.fasta" % famid)
        util.logger(seqfile)
        seqs2 = seqs.get(part)
        seqs2.write(seqfile)


main(conf)

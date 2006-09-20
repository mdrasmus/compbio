#!/usr/bin/env python

import sys, time, os
from rasmus import synteny, synphylweb, bionj, ensembl, phyloutil, env
from rasmus import fasta, util, genomeio, genomeutil
from rasmus import muscle, phylip, algorithms, clustalw, alignlib



options = [
    ["d:", "dna=", "dna", "<dna fasta>"],
    ["g:", "genomes=", "genomes", "<genome1>,<genome2>,..."],
    ["o:", "oldext=", "oldext", "<old extension>",
        {"default": ".align",
         "single": True}],
    ["n:", "newext=", "newext", "<new extension>",
        {"default": ".nt.align",
         "single": True}]
]


param = util.parseOptions(sys.argv, options, 
                          resthelp="<aa alignments>...", quit=True)



def main(param):
    env.addEnvPaths("DATAPATH")

    seqs = fasta.FastaDict()

    # read dna seqs
    util.tic("read DNA sequences")
    if "dna" in param:
        for f in param["dna"]:
            util.tic("read fasta '%s'" % f)
            seqs.read(env.findFile(f))
            util.toc()
    
    if "genomes" in param:
        genomes = param["genomes"][-1].split(",")
        for genome in genomes:
            util.tic("read fasta '%s'" % genomeio.codingfile(genome))
            seqs.read(env.findFile(genomeio.codingfile(genome)))
            util.toc()
    util.toc()
    
    
    # process alignments
    util.tic("process alignments")
    for alnFile in param[""]:
        if not os.path.exists(alnFile):
            util.log("skipping '%s', does not exist" % alnFile)
        
        newfile = alnFile.replace(param["oldext"], param["newext"])
        util.log(alnFile, "===>", newfile)
        aln = fasta.readFasta(alnFile)
        alnDna = alignlib.revtranslateAlign(aln, seqs)
        
        alnDna.write(newfile)
    util.toc()

main(param)

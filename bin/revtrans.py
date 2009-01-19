#!/usr/bin/env python

import sys, time, os
from rasmus import env, util
from rasmus.bio import fasta, alignlib



options = [
    ["c", "check", "check", "",
        {"single": True,
         "help": "check that the dna translates to the peptide"}],
    ["d:", "dna=", "dna", "<dna fasta>"],
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
    util.toc()
    
    
    # process alignments
    util.tic("process alignments")
    for alnFile in param["REST"]:
        if not os.path.exists(alnFile):
            util.log("skipping '%s', does not exist" % alnFile)
            continue
        
        newfile = alnFile.replace(param["oldext"], param["newext"])
        util.log(alnFile, "===>", newfile)
        aln = fasta.readFasta(alnFile)
        
        try:
            alnDna = alignlib.revtranslateAlign(aln, seqs, check=param["check"])
        except alignlib.TranslateError, e:
            print e
            print "%s  ! ===> %s" % (e.a, e.codon)
            print e.dna
            print e.aa
        else:    
            alnDna.write(newfile)
            
    util.toc()

main(param)

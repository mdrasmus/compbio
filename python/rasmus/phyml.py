import os
import sys

import phylip
import fasta
import util




def phyml(seqs, verbose=True, force = False, args=None, 
          usertree=None, seqtype="prot", saveOutput="", bootiter=1):
    phylip.validateSeq(seqs)
    cwd = phylip.createTempDir()
    
    util.tic("phyml on %d of length %d" % (len(seqs), len(seqs.values()[0])))

    # create input
    labels = phylip.fasta2phylip(file("infile", "w"), seqs)
    
    options = "y"
    
    if usertree != None:
        phylip.writeInTree("intree", usertree, labels)
        treefile = "intree"
        opttree = "n y"
    else:
        treefile = "BIONJ"
        opttree = "y y"
    
    if args == None:
        if seqtype == "dna":
            args = "infile 0 s 1 %d HKY e e 2 e %s %s" % \
                (bootiter, treefile, opttree)
        elif seqtype == "pep":
            args = "infile 1 s 1 %d JTT e 2 e %s %s" % \
                (bootiter, treefile, opttree)
        else:
            assert False, "unknown sequence type '%s'" % seqtype
        
    phylip.execPhylip("phyml %s" % args, options, verbose)
    
    # parse tree
    tree = phylip.readOutTree("infile_phyml_tree.txt", labels)
    
    # parse likelihood
    tree.data["logl"] = float(file("infile_phyml_lk.txt").read())
    
    if saveOutput != "":
        phylip.saveTempDir(cwd, saveOutput)
    else:
        phylip.cleanupTempDir(cwd)
    util.toc()
    
    return tree


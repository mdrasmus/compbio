import os
import sys

import phylip
import fasta
import util
import treelib



def phyml(seqs, verbose=True, args=None, 
          usertree=None, seqtype="pep", saveOutput="", bootiter=0,
          opttree=True, optbranches=True):
    
    phylip.validateSeq(seqs)
    cwd = phylip.createTempDir()
    
    util.tic("phyml on %d of length %d" % (len(seqs), len(seqs.values()[0])))

    # create input
    labels = phylip.fasta2phylip(file("infile", "w"), seqs)
    util.writeVector(file("labels", "w"), labels)
    
    options = "y"
    
    # only bootstrap when iterations are above 1
    if bootiter == 1:
        bootiter = 0
    
    if usertree != None:
        usertree = treelib.unroot(usertree)
        phylip.writeInTree("intree", usertree, labels)
        treefile = "intree"
    else:
        treefile = "BIONJ"
    
    
    optimize = ""
    if opttree:
        optimize += "y "
    else:
        optimize += "n "
    
    if optbranches:
        optimize += "y "
    else:
        optimize += "n "
    
    
    if args == None:
        if seqtype == "dna":
            args = "infile 0 s 1 %d HKY e e 2 e %s %s" % \
                (bootiter, treefile, optimize)
        elif seqtype == "pep":
            args = "infile 1 s 1 %d JTT e 2 e %s %s" % \
                (bootiter, treefile, optimize)
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


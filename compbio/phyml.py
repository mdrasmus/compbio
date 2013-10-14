"""

    Python wrapper around the PHYML phylogeny program

"""

import os
import sys

from rasmus import util
from rasmus import treelib

from . import phylip
from . import fasta




def phyml(seqs, verbose=True, args=None, 
          usertree=None, seqtype="pep", saveOutput="", bootiter=0,
          opttree=True, optbranches=True, nrates=4):
    
    phylip.validate_seqs(seqs)
    cwd = phylip.create_temp_dir()
    
    util.tic("phyml on %d of length %d" % (len(seqs), len(seqs.values()[0])))

    # create input
    labels = phylip.write_phylip_align(file("infile", "w"), seqs)
    util.write_list(file("labels", "w"), labels)
    
    options = "y"
    
    # only bootstrap when iterations are above 1
    if bootiter == 1:
        bootiter = 0
    
    if usertree != None:
        usertree = treelib.unroot(usertree)
        phylip.write_in_tree("intree", usertree, labels)
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
            args = "infile 0 s 1 %d HKY e e %d e %s %s" % \
                (bootiter, nrates, treefile, optimize)
        elif seqtype == "pep":
            args = "infile 1 s 1 %d JTT e %d e %s %s" % \
                (bootiter, nrates, treefile, optimize)
        else:
            assert False, "unknown sequence type '%s'" % seqtype
    
    
    phylip.exec_phylip("phyml %s" % args, options, verbose)
    
    # parse tree
    tree = phylip.read_out_tree("infile_phyml_tree.txt", labels)
    
    # parse likelihood
    tree.data["logl"] = float(file("infile_phyml_lk.txt").read())
    
    if saveOutput != "":
        phylip.save_temp_dir(cwd, saveOutput)
    else:
        phylip.cleanup_temp_dir(cwd)
    util.toc()
    
    return tree


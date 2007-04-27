# python libs
import os

# rasmus libs
from rasmus import util
from rasmus.bio import fasta
from rasmus.bio import phylip



def getDistMatrix(seqs, output=None, verbose=True, force = False, args=None):
    if args == None:
        args = "y"
    
    # only produce pair-wise distance estimation
    args = "k\nk\nk\n" + args
    

    phylip.validateSeq(seqs)
    cwd = phylip.createTempDir()
    util.tic("puzzle on %d of length %d" % (len(seqs), len(seqs.values()[0])))
    
    # create input
    labels = phylip.writePhylipAlign(file("infile", "w"), seqs)
    util.writeVector(file("labels", "w"), labels)
    
    # run phylip
    phylip.execPhylip("puzzle infile", args, verbose)
    
    util.toc()
    
    # parse output
    if output != None:
        os.rename("infile.dist", "../" + output)
        phylip.cleanupTempDir(cwd)
        return labels
    else:
        name, mat = phylip.readDistMatrix("infile.dist")
        phylip.cleanupTempDir(cwd)
        return labels, mat

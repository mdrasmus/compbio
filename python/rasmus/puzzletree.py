# python libs
import os

# rasmus libs
import fasta
import phylip
import util


def getDistMatrix(seqs, output=None, verbose=True, force = False, args=None):
    if args == None:
        args = "y"
    
    # only produce pair-wise distance estimation
    args = "k\nk\nk\n" + args
    

    phylip.validateSeq(seqs)
    cwd = phylip.createTempDir()
    util.tic("puzzle on %d of length %d" % (len(seqs), len(seqs.values()[0])))
    
    # create input
    labels = phylip.fasta2phylip(file("infile", "w"), seqs)
    
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

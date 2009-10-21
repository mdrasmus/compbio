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
    

    phylip.validate_seqs(seqs)
    cwd = phylip.create_temp_dir()
    util.tic("puzzle on %d of length %d" % (len(seqs), len(seqs.values()[0])))
    
    # create input
    labels = phylip.write_phylip_align(file("infile", "w"), seqs)
    util.write_list(file("labels", "w"), labels)
    
    # run phylip
    phylip.exec_phylip("puzzle infile", args, verbose)
    
    util.toc()
    
    # parse output
    if output != None:
        os.rename("infile.dist", "../" + output)
        phylip.cleanup_temp_dir(cwd)
        return labels
    else:
        name, mat = phylip.read_dist_matrix("infile.dist")
        phylip.cleanup_temp_dir(cwd)
        return labels, mat

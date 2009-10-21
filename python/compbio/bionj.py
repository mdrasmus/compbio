#
# python wrapper for BIONJ
#


# python imports
import os

# rasmus imports
from rasmus import util
from rasmus import treelib

# compbio imports
from . import phylip



def bionj(aln=None, labels=None, distmat=None, seqtype="pep", verbose=True):
    # make temp files
    distfile = util.tempfile(".", "bionj-in", ".dist")
    treefile = util.tempfile(".", "bionj-out", ".tree")
    
    # find distances and then NJ tree
    if distmat != None:
        phylip.write_dist_matrix(distmat, out=distfile)
        
        if labels == None:
            labels = aln.keys()
    else:
        if seqtype == "pep":
            labels = phylip.protdist(aln, distfile, verbose=verbose)
        else:
            labels = phylip.dnadist(aln, distfile, verbose=verbose)
    
    os.system("echo -n '%s\n%s' | bionj > /dev/null" % (distfile, treefile))
    tree = treelib.read_tree(treefile)
    phylip.rename_tree_with_names(tree, labels)
    
    # clean up
    os.remove(distfile)
    os.remove(treefile)
    
    return tree



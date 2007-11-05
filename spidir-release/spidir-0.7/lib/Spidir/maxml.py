
from rasmus.bio import phylip
from rasmus import treelib
from rasmus import util


def maxml(seqs, verbose=True, args=None, 
          usertree=None, seqtype="pep", saveOutput="", bootiter=0,
          opttree=True, optbranches=True):
    
    phylip.validateSeq(seqs)
    cwd = phylip.createTempDir()
    
    util.tic("maxml on %d of length %d" % (len(seqs), len(seqs.values()[0])))

    # create input
    seqs.write("infile.fa")
        
    if args == None:
        args = ""    
    
    phylip.execPhylip("maxml -o out -a infile.fa %s" % args, "", verbose)
    
    # parse tree
    tree = treelib.readTree("out.tree")
        
    if saveOutput != "":
        phylip.saveTempDir(cwd, saveOutput)
    else:
        phylip.cleanupTempDir(cwd)
    util.toc()
    
    return tree
    


# python imports
import os

# rasmus imports
from rasmus import util
from rasmus import treelib

# compbio imports
from . import fasta



def muscle(seqs, verbose = True, removetmp = True, options = ""):
    if len(seqs) < 2:
        return seqs

    # make input file for muscle
    infilename = util.tempfile(".", "muscle-in", ".fa")
    fasta.write_fasta(infilename, seqs)
    
    if not verbose:
        options += " -quiet "
    
    # run muscle
    outfilename = util.tempfile(".", "muscle-out", ".aln")
    cmd = "muscle " + options + " -in " + infilename + \
          " -out " + outfilename
    os.system(cmd)
    
    # parse output
    aln = fasta.read_fasta(outfilename)
    
    # cleanup tempfiles
    if removetmp:
        os.remove(infilename)
        os.remove(outfilename)
    
    return aln

    
def muscleFast(seqs, verbose = True, removetmp = True, options = ""):
    return muscle(seqs, verbose, removetmp, options + " -diags1 -sv -maxiters 1 ")


def buildAlignTree(seqs, verbose = True, removetmp = True, options = ""):
    if len(seqs) < 2:
        return seqs

    # make input file for muscle
    infilename = util.tempfile(".", "muscle-in", ".fa")
    fasta.write_fasta(infilename, seqs)
    
    # run muscle
    outfilename = util.tempfile(".", "muscle-out", ".aln")
    outfilename2 = util.tempfile(".", "muscle-out", ".tree")    
    cmd = "muscle " + options + " -in " + infilename + \
          " -out " + outfilename + " -tree2 " + outfilename2
    os.system(cmd)
    
    # parse output
    aln = fasta.read_fasta(outfilename)
    tree = treelib.Tree()
    tree.read_newick(outfilename2)
    
    # cleanup tempfiles
    if removetmp:
        os.remove(infilename)
        os.remove(outfilename)
        os.remove(outfilename2)
    
    return (aln, tree)

def buildAlignBigTree(seqs, verbose = True, removetmp = True, options = ""):
    if len(seqs) < 2:
        return seqs

    # make input file for muscle
    infilename = util.tempfile(".", "muscle-in", ".fa")
    fasta.write_fasta(infilename, seqs)
    
    # run muscle
    outfilename = util.tempfile(".", "muscle-out", ".aln")
    outfilename2 = util.tempfile(".", "muscle-out", ".tree")    
    cmd = "muscle -diags1 -sv -maxiters 1 " + options + " -in " + infilename + \
          " -out " + outfilename + " -tree1 " + outfilename2
    os.system(cmd)
    
    # parse output
    aln = fasta.read_fasta(outfilename)
    tree = treelib.Tree()
    tree.read_newick(outfilename2)
    
    # cleanup tempfiles
    if removetmp:
        os.remove(infilename)
        os.remove(outfilename)
        os.remove(outfilename2)
    
    return (aln, tree)    

def buildTree(seqs, verbose = True, removetmp = True, options = ""):

    # make input file for muscle
    infilename = util.tempfile(".", "muscle-in", ".fa")
    fasta.write_fasta(infilename, seqs)
    
    # run muscle
    outfilename = util.tempfile(".", "muscle-out", ".tree")
    cmd = "muscle " + options + " -in " + infilename + \
          " -cluster -tree1 " + outfilename
    
    if not verbose:
        cmd += " 2>/dev/null"
    
    os.system(cmd)
    
    tree = treelib.Tree()
    tree.read_newick(outfilename)
    
    if removetmp:
        os.remove(infilename)
        os.remove(outfilename)
    
    return tree


def buildBigTree(seqs, verbose = True, removetmp = True, options = ""):
    return buildTree(seqs, verbose, removetmp, options + 
        " -diags1 -sv -maxiters 1")    


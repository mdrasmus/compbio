"""

    CLUSTALW wrapper for python

    author: Matt Rasmussen
    date:   2/4/2007

"""



# python libs
import math
import os
import sys

# rasmus libs
from rasmus import treelib
from rasmus import util

# compbio imports
from . import fasta



# TODO: change removetmp to saveOutput
def clustalw(seqs, verbose=True, removetmp=True, options=""):
    """Align sequences 'seqs' with clustalw"""

    if len(seqs) < 2:
        return seqs

    # make input file for clustalw
    infilename = util.tempfile(".", "clustalw-in", ".fa")
    fasta.writeFasta(infilename, seqs)
    
    # run clustalw
    outfilename = util.tempfile(".", "clustalw-out", ".aln")
    cmd = "clustalw " + options + " -quicktree -infile=" + infilename + \
          " -outfile=" + outfilename
    if not verbose:
        cmd += " > /dev/null"
    os.system(cmd)
    
    # parse output
    aln = readClustalwAlign(outfilename)
    
    # cleanup tempfiles
    if removetmp:
        os.remove(infilename)
        os.remove(outfilename)
        os.remove(infilename.replace(".fa", ".dnd"))

    # convert output    
    return aln


def buildTree(seqs, verbose=True, removetmp=True, options=""):
    # make input file for clustalw
    infilename = util.tempfile(".", "clustalw-in", ".fa")
    fasta.writeFasta(infilename, seqs)
    
    # run clustalw
    outfilename = infilename.replace(".fa", ".ph")
    cmd = "clustalw " + options + " -tree -infile=" + infilename + \
          " -outfile=" + outfilename
    if not verbose:
        cmd += " > /dev/null"
    os.system(cmd)
    
    # parse output
    tree = treelib.Tree()
    tree.read_newick(outfilename)
    
    # cleanup tempfiles
    if removetmp:
        os.remove(infilename)
        os.remove(outfilename)
    
    return tree   


def clustalwProfiles(aln1, aln2, verbose=True, removetmp=True, options=""):
    # make input file for clustalw
    infilename1 = util.tempfile(".", "clustalw-in", ".fa")
    infilename2 = util.tempfile(".", "clustalw-in", ".fa")
    fasta.writeFasta(infilename1, aln1)
    fasta.writeFasta(infilename2, aln2)
    
    # run clustalw
    outfilename = util.tempfile(".", "clustalw-out", ".aln")
    cmd = "clustalw " + options + " -quicktree -profile1=" + infilename1 + \
          " -profile2=" + infilename2 + " -outfile=" + outfilename
    if not verbose:
        cmd += " > /dev/null"
    os.system(cmd)
    
    # parse output
    aln = readClustalwAlign(outfilename)
    
    # cleanup tempfiles
    if removetmp:
        os.remove(infilename1)
        os.remove(infilename2)
        os.remove(outfilename)
        try:
            os.remove(infilename1.replace(".fa", ".dnd"))
        except:
            pass
        try:        
            os.remove(infilename2.replace(".fa", ".dnd"))
        except:
            pass

    # convert output    
    keys = aln.keys()
    return aln



def readClustalwAlign(filename):
    infile = util.open_stream(filename)
    
    seqs = fasta.FastaDict()
    
    # skip first three lines
    infile.next()
    infile.next()
    infile.next()
    
    # parse remaining lines
    for line in infile:
        if line[0].isdigit() or line[0].isalpha():
            (name, seq) = line.split()[:2]
            if not name in seqs:
                seqs[name] = seq
            else:
                seqs[name] += seq
    return seqs


#=============================================================================
# Alignment stats
# TODO: either move to alignlib or get rid of it
#

def alignInfo(aln):
    score = 0
    
    entropyBefore = - len(aln) * .2 * math.log(.2, 2)
    
    for i in xrange(len(aln[0])):
        # count characters in a column
        charsums = {"A":0, "C":0, "G":0, "T":0, "-":0}
        for seq in aln:
            charsums[seq[i]] += 1
        
        # calc entropy of column
        entropyAfter = 0
        for char in "ACGT-":
            p = charsums[char] / float(len(aln))
            if p != 0:
                entropyAfter += - p * math.log(p, 2)
        
        score += entropyBefore - entropyAfter
    return score


def alignAvgInfo(aln):
    return alignInfo(aln) / len(aln[0])


def alignScore(seqs, match=1, mismatch=-1, gapopen=-5, gapext=-1):
    score = 0
    for i in xrange(len(seqs)):
        for j in xrange(0,i):
            for k in xrange(len(seqs[i])):
                if seqs[i][k] == "-" or seqs[j][k] == "-":
                    continue
                if seqs[i][k] == seqs[j][k]:
                    score += match
                else:
                    score += mismatch
        
        # find end of sequence
        size = len(seqs[i])
        while size > 0 and seqs[i][size-1] == "-":
            size -= 1
        
        # skip first gap
        k = 0 
        while k < size and seqs[i][k] == "-":
            k += 1
        
        # count gaps
        while k < size:
            if seqs[i][k] == "-":
                score += gapopen
            k += 1                
            while k < size and seqs[i][k] == "-":
                score += gapext
                k += 1
    return score


#=============================================================================
# testing
#

if __name__ == "__main__":
    pass

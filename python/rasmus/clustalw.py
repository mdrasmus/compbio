import math
import util
import os
import sys
import algorithms
import fasta

def clustalw(seqs, verbose = True, removetmp = True, options = ""):
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
    aln = parseAlignment(outfilename)
    
    # cleanup tempfiles
    if removetmp:
        os.remove(infilename)
        os.remove(outfilename)
        os.remove(infilename.replace(".fa", ".dnd"))

    # convert output    
    return aln


def buildTree(seqs, verbose = True, removetmp = True, options = ""):
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
    tree = algorithms.Tree()
    tree.readNewick(outfilename)
    
    # cleanup tempfiles
    if removetmp:
        os.remove(infilename)
        os.remove(outfilename)
    
    return tree   


def clustalw_profiles(aln1, aln2, verbose = True, removetmp = True, options = ""):
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
    aln = parseAlignment(outfilename)
    
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



def parseAlignment(filename):
    infile = file(filename)
    
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


def printAlign(aln, width = 59, out=sys.stdout, order=None):
    if order == None:
        order = aln.keys()

    # get basic info
    length = len(aln.values()[0])
    seqs = aln.values()
    
    # find identity positions
    identity = ""
    for i in xrange(length):
        chars = {}
        for j in xrange(len(aln)):
            char = seqs[j][i]
            if char != "-":
                if not char in chars:
                    chars[char] = 1
                else:
                    chars[char] += 1
        pid = max(chars.values()) / float(len(aln))
        if pid == 1:
            identity += "*"
        elif pid > .5:
            identity += "."
        else:
            identity += " "
    
    # print alignment
    for i in xrange(0, len(aln.values()[0]), width):
        for name in order:
            print >>out, "%20s %s" % (name, aln[name][i:i+width])
        print >>out, (" "*21) + identity[i:i+width]
        print >>out


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

def alignIdentity(aln):
    return alignScore(aln, 1, 0, 0, 0) / float(len(aln[0]))


def alignGroup(matching, groupid):
    hgroup = matching.homology[groupid]
    seqs = map(lambda gene: gene.sequence(), hgroup.genes)
    alnseqs = clustalw(seqs)
    
    aln = {}
    for i in xrange(len(alnseqs)):
        aln[hgroup.genes[i]] = alnseqs[i]
    return aln


def seqAlignScore(seqs):
    return alignScore(clustalw(seqs, True))


def seqIdentity(seqs):
    return alignIdentity(clustalw(seqs, False))


def geneAlignScore(genes):
    seqs = map(lambda gene: gene.sequence(), hgroup.genes)
    return seqSimilarityScore(seqs)


# testing
if __name__ == "__main__":
    seqs = [ "CTACACCGGCCCCA",
             "GGTTCTGGTGTTCAG",
             "CACTGCTCTTTAA",
             "CTGAGGAGATGAT"]
    seqs2 = [ "CTCCACCGGGGAACAA",
             "CTCCCCGGGGAACAA",
             "CTCCGGGGATAAA",
             "CTCCGGGGATAAA"]

    for i in xrange(4):
        for i in xrange(len(seqs)):
            seqs[i] += seqs[i]

    alnseqs = clustalw(seqs)
    print alignInfo(alnseqs)
    printAlign(alnseqs)

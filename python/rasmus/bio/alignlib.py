#
# alignlib.py
# Sequence and alignment functions
#


# python libs
import sys
import copy

# rasmus libs
from rasmus.bio import fasta
from rasmus.bio.seqlib import *



#=============================================================================
# Alignment functions
#

# TODO: maybe ignore dict alignments altogether?

def newAlign(aln=None):
    """Makes a new alignment object based on the given object
        
       given      return
       -----      ------
       dict       FastaDict
       other      other
    """
    
    if aln == None:
        return fasta.FastaDict()
    elif isinstance(aln, SeqDict):
        return type(aln)()
    else:
        return fasta.FastaDict()


def mapalign(aln, keyfunc=lambda x: x, valfunc=lambda x: x):
    """Maps the keys and values of an alignment
       
       very similar to util.mapdict()
    """
    
    aln2 = newAlign(aln)

    for key, val in aln.iteritems():
        aln2[keyfunc(key)] = valfunc(val)
    return aln2

            
def subalign(aln, cols):
    """Returns an alignment with a subset of the columns (cols)"""
    
    return mapalign(aln, valfunc=lambda x: "".join(util.mget(x, cols)))


def removeEmptyColumns(aln):
    """Removes any column from an alignment 'aln' that contains only gaps
    
       A new alignment is returned
    """
       
    dels = {}
    aln2 = newAlign(aln)
    
    for i in range(aln.alignlen()):
        col = {}
        for name in aln:
            col[aln[name][i]] = 1
        if len(col) == 1 and col.keys()[0] == '-':
            dels[i] = 1
    
    for name in aln:
        val = ""
        for i in range(len(aln[name])):
            if not i in dels:
                val += aln[name][i]
        aln2[name] = val
    aln2.orderNames(aln)
    
    return aln2


def removeGappedColumns(aln):
    """Removes any column form an alignment 'aln' that contains a gap
    
       A new alignment is returned
    """
    cols = zip(* aln.values())
    ind = util.find(lambda col: "-" not in col, cols)
    return subalign(aln, ind)


def calcConservationString(aln):
    """Returns a string of stars representing the conservation of an alignment"""
    
    percids = calcConservation(aln)
    
    # find identity positions
    identity = ""
    for pid in percids:
        if pid == 1:
            identity += "*"
        elif pid > .5:
            identity += "."
        else:
            identity += " "
    
    return identity


def calcConservation(aln):
    """Returns a list of percent matching in each column of an alignment"""

    length = len(aln.values()[0])
    seqs = aln.values()
    percids = []
    
    # find identity positions
    identity = ""
    for i in xrange(length):
        chars = util.histDict(util.cget(seqs, i))
        if "-" in chars: del chars["-"]
        
        if len(chars) == 0:
            percids.append(0)
        else:
            pid = max(chars.values()) / float(len(aln))
            percids.append(pid)
    return percids



def printAlign(aln, seqwidth = 59, spacing=2, extra=fasta.FastaDict(), 
               out=sys.stdout, order=None):
    """Pretty print an align"""
               
    if order == None:
        order = aln.keys()
    
    namewidth = max(map(len, order)) + spacing
    
    def mkname(name, namewidth):
        name2 = name[:namewidth]
        name2 += " " * (namewidth - len(name2))
        return name2

    identity = calcConservationString(aln)
    
    # print alignment
    for i in xrange(0, len(aln.values()[0]), seqwidth):
        # print sequences
        for name in order:
            print >>out, "%s %s" % (mkname(name, namewidth), 
                                    aln[name][i:i+seqwidth])
        
        # print extra
        for name in extra.keys():
            print >>out, "%s %s" % (mkname(name, namewidth), 
                                    extra[name][i:i+seqwidth])
        
        # print identity
        print >>out, (" "*namewidth) + " " + identity[i:i+seqwidth]
        print >>out


def revtranslateAlign(aaseqs, dnaseqs):
    """Reverse translates aminoacid alignment into DNA alignment
    
       Must supply original ungapped DNA.
    """
    
    align = newAlign(aaseqs)
    
    for name, seq in aaseqs.iteritems():
        align[name] = revtranslate(seq, dnaseqs[name])
    align.orderNames(aaseqs)
    
    return align



def substitutionType(char1, char2):
    """Determine the subsitution type of a pair of bases"""

    return SUBSTITUTION_TYPES[char1.upper() + char2.upper()]


def tsitTver(seq1, seq2):
    """Find the transition and transversion ratio of two sequences"""
    
    assert len(seq1) == len(seq2), "sequences are not same length"
    
    subs = map(substitutionType, seq1, seq2)
    counts = util.histInt(subs)
    
    return counts[SUB_TSIT] / float(counts[SUB_TVER])


def calcTransitionMatrix(seq1, seq2):
    """Produce a transition matrix based on two sequences"""

    assert len(seq1) == len(seq2), "sequences are not same length"
    seq1, seq2 = seq1.upper(), seq2.upper()
    
    c = util.histDict(map("".join, zip(seq1, seq2)))
    keys = filter(lambda x: "-" not in x, c.keys())
    total = float(sum(util.mget(c, keys)))
    c = util.mapdict(c, valfunc=lambda x: x/total)
    
    mat = [
        [" ",    "A",     "T",     "C",     "G"],
        ["A", c["AA"], c["AT"], c["AC"], c["AG"]],
        ["T", c["TA"], c["TT"], c["TC"], c["TG"]],
        ["C", c["CA"], c["CT"], c["CC"], c["CG"]],
        ["G", c["GA"], c["GT"], c["GC"], c["GG"]]
    ]
    
    return mat
    
def countTransitions(seq1, seq2, counts=None):
    """Count the substitution types between two sequences"""
    chars = "ATCGN-"
    
    if counts == None:
        counts = {}
    
        for a in chars:
            for b in chars:
                counts[a+b] = 0
    
    for a, b in zip(seq1, seq2):
        counts[a+b] += 1
    
    return counts


# TODO: add documentation
def checkAlignOverlap(aln, overlap):
    mat = aln.values()
    
    for i in range(len(mat)):
        for j in range(i+1, len(mat)):
            overlaps = 0
            
            for k in range(aln.alignlen()):
                gap1 = (mat[i][k] == '-')
                gap2 = (mat[j][k] == '-')
                
                if gap1 == gap2 == False:
                    overlaps += 1
            
            if overlaps < overlap * aln.alignlen():
                print "%s\t%s\toverlap %d (%2.0f%%) out of %d" % \
                    (aln.keys()[i], aln.keys()[j], overlaps, 
                     100 * (overlaps / float(aln.alignlen())),
                     aln.alignlen())


#=============================================================================
# Ka, Ks, four fold degeneracy
#

def markCodonPos(seq, pos=0):
    """
    return the codon position for each base in a gapped sequence

    codon
    ATG
    012

    gaps are given codon pos -1
    Ns are counted as bases
    """
    
    codons = []

    for base in seq:
        if base != "-":
            codons.append(pos)
            pos = (pos + 1) % 3
        else:
            codons.append(-1)

    return codons


def makeCodonPosAlign(aln):
    """Get the codon position of every base in an alignment"""
    
    def func(seq):
        dct = {-1: "-",
               0: "0",
               1: "1",
               2: "2"}
        return "".join(util.mget(dct, markCodonPos(seq)))
    return mapalign(aln, valfunc=func)


def findAlignedCodons(aln, ref=None):
    # throw out cols with gap in reference species
    if ref != None:
        ind = util.find(util.neqfunc("-"), aln[ref])
    else:
        ind = range(aln.alignlen())
    
    # throw out codons with non mod 3 gaps
    ind2 = []
    for i in range(0, len(ind), 3):
        bad = False
        
        for key, val in aln.iteritems():
            codon = util.mget(val, ind[i:i+3])
            if "-" in codon and \
               codon != ["-", "-", "-"]:
                bad = True
                break

        if not bad:
            ind2.extend(ind[i:i+3])

    return ind2




def filterAlignCodons(aln, ref=None):
    """filter an alignment for only aligned codons"""

    ind = findAlignCodons(aln, ref=ref)
    return subalign(aln, ind)


'''
def findAlignCodons(aln):
    """find all columns of aligned codons"""
    
    codonAln = mapalign(aln, valfunc=markCodonPos)
    cols = map(util.histDict, zip(* codonAln.values()))

    ind = []
    codon = []
    gaps = util.Dict(default=0)
    for i in range(len(cols)):

        if len(cols[i]) == 1:
            codon.append(i)
        elif len(cols[i]) == 2 and -1 in cols[i]:
            for key, val in aln.iteritems():
                if val[i] == "-":
                    gaps[key] += 1 
            codon.append(i)
        else:
            codon = []
        if len(codon) == 3:
            if len(gaps) == 0 or \
               util.unique([x % 3 for x in gaps.values()]) == [0]:
                ind.extend(codon)
            codon = []
            for key in gaps:
                gaps[key] = 0

    return ind
'''



def findFourFold(aln):
    """Returns index of all columns in alignment that are completely 
       fourfold degenerate
       
       Assumes that columns are already filtered for aligned codons
    """
    
    # create peptide alignment
    pepAln = mapalign(aln, valfunc=translate)
    
    # find peptide conservation
    pepcons = []
    for i in xrange(pepAln.alignlen()):
        # get a column from the peptide alignment
        col = [seq[i] for seq in pepAln.itervalues()]
        
        # compute the histogram of the column.
        # ignore gaps '-' and non-translated 'X'
        hist = util.histDict(col)
        if "-" in hist:
            del hist["-"]
        if "X" in hist:
            del hist["X"]
        
        # column is conserved if only one AA appears
        pepcons.append(len(hist) == 1 and "X" not in hist)
        
    
    ind = []
    
    # get peptides of 1st sequence
    pep = pepAln.values()[0]
    
    for i in range(0, len(aln.values()[0]), 3):
        # process only those columns that are conserved at the peptide level
        if pepcons[i//3]:
            degen = AA_DEGEN[pep[i//3]]
            
            for j in range(3):
                if degen[j] == 4:
                    ind.append(i+j)
    return ind



'''
def findFourFold(aln):
    """Returns index of all columns in alignment that are completely 
       fourfold degenerate
    """
    
    aln = filterAlignCodons(aln)
    pepAln = mapalign(aln, valfunc=translate)
    pep = pepAln.values()[0]
    
    # pep conservation
    pepcons = []
    for i in xrange(pepAln.alignlen()):
        col = [seq[i] for seq in pepAln.itervalues()]
        hist = util.histDict(col)
        if "-" in hist:
            del hist["-"]
        if "X" in hist:
            del hist["X"]
        pepcons.append(len(hist) == 1)
        

    ind = []

    for i in range(0, len(aln.values()[0]), 3):
        if pepcons[i//3]:
            degen = AA_DEGEN[pep[i//3]]
            
            for j in range(3):
                if degen[j] == 4:
                    ind.append(i+j)
    return ind
'''

def calcFourFoldDistMatrix(aln):
    names = aln.keys()

    mat = []
    # calc upper triangular
    for i in range(len(names)):
        mat.append([0.0] * (i+1))
        for j in range(i+1, len(names)):
            ind = findFourFold(aln.get([names[i], names[j]]))

            mismatches = 0
            for k in ind:
                if aln[names[i]][k] != aln[names[j]][k]:
                    mismatches += 1
            
            if len(ind) == 0:
                mat[-1].append(1.0)
            else:            
                mat[-1].append(mismatches / float(len(ind)))

    # make symmetric
    for j in range(len(names)):
        for i in range(j):
            mat[j][i] = mat[i][j]

    return mat


def findDegen(aln):
    """Determine the degeneracy of each column in an alignment"""

    codonInd = findAlignCodons(aln)
    aln2 = filterAlignCodons(aln)
    
    pepAln = mapalign(aln2, valfunc=translate)
    pep = pepAln.values()[0]
    identies = calcConservation(pepAln)
    
    degens = [-1] * len(aln.values()[0])
    
    for i in range(0, len(codonInd), 3):
        if pep[i/3] == "X":
            continue
        degen = AA_DEGEN[pep[i/3]]
        if identies[i/3] == 1.0:
            for j in range(3):
                degens[codonInd[i+j]] = degen[j]
                    
    return degens


def makeDegenStr(aln):
    """Returns a string containing the degeneracy for each column 
       in an alignment
    """

    degens = findDegen(aln)
    degenmap = {-1: " ",
                 0: "0",
                 1: "1",
                 2: "2",
                 3: "3",
                 4: "4"}
    
    return "".join(util.mget(degenmap, degens))
    

def printDegen(aln, **args):
    """Pretty print an alignment with its degeneracy for each column"""

    extra = fasta.FastaDict()
    extra["DEGEN"] = makeDegenStr(aln)
    
    printAlign(aln, extra=extra, **args)


#-------------------------------------------------------------------------------
# Position Specific Scoring Matrix (PSSM)
#-------------------------------------------------------------------------------

def align2pssm(aln, pseudocounts = {}):
    pssm = []
    denom = float(len(aln)) + sum(pseudocounts.values())
    
    for i in xrange(len(aln[0])):
        freqs = util.Dict(1, 0)
        for j in xrange(len(aln)):
            freqs[aln[j][i]] += 1
        
        for key in pseudocounts:
            freqs[key] += pseudocounts[key]
        
        for key in freqs:
            freqs[key] = math.log(freqs[key] / denom, 2)
        pssm.append(freqs)
    return pssm


def pssmSeq(pssm, seq):
    score = 0.0
    for i in xrange(len(seq)):
        score += pssm[i][seq[i]]
    return score



#--------------------------------------------------------------------------------
# Coordinate conversions
#
# Coordinate systems
# 
#   1. local
#       01234567
#       ATGCTGCG
# 
#   2. align
#       012222345567
#       ATG---CTG-CG
#
#   3. global
#       coordinate on chromosome on positive strand
#
# There should only be two kinds of indexing
# 1. 0-based, end exclusive (local/align coordinates)
# 2. 1-based, end inclusive (global coordinates)
#
#--------------------------------------------------------------------------------


def local2align(seq):
    """
    Returns list of indices of non-gap characters
    
    'ATG---CTG-CG' ==> [0,1,2,6,7,8,10,11]
    
    Used to go from local -> align space
    """
    
    lookup = []
    for i in xrange(len(seq)):
        if seq[i] == "-": continue
        lookup.append(i)
    return lookup


def align2local(seq):
    """
    Returns list such that 
    
    'ATG---CTG-CG' ==> [0,1,2,2,2,3,4,5,5,6,7]
    
    Used to go from align -> local space
    """

    i = -1
    lookup = []
    for c in seq:
        if c != "-":
            i += 1
        lookup.append(i)
    return lookup



def global2local(gobal_coord, start, end, strand):
    """Returns local coordinate in a global region"""

    # swap if strands disagree
    if strand == 1:
        return gobal_coord - start
    else:
        return end - gobal_coord


def local2global(local_coord, start, end, strand):
    """Return global coordinate within a region from a local coordinate"""
    
    # swap if strands disagree
    if strand == 1:
        return local_coord + start
    else:
        return end - local_coord


def global2align(global_coord, start, end, strand, alignLookup):
    local_coord = global2local(global_coord, start, end, strand)
    
    # throw exception for out of bounds
    if local_coord < 0 or \
       local_coord >= len(alignLookup):
        raise Exception("coordinate outside [start, end]")
    
    return alignLookup[local_coord]


def align2global(align_coord, start, end, strand, localLookup):
    local_coord = localLookup[align_coord]
    return local2global(local_coord, start, end, strand)





'''

# old code

def getAlignLookup(seq):
    """
    Returns list of indices of non-gap characters
    
    'ATG---CTG-CG' ==> [0,1,2,6,7,8,10,11]
    
    Used to go from local -> align space
    """
    
    lookup = []
    for i in xrange(len(seq)):
        if seq[i] == "-": continue
        lookup.append(i)
    return lookup


def getLocalLookup(seq):
    """
    Returns list such that 
    
    'ATG---CTG-CG' ==> [0,1,2,2,2,3,4,5,5,6,7]
    
    Used to go from align -> local space
    """

    i = -1
    lookup = []
    for c in seq:
        if c != "-":
            i += 1
        lookup.append(i)
    return lookup


def global2local(coord, start, end, strand):
    """Returns local coordinate in a global region"""

    # swap if strands disagree
    if strand == 1:
        return coord - start
    else:
        return end - coord


def local2global(coord, start, end, strand):
    """Return global coordinate within a region from a local coordinate"""
    
    # swap if strands disagree
    if strand == 1:
        return coord + start
    else:
        return end - coord


def global2align(coord, start, end, strand, alignLookup):
    coord = global2local(coord, start, end, strand)
    
    # maybe throw exception for out of bounds
    if coord < 0: coord = 0
    if coord >= len(alignLookup): coord = len(alignLookup) - 1
    
    return alignLookup[coord]


def align2global(coord, start, end, strand, localLookup):
    coord = localLookup[coord]
    return local2global(coord, start, end, strand)
'''



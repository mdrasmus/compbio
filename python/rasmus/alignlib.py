#
# alignlib.py
# Sequence and alignment functions
#


# python libs
import sys

# rasmus libs
import fasta
from seqlib import *


#--------------------------------------------------------------------------------
# Constants
#--------------------------------------------------------------------------------

CODON_TABLE = {
    "TTT": "F",  "CTT": "L",  "ATT": "I",  "GTT": "V",
    "TTC": "F",  "CTC": "L",  "ATC": "I",  "GTC": "V",
    "TTA": "L",  "CTA": "L",  "ATA": "I",  "GTA": "V",
    "TTG": "L",  "CTG": "L",  "ATG": "M",  "GTG": "V",
    "TCT": "S",  "CCT": "P",  "ACT": "T",  "GCT": "A",
    "TCC": "S",  "CCC": "P",  "ACC": "T",  "GCC": "A",
    "TCA": "S",  "CCA": "P",  "ACA": "T",  "GCA": "A",
    "TCG": "S",  "CCG": "P",  "ACG": "T",  "GCG": "A",
    "TAT": "Y",  "CAT": "H",  "AAT": "N",  "GAT": "D",
    "TAC": "Y",  "CAC": "H",  "AAC": "N",  "GAC": "D",
    "TAA": "*",  "CAA": "Q",  "AAA": "K",  "GAA": "E",
    "TAG": "*",  "CAG": "Q",  "AAG": "K",  "GAG": "E",
    "TGT": "C",  "CGT": "R",  "AGT": "S",  "GGT": "G",
    "TGC": "C",  "CGC": "R",  "AGC": "S",  "GGC": "G",
    "TGA": "*",  "CGA": "R",  "AGA": "R",  "GGA": "G",
    "TGG": "W",  "CGG": "R",  "AGG": "R",  "GGG": "G",
    
    "---": "-"
}


# make reverse codon table
REV_CODON_TABLE = {}
for key,val in CODON_TABLE.items():
    REV_CODON_TABLE.setdefault(val, []).append(key)


# make degenerate counts
CODON_DEGEN = {}
AA_DEGEN = {}
for aa, lst in REV_CODON_TABLE.items():
    folds = map(lambda x: len(util.unique(x)), zip(* lst))
    for codon in lst:
        AA_DEGEN[aa] = folds
        CODON_DEGEN[codon] = folds


# substitution types
SUB_NONE = 0  # none
SUB_TSIT = 1  # tranSition
SUB_TVER = 2  # transVersion
SUB_INS  = 3  # insert
SUB_DEL  = 4  # del
SUBSITUTION_TYPES = {
    "AA": SUB_NONE, "AC": SUB_TVER, "AG": SUB_TSIT, "AT": SUB_TVER,
    "CA": SUB_TVER, "CC": SUB_NONE, "CG": SUB_TVER, "CT": SUB_TSIT,
    "GA": SUB_TSIT, "GC": SUB_TVER, "GG": SUB_NONE, "GT": SUB_TVER,
    "TA": SUB_TVER, "TC": SUB_TSIT, "TG": SUB_TVER, "TT": SUB_NONE,
    
    "A-": SUB_DEL, "C-": SUB_DEL, "G-": SUB_DEL, "T-": SUB_DEL,
    "-A": SUB_INS, "-C": SUB_INS, "-G": SUB_INS, "-T": SUB_INS,
    
    "--": SUB_NONE, "NN": SUB_NONE, 
    "NA": SUB_NONE, "NC": SUB_NONE, "NT": SUB_NONE, "NG": SUB_NONE,    
    "AN": SUB_NONE, "CN": SUB_NONE, "TN": SUB_NONE, "GN": SUB_NONE,    
    "N-": SUB_NONE, "N-": SUB_NONE, "N-": SUB_NONE, "N-": SUB_NONE,    
    "-N": SUB_NONE, "-N": SUB_NONE, "-N": SUB_NONE, "-N": SUB_NONE
}


# hydrophobic / hydrophilic

def hydrophobic(aa):
    if aa in 'VILMFWC': return 2
    if aa in 'AYHTSPG': return 1
    if aa in 'RK': return 0.5
    return 0

def aa2property(aa): 
    aa2prop = {'A': 'weakly hydrophobic',
               'R': 'charged',
               'N': 'polar',
               'D': 'charged',
               'C': 'polar',
               'E': 'charged',
               'Q': 'polar',
               'G': 'turn',
               'H': 'charged',
               'I': 'hydrophobic',
               'L': 'hydrophobic',
               'K': 'polar',
               'M': 'met',
               'F': 'hydrophobic',
               'P': 'hydrophobic',
               'S': 'polar',
               'T': 'polar',
               'W': 'hydrophobic',
               'Y': 'polar',
               'V': 'hydrophobic',
               'U': 'polar',
               '*': 'stop',
               '-': 'gap'}
    return aa2prop[aa]



blosum62 = \
       {'A': {'A': 4, 'R':-1, 'N':-2, 'D':-2, 'C': 0, 'Q':-1, 'E':-1, 'G': 0, 'H':-2, 'I':-1, 'L':-1, 'K':-1,
              'M':-1, 'F':-2, 'P':-1, 'S': 1, 'T': 0, 'W':-3, 'Y':-2, 'V': 0, 'B':-2, 'Z':-1, 'X': 0, '*':-4},
        'R': {'A':-1, 'R': 5, 'N': 0, 'D':-2, 'C':-3, 'Q': 1, 'E': 0, 'G':-2, 'H': 0, 'I':-3, 'L':-2, 'K': 2,
              'M':-1, 'F':-3, 'P':-2, 'S':-1, 'T':-1, 'W':-3, 'Y':-2, 'V':-3, 'B':-1, 'Z': 0, 'X':-1, '*':-4},
        'N': {'A':-2, 'R': 0, 'N': 6, 'D': 1, 'C':-3, 'Q': 0, 'E': 0, 'G': 0, 'H': 1, 'I':-3, 'L':-3, 'K': 0,
              'M':-2, 'F':-3, 'P':-2, 'S': 1, 'T': 0, 'W':-4, 'Y':-2, 'V':-3, 'B': 3, 'Z': 0, 'X':-1, '*':-4},
        'D': {'A':-2, 'R':-2, 'N': 1, 'D': 6, 'C':-3, 'Q': 0, 'E': 2, 'G':-1, 'H':-1, 'I':-3, 'L':-4, 'K':-1,
              'M':-3, 'F':-3, 'P':-1, 'S': 0, 'T':-1, 'W':-4, 'Y':-3, 'V':-3, 'B': 4, 'Z': 1, 'X':-1, '*':-4},
        'C': {'A': 0, 'R':-3, 'N':-3, 'D':-3, 'C': 9, 'Q':-3, 'E':-4, 'G':-3, 'H':-3, 'I':-1, 'L':-1, 'K':-3,
              'M':-1, 'F':-2, 'P':-3, 'S':-1, 'T':-1, 'W':-2, 'Y':-2, 'V':-1, 'B':-3, 'Z':-3, 'X':-2, '*':-4},
        'Q': {'A':-1, 'R': 1, 'N': 0, 'D': 0, 'C':-3, 'Q': 5, 'E': 2, 'G':-2, 'H': 0, 'I':-3, 'L':-2, 'K': 1,
              'M': 0, 'F':-3, 'P':-1, 'S': 0, 'T':-1, 'W':-2, 'Y':-1, 'V':-2, 'B': 0, 'Z': 3, 'X':-1, '*':-4},
        'E': {'A':-1, 'R': 0, 'N': 0, 'D': 2, 'C':-4, 'Q': 2, 'E': 5, 'G':-2, 'H': 0, 'I':-3, 'L':-3, 'K': 1,
              'M':-2, 'F':-3, 'P':-1, 'S': 0, 'T':-1, 'W':-3, 'Y':-2, 'V':-2, 'B': 1, 'Z': 4, 'X':-1, '*':-4},
        'G': {'A': 0, 'R':-2, 'N': 0, 'D':-1, 'C':-3, 'Q':-2, 'E':-2, 'G': 6, 'H':-2, 'I':-4, 'L':-4, 'K':-2,
              'M':-3, 'F':-3, 'P':-2, 'S': 0, 'T':-2, 'W':-2, 'Y':-3, 'V':-3, 'B':-1, 'Z':-2, 'X':-1, '*':-4},
        'H': {'A':-2, 'R': 0, 'N': 1, 'D':-1, 'C':-3, 'Q': 0, 'E': 0, 'G':-2, 'H': 8, 'I':-3, 'L':-3, 'K':-1,
              'M':-2, 'F':-1, 'P':-2, 'S':-1, 'T':-2, 'W':-2, 'Y': 2, 'V':-3, 'B': 0, 'Z': 0, 'X':-1, '*':-4},
        'I': {'A':-1, 'R':-3, 'N':-3, 'D':-3, 'C':-1, 'Q':-3, 'E':-3, 'G':-4, 'H':-3, 'I': 4, 'L': 2, 'K':-3,
              'M': 1, 'F': 0, 'P':-3, 'S':-2, 'T':-1, 'W':-3, 'Y':-1, 'V': 3, 'B':-3, 'Z':-3, 'X':-1, '*':-4},
        'L': {'A':-1, 'R':-2, 'N':-3, 'D':-4, 'C':-1, 'Q':-2, 'E':-3, 'G':-4, 'H':-3, 'I': 2, 'L': 4, 'K':-2,
              'M': 2, 'F': 0, 'P':-3, 'S':-2, 'T':-1, 'W':-2, 'Y':-1, 'V': 1, 'B':-4, 'Z':-3, 'X':-1, '*':-4},
        'K': {'A':-1, 'R': 2, 'N': 0, 'D':-1, 'C':-3, 'Q': 1, 'E': 1, 'G':-2, 'H':-1, 'I':-3, 'L':-2, 'K': 5,
              'M':-1, 'F':-3, 'P':-1, 'S': 0, 'T':-1, 'W':-3, 'Y':-2, 'V':-2, 'B': 0, 'Z': 1, 'X':-1, '*':-4},
        'M': {'A':-1, 'R':-1, 'N':-2, 'D':-3, 'C':-1, 'Q': 0, 'E':-2, 'G':-3, 'H':-2, 'I': 1, 'L': 2, 'K':-1,
              'M': 5, 'F': 0, 'P':-2, 'S':-1, 'T':-1, 'W':-1, 'Y':-1, 'V': 1, 'B':-3, 'Z':-1, 'X':-1, '*':-4},
        'F': {'A':-2, 'R':-3, 'N':-3, 'D':-3, 'C':-2, 'Q':-3, 'E':-3, 'G':-3, 'H':-1, 'I': 0, 'L': 0, 'K':-3,
              'M': 0, 'F': 6, 'P':-4, 'S':-2, 'T':-2, 'W': 1, 'Y': 3, 'V':-1, 'B':-3, 'Z':-3, 'X':-1, '*':-4},
        'P': {'A':-1, 'R':-2, 'N':-2, 'D':-1, 'C':-3, 'Q':-1, 'E':-1, 'G':-2, 'H':-2, 'I':-3, 'L':-3, 'K':-1,
              'M':-2, 'F':-4, 'P': 7, 'S':-1, 'T':-1, 'W':-4, 'Y':-3, 'V':-2, 'B':-2, 'Z':-1, 'X':-2, '*':-4},
        'S': {'A': 1, 'R':-1, 'N': 1, 'D': 0, 'C':-1, 'Q': 0, 'E': 0, 'G': 0, 'H':-1, 'I':-2, 'L':-2, 'K': 0,
              'M':-1, 'F':-2, 'P':-1, 'S': 4, 'T': 1, 'W':-3, 'Y':-2, 'V':-2, 'B': 0, 'Z': 0, 'X': 0, '*':-4},
        'T': {'A': 0, 'R':-1, 'N': 0, 'D':-1, 'C':-1, 'Q':-1, 'E':-1, 'G':-2, 'H':-2, 'I':-1, 'L':-1, 'K':-1,
              'M':-1, 'F':-2, 'P':-1, 'S': 1, 'T': 5, 'W':-2, 'Y':-2, 'V': 0, 'B':-1, 'Z':-1, 'X': 0, '*':-4},
        'W': {'A':-3, 'R':-3, 'N':-4, 'D':-4, 'C':-2, 'Q':-2, 'E':-3, 'G':-2, 'H':-2, 'I':-3, 'L':-2, 'K':-3,
              'M':-1, 'F': 1, 'P':-4, 'S':-3, 'T':-2, 'W':11, 'Y': 2, 'V':-3, 'B':-4, 'Z':-3, 'X':-2, '*':-4},
        'Y': {'A':-2, 'R':-2, 'N':-2, 'D':-3, 'C':-2, 'Q':-1, 'E':-2, 'G':-3, 'H': 2, 'I':-1, 'L':-1, 'K':-2,
              'M':-1, 'F': 3, 'P':-3, 'S':-2, 'T':-2, 'W': 2, 'Y': 7, 'V':-1, 'B':-3, 'Z':-2, 'X':-1, '*':-4},
        'V': {'A': 0, 'R':-3, 'N':-3, 'D':-3, 'C':-1, 'Q':-2, 'E':-2, 'G':-3, 'H':-3, 'I': 3, 'L': 1, 'K':-2,
              'M': 1, 'F':-1, 'P':-2, 'S':-2, 'T': 0, 'W':-3, 'Y':-1, 'V': 4, 'B':-3, 'Z':-2, 'X':-1, '*':-4},
        'B': {'A':-2, 'R':-1, 'N': 3, 'D': 4, 'C':-3, 'Q': 0, 'E': 1, 'G':-1, 'H': 0, 'I':-3, 'L':-4, 'K': 0,
              'M':-3, 'F':-3, 'P':-2, 'S': 0, 'T':-1, 'W':-4, 'Y':-3, 'V':-3, 'B': 4, 'Z': 1, 'X':-1, '*':-4},
        'Z': {'A':-1, 'R': 0, 'N': 0, 'D': 1, 'C':-3, 'Q': 3, 'E': 4, 'G':-2, 'H': 0, 'I':-3, 'L':-3, 'K': 1,
              'M':-1, 'F':-3, 'P':-1, 'S': 0, 'T':-1, 'W':-3, 'Y':-2, 'V':-2, 'B': 1, 'Z': 4, 'X':-1, '*':-4},
        'X': {'A': 0, 'R':-1, 'N':-1, 'D':-1, 'C':-2, 'Q':-1, 'E':-1, 'G':-1, 'H':-1, 'I':-1, 'L':-1, 'K':-1,
              'M':-1, 'F':-1, 'P':-2, 'S': 0, 'T': 0, 'W':-2, 'Y':-1, 'V':-1, 'B':-1, 'Z':-1, 'X':-1, '*':-4},
        '*': {'A':-4, 'R':-4, 'N':-4, 'D':-4, 'C':-4, 'Q':-4, 'E':-4, 'G':-4, 'H':-4, 'I':-4, 'L':-4, 'K':-4,
              'M':-4, 'F':-4, 'P':-4, 'S':-4, 'T':-4, 'W':-4, 'Y':-4, 'V':-4, 'B':-4, 'Z':-4, 'X':-4, '*': 1}}


#--------------------------------------------------------------------------------
# Sequence functions
#--------------------------------------------------------------------------------

def translate(dna):
    """Translates DNA (with gaps) into amino-acids"""
    
    aa = []
    
    assert len(dna) % 3 == 0, "dna sequence length is not a multiple of 3"
    
    for i in xrange(0, len(dna), 3):
        if "N" in dna[i:i+3]:
            aa.append("X")     # unkown aa
        else:
            aa.append(CODON_TABLE[dna[i:i+3]])
    return "".join(aa)


def revtranslate(aa, dna):
    """Reverse translates aminoacids (with gaps) into DNA
    
       Must supply original ungapped DNA.
    """

    seq = []
    i = 0
    for a in aa:
        if a == "-":
            seq.append("---")
        else:
            seq.append(dna[i:i+3])
            i += 3
    return "".join(seq)


def revcomp(seq):
    """Reverse complement a sequence"""

    comp = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}
    
    seq2 = []
    for i in range(len(seq)-1, -1, -1):
        seq2.append(comp[seq[i]])
    return "".join(seq2)


def cgContent(seq):
    hist = util.histDict(seq)
    total = hist["A"] + hist["C"] + hist["T"] + hist["G"]
    
    return (hist["C"] + hist["G"]) / float(total)


#--------------------------------------------------------------------------------
# Alignment functions
#--------------------------------------------------------------------------------

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
    aln2 = newAlign(aln)

    for key, val in aln.iteritems():
        aln2[keyfunc(key)] = valfunc(val)
    return aln2

            
def subalign(aln, cols):
    return mapalign(aln, valfunc=lambda x: "".join(util.mget(x, cols)))


def removeEmptyColumns(seqs):
    dels = {}
    seqs2 = newAlign(seqs)
    
    for i in range(len(seqs.values()[0])):
        col = {}
        for name in seqs:
            col[seqs[name][i]] = 1
        if len(col) == 1 and col.keys()[0] == '-':
            dels[i] = 1
    
    for name in seqs:
        val = ""
        for i in range(len(seqs[name])):
            if not i in dels:
                val += seqs[name][i]
        seqs2[name] = val
    seqs2.orderNames(seqs)
    
    return seqs2


def calcConservationString(aln):
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



def subType(char1, char2):
    return SUBSITUTION_TYPES[char1 + char2]


def tsitTver(seq1, seq2):
    """Find the transition and transversion ratio of two sequences"""
    
    assert len(seq1) == len(seq2), "sequences are not same length"
    
    subs = map(subType, seq1, seq2)
    counts = util.histInt(subs)
    
    return counts[SUB_TSIT] / float(counts[SUB_TVER])


def transitionMatrix(seq1, seq2):
    assert len(seq1) == len(seq2), "sequences are not same length"
    
    c = util.histDict(map("".join, zip(seq1, seq2)))
    keys = filter(lambda x: "-" not in x, c.keys())
    total = float(sum(util.sublist(c, keys)))
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
    chars = "ATCGN-"
    
    if counts == None:
        counts = {}
    
        for a in chars:
            for b in chars:
                counts[a+b] = 0
    
    for a, b in zip(seq1, seq2):
        counts[a+b] += 1
    
    return counts


#-------------------------------------------------------------------------------
# Ka, Ks, four fold degeneracy
#-------------------------------------------------------------------------------


def markCodonPos(seq, pos=0):
    """
    return the codon position for each base in a gapped sequence

    codon
    ATG
    012

    gaps are given codon pos -1
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


def findAlignCodons(aln):
    """find all columns of aligned codons"""
    
    codonAln = mapalign(aln, valfunc=markCodonPos)
    cols = map(util.histDict, zip(* codonAln.values()))

    ind = []
    codon = []
    for i in range(len(cols)):

        if len(cols[i]) == 1 or \
           len(cols[i]) == 2 and -1 in cols[i]:
            codon.append(i)
        else:
            codon = []
        if len(codon) == 3:
            ind.extend(codon)
            codon = []

    return ind


def filterAlignCodons(aln):
    """filter an alignment for only aligned codons"""

    ind = findAlignCodons(aln)
    return mapalign(aln, valfunc=lambda x: "".join(util.mget(x, ind)))


def findFourFold(aln):
    """Returns index of all columns in alignment that are completely 
       fourfold degenerate"""
    
    aln = filterAlignCodons(aln)
    pepAln = mapalign(aln, valfunc=translate)
    pep = pepAln.values()[0]
    identies = calcConservation(pepAln)

    ind = []

    for i in range(0, len(aln.values()[0]), 3):
        if pep[i/3] == "X":
            continue
        degen = AA_DEGEN[pep[i/3]]
        if identies[i/3] == 1.0:
            for j in range(3):
                if degen[j] == 4:
                    ind.append(i+j)
    return ind


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
    degens = findDegen(aln)
    degenmap = {-1: " ",
                 0: "0",
                 1: "1",
                 2: "2",
                 3: "3",
                 4: "4"}
    
    return "".join(util.mget(degenmap, degens))
    

def printDegen(aln, **args):
    extra = fasta.FastaDict()
    extra["DEGEN"] = makeDegenStr(aln)
    
    printAlign(aln, extra=extra, **args)


#-------------------------------------------------------------------------------
# Position Specific Scoring Matrix (PSSM)
#-------------------------------------------------------------------------------

def align2pssm(seqs, pseudocounts = {}):
    pssm = []
    denom = float(len(seqs)) + sum(pseudocounts.values())
    
    for i in xrange(len(seqs[0])):
        freqs = util.Dict(1, 0)
        for j in xrange(len(seqs)):
            freqs[seqs[j][i]] += 1
        
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
#       02134567
#       ATGCTGCG
# 
#   2. align
#       012222345567
#       ATG---CTG-CG
#
#   3. global
#       coordinate on chromosome on positive strand
#
# Notes: end is inclusive
#
#--------------------------------------------------------------------------------

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

    i = 0
    lookup = []
    for c in seq:
        lookup.append(i)
        if c != "-":
            i += 1
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
    
    if coord < 0: coord = 0
    if coord >= len(alignLookup): coord = len(alignLookup) - 1
    
    return alignLookup[coord]


def align2global(coord, start, end, strand, localLookup):
    coord = localLookup[coord]
    return local2global(coord, start, end, strand)




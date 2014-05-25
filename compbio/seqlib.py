
# python imports
import copy
import math
import random

# rasmus imports
from rasmus import util


class SeqDict (dict):
    """
    A dictionary for molecular sequences.  Also keeps track of their order,
    useful for reading and writing sequences from fasta's.  See fasta.FastaDict
    for subclass that implements FASTA reading and writing.
    """

    def __init__(self):
        dict.__init__(self)

        self.names = []

    def order_names(self, aln):
        """Orders the names in the same order they appear in aln"""

        lookup = util.list2lookup(aln.keys())
        self.names.sort(key=lambda x: lookup[x])

    # add a key, value pair
    def add(self, key, value, errors=False):
        if key in self:
            if errors:
                util.logger("duplicate key", key)

            # keep the longest value, by default
            if len(value) >= len(self[key]):
                dict.__setitem__(self, key, value)
        else:
            self.names.append(key)
            dict.__setitem__(self, key, value)

    def get(self, keys, new=None):
        """Return a subset of the sequences"""

        if new is None:
            new = type(self)()

        for key in keys:
            if key in self:
                new[key] = self[key]

        return new

    def alignlen(self):
        """
        If this SeqDict is an alignment, this function
        will return its length
        """
        return len(self.values()[0])

    # The following methods keep names in sync with dictionary keys
    def __setitem__(self, key, value):
        if key not in self:
            self.names.append(key)
        dict.__setitem__(self, key, value)

    def __delitem__(self, key):
        self.names.remove(key)

    def update(self, dct):
        for key in dct:
            if key not in self.names:
                self.names.append(key)
        dict.update(self, dct)

    def setdefault(self, key, value):
        if key not in self.names:
            self.names.append(key)
        dict.setdefault(self, key, value)

    def clear(self):
        self.names = []
        dict.clear(self)

    # keys are always sorted in order added
    def keys(self):
        return list(self.names)

    def iterkeys(self):
        return iter(self.names)

    def values(self):
        return [self[key] for key in self.iterkeys()]

    def itervalues(self):
        def func():
            for key in self.iterkeys():
                yield self[key]
        return func()

    def iteritems(self):
        def func():
            for key in self.iterkeys():
                yield (key, self[key])
        return func()

    def __iter__(self):
        return iter(self.names)

    def __len__(self):
        return len(self.names)


#-----------------------------------------------------------------------------
# Constants
#-----------------------------------------------------------------------------


# standard codon table
CODON_TABLE = {
    "TTT": "F", "CTT": "L", "ATT": "I", "GTT": "V",
    "TTC": "F", "CTC": "L", "ATC": "I", "GTC": "V",
    "TTA": "L", "CTA": "L", "ATA": "I", "GTA": "V",
    "TTG": "L", "CTG": "L", "ATG": "M", "GTG": "V",

    "TCT": "S", "CCT": "P", "ACT": "T", "GCT": "A",
    "TCC": "S", "CCC": "P", "ACC": "T", "GCC": "A",
    "TCA": "S", "CCA": "P", "ACA": "T", "GCA": "A",
    "TCG": "S", "CCG": "P", "ACG": "T", "GCG": "A",

    "TAT": "Y", "CAT": "H", "AAT": "N", "GAT": "D",
    "TAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D",
    "TAA": "*", "CAA": "Q", "AAA": "K", "GAA": "E",
    "TAG": "*", "CAG": "Q", "AAG": "K", "GAG": "E",

    "TGT": "C", "CGT": "R", "AGT": "S", "GGT": "G",
    "TGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
    "TGA": "*", "CGA": "R", "AGA": "R", "GGA": "G",
    "TGG": "W", "CGG": "R", "AGG": "R", "GGG": "G",

    "---": "-"
}

# codon table specific to the Candida species
CANDIDA_CODON_TABLE = copy.copy(CODON_TABLE)
CANDIDA_CODON_TABLE["CTG"] = "S"  # originally L


# make reverse codon table
REV_CODON_TABLE = {}
for key, val in CODON_TABLE.items():
    REV_CODON_TABLE.setdefault(val, []).append(key)


# make degenerate counts
#
# example:
#
# CGT => "R"
# CGC => "R"
# CGA => "R"
# CGG => "R"
#
# CODON_DEGEN["R"] = [1, 1, 4]
# CODON_DEGEN["CGT"] = [1, 1, 4]
#
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
SUB_INS = 3  # insert
SUB_DEL = 4  # del
SUBSTITUTION_TYPES = {
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
    if aa in 'VILMFWC':
        return 2.0
    elif aa in 'AYHTSPG':
        return 1.0
    elif aa in 'RK':
        return 0.5
    else:
        return 0.0


AA_PROPERTY = {'A': 'weakly hydrophobic',
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


BASE2INT = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3
}

INT2BASE = ["A", "C", "G", "T"]


#=============================================================================
# Sequence functions
#

class TranslateError (Exception):
    def __init__(self, msg, aa=None, dna=None, a=None, codon=None):
        Exception.__init__(self, msg)
        self.aa = aa
        self.dna = dna
        self.a = a
        self.codon = codon


def translate(dna, table=CODON_TABLE):
    """Translates DNA (with gaps) into amino-acids"""

    aa = []

    assert len(dna) % 3 == 0, "dna sequence length is not a multiple of 3"

    for i in xrange(0, len(dna), 3):
        codon = dna[i:i+3].upper()
        if "N" in codon:
            aa.append("X")     # unkown aa
        else:
            aa.append(table[codon])
    return "".join(aa)


def revtranslate(aa, dna, check=False):
    """Reverse translates aminoacids (with gaps) into DNA

       Must supply original ungapped DNA.
    """

    # trim stop codon
    #if aa[-1] in "*X" and CODON_TABLE.get(dna[-3:], "") == "*":
    #    aa = aa[:-3]
    #    dna = dna[:-3]

    a = len(aa.replace("-", "")) * 3
    b = len(dna.replace("-", ""))

    if a != b:
        raise TranslateError(
            "sequences have wrong lengths (pep %d != dna %d)" %
            (a, b), aa, dna, None, None)

    seq = []
    i = 0
    for a in aa:
        if a == "-":
            seq.append("---")
        else:
            codon = dna[i:i+3]
            if check and a != CODON_TABLE.get(codon, "X"):
                raise TranslateError("bad translate", aa, dna, a, codon)
            seq.append(codon)
            i += 3
    return "".join(seq)


_comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N",
         "a": "t", "c": "g", "g": "c", "t": "a", "n": "n",
         "R": "Y", "Y": "R", "S": "W", "W": "S", "K": "M", "M": "K",
         "r": "y", "y": "r", "s": "w", "w": "s", "k": "m", "m": "k",
         "B": "V", "V": "B", "D": "H", "H": "D",
         "b": "v", "v": "b", "d": "h", "h": "d"}


def revcomp(seq):
    """Reverse complement a sequence"""

    seq2 = []
    for i in xrange(len(seq)-1, -1, -1):
        seq2.append(_comp[seq[i]])
    return "".join(seq2)


def gcContent(seq):
    hist = util.hist_dict(seq)
    total = hist["A"] + hist["C"] + hist["T"] + hist["G"]

    return (hist["C"] + hist["G"]) / float(total)


#=============================================================================
# Kimura sequence mutation model
#
# TODO: maybe move to phylo module


KIMURA_MATRIX = [
    ['r', 's', 'u', 's'],
    ['s', 'r', 's', 'u'],
    ['u', 's', 'r', 's'],
    ['s', 'u', 's', 'r']
]


def evolveKimuraSeq(seq, time, alpha=1, beta=1):
    probs = {
        's': .25 * (1 - math.e**(-4 * beta * time)),
        'u': .25 * (1 + math.e**(-4 * beta * time)
                    - 2*math.e**(-2*(alpha+beta)*time))
    }
    probs['r'] = 1 - 2*probs['s'] - probs['u']

    seq2 = []

    for base in seq:
        cdf = 0
        row = KIMURA_MATRIX[BASE2INT[base]]
        pick = random.random()

        for i in range(4):
            cdf += probs[row[i]]
            if cdf >= pick:
                seq2.append(INT2BASE[i])
                break

    assert len(seq2) == len(seq), "probabilities do not add to one"

    return "".join(seq2)


def evolveKimuraBase(base, time, alpha, beta):
    probs = {
        's': .25 * (1 - math.e**(-4 * beta * time)),
        'u': .25 * (1 + math.e**(-4 * beta * time)
                    - 2*math.e**(-2*(alpha+beta)*time))
    }
    probs['r'] = 1 - 2*probs['s'] - probs['u']

    cdf = 0
    row = KIMURA_MATRIX[BASE2INT[base]]
    pick = random.random()

    for i in range(4):
        cdf += probs[row[i]]
        if cdf >= pick:
            return INT2BASE[i]

    assert False, "probabilities do not add to one"

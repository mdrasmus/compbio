import sys
import os

from rasmus import util
from rasmus.seqlib import SeqDict
#from alignlib import printAlign



def removestar(value):
    return value.replace("*", "")

def firstword(key):
    return key.split()[0]


class FastaDict (SeqDict):
    """Store a FASTA file as a dictionary-like object
       
       FastaDict works exactly like a python dict except keys are guaranteed to
       in the same order as they appear in the file.  If a FASTA has *.psd and
       *.psi index files direct indexing with 'fastacmd' will be used.  This 
       feature can be disabled with useIndex=False when reading a FASTA file.
    """
       

    def __init__(self, *args, ** keywords):
        SeqDict.__init__(self)
        
        self.filenames = []
        
        if len(args) > 0:
            self.read(* args, **keywords)
    
    
    def read(self, filename, keyfunc=firstword, valuefunc = lambda x: x, 
              errors=True, useIndex=True):
        key = ""
        value = ""
        
        if isinstance(filename, str) and useIndex and hasFastaIndex(filename):
            self.filenames.append(filename)
            
            # store None's for when indexing should be used
            infile = os.popen("grep '>' '%s'" % filename)
            for line in infile:
                key = line[1:].rstrip()
                self.names.append(key)
        else:
            for line in util.openStream(filename):
                if line[0] == ">":
                    if key != "":
                        self.add(key, valuefunc(value), errors)
                    key = keyfunc(line[1:].rstrip())
                    value = ""
                else:
                    assert key != ""
                    value += line.rstrip()
            if key != "":
                self.add(key, valuefunc(value), errors)
    
    
    def write(self, filename=sys.stdout, names = None, width=80):
        out = util.openStream(filename, "w")
        
        if names == None:
            names = self.names
        
        for key in names:
            print >>out, ">"+ key
            util.printwrap(self[key], width, out=out)
    
    
    def __getitem__(self, key):
        if not SeqDict.__contains__(self, key):
            # if val == None, then we are using fasta indexing
            for filename in self.filenames:
                try:
                    return fastaGet(filename, key)
                except: pass
            
        else:
            return SeqDict.__getitem__(self, key)
    
    def __contains__(self, key):
        if not SeqDict.__contains__(self, key):
            # if val == None, then we are using fasta indexing
            for filename in self.filenames:
                try:
                    fastaGet(filename, key)
                    return True
                except: pass
            return False
        else:
            return True
    
    
    def alignlen(self):
        """If this fasta is an alignment, this function will return its length"""
        
        return len(self.values()[0])


def readFasta(filename, keyfunc=firstword, valuefunc = lambda x: x, 
              errors=True):   
    """Read a FASTA file into a sequence dictionary"""
    
    fa = FastaDict()
    fa.read(filename, keyfunc, valuefunc, errors)
    return fa



def readFastaOrdered(filename, keyfunc=firstword, valuefunc=lambda x:x):
    """Read a FASTA file into a 'keys' and 'values' lists"""
    
    infile = util.openStream(filename)
    seqs = []
    names = []
    key = ""
    value = ""
        
    for line in infile:
        if line[0] == ">":
            if key != "":
                names.append(keyfunc(key))
                seqs.append(valuefunc(value))
            key = line[1:].rstrip()
            value = ""
        else:
            assert key != ""
            value += line.rstrip()
    if key != "":
        names.append(keyfunc(key))
        seqs.append(valuefunc(value))
    return (names, seqs)



def _revcomp(seq):
    """Reverse complement a sequence"""

    comp = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N", 
            "a":"t", "c":"g", "g":"c", "t":"a", "n":"n"}
    
    seq2 = []
    for i in xrange(len(seq)-1, -1, -1):
        seq2.append(comp[seq[i]])
    return "".join(seq2)



def fastaGet(fastaFile, key, start=0, end=0, strand=1):
    """Get a sequence from a fasta file that has been indexed by 'formatdb'"""
    
    stream = os.popen("fastacmd -d %s -s %s -L %d,%d 2>/dev/null" % 
                      (fastaFile, key, start, end))
    
    # remove key
    val = stream.read()
    if val == "":
        raise Exception("no such sequence")
    else:
        seq = val.split("\n")[1:]
        seq = "".join(seq)
    
    if strand == -1:
        seq = _revcomp(seq)
    
    return seq


def hasFastaIndex(fastaFile):
    """Check to see if fastaFile has a formatdb fasta index"""

    return os.path.exists(fastaFile + ".psd") and \
           os.path.exists(fastaFile + ".psi")



def writeFasta(filename, seqs, order = None, width=None):
    """Write a FASTA dictionary into a file"""
    
    out = util.openStream(filename, "w")
    
    if type(seqs) == list:
        names = map(str, range(len(seqs)))
        writeFastaOrdered(out, names, seqs, width)
    else:
        if order == None:
            order = seqs.keys()
            
        for key in order:
            print >>out, ">"+ key
            util.printwrap(seqs[key], width, out=out)


def writeFastaOrdered(filename, names, seqs, width=None):
    out = util.openStream(filename, "w")
    
    for name,seq in zip(names,seqs):
        print >>out, ">%s" % name
        util.printwrap(seq, width, out=out)


#def getkey(key, field, delim="|"):
#    return key.split(delim)[field]


def array2dict(names, seqs):
    fa = FastaDict()
    for name, seq in zip(names, seqs):
        fa.add(name, seq)
    return fa


def dict2array(fa, order = None):
    names = []
    seqs = []
    
    if order == None:
        order = fa.keys()
    
    for name in order:
        names.append(name)
        seqs.append(fa[name])
    return (names, seqs)



#
# Special alignment format
#


class Alignment:
    def __init__(self, key, value):
        tokens = key.split("|")
        
        self.genome = tokens[0]        
        self.chrom, self.start, self.end = chromParse(tokens[1])
        self.kind   = tokens[2]
        self.seq    = value
        
        if self.kind == "_aligned":
            self.strand = 1
        elif self.kind == " revcomp_aligned":
            self.strand = -1
        else:
            raise Exception("unknown kind")
    
    def makekey(self):
        return "%s|%s|%s" % \
            (self.genome, 
             chromFormat(self.chrom, self.start, self.end), 
             self.kind)


def chromParse(string):
    tokens = string.split(":")
    chrom  = tokens[0].replace("chr", "")
    
    tokens = tokens[1].split("-")
    start  = util.pretty2int(tokens[0])
    end    = util.pretty2int(tokens[1])
    
    return (chrom, start, end)


def chromFormat(chrom, start, end):
    if chrom.isdigit():
        prefix = "chr"
    else:
        prefix = ""
    return "%s%s:%s-%s" % \
        (prefix, chrom, util.int2pretty(start), util.int2pretty(end))

    
def readAlignment(filename):
    aln = []
    names, seqs = readFastaOrdered(filename)
    for name, seq in zip(names, seqs):
        aln.append(Alignment(name, seq))
    return aln


def writeAlignment(out, aln):
    writeFasta(out, alignment2fasta(aln), order = keys)


def alignment2fasta(aln):
    fa = {}
    keys = []
    for seq in aln:
        key = seq.makekey()
        keys.append(key)
        fa[key] = seq.seq
    return fa    




def makeFastaIndex(filename):
    infile = util.openStream(filename)
    
    index = {}
    
    for line in util.SafeReadIter(infile):
        if line.startswith(">"):
            index[line[1:].rstrip()] = infile.tell()
    
    return index


class FastaIndex:
    def __init__(self, filename):
        # open fasta
        self.fafile = file(filename, "rb")
        
        # estimate column width
        numlines = 5
        widths = []
        for line in self.fafile:
            if len(line) != 0  and line[0] != ">":
                widths.append(len(line.rstrip()))
                if len(widths) == numlines:
                    break
        if not util.equal(* widths):
            raise Exception("lines do not have consistent width")
        self.width = widths[0]
        
        # read index
        self.index = {}
        for key, ind in util.DelimReader(filename + ".index"):
            self.index[key] = int(ind)
    
    
    def get(self, key, start=1, end=None, strand=1):
        assert start > 0, Exception("must specify coordinates one-based")
        assert key in self.index, Exception("key '%s' not in index" % key)
        
        # must translate from one-based to zero-based
        # must account for newlines
        start2 = self.index[key] + start - 1 + ((start-1) / self.width)
        
        # seek to beginning
        self.fafile.seek(start2)
        
        # read until end
        if end != None:
            end2 = self.index[key] + end + ((end-1) / self.width)
            readlen = end2 - start2
            seq = self.fafile.read(readlen).replace("\n", "")
        else:
            # read until end of sequence
            seq = []
            while True:
                line = self.fafile.readline()
                if line.startswith(">") or len(line) == 0:
                    break
                seq.append(line.rstrip())
            seq = "".join(seq)
        
        
        # reverse complement if needed
        if strand == -1:
            seq = _revcomp(seq)
        
        return seq



"""


def makeFastaIndex(filename, keyfunc = lambda x: x):
    def addKey(index, key, locs):
        if key in index:
            raise "duplicate key '%s'" % key
        index[keyfunc(key)] = locs

    infile = util.openStream(filename)
    index = {}
    key = ""
    locs = []
    ind = 0
    
    for line in util.SafeReadIter(infile):
        if line[0] == ">":
            if key != "":
                addKey(index, key, locs)
            key = line[1:].rstrip()
            locs = []
            ind = 0
        else:
            assert key != ""
            locs.append((ind, pos))
            ind += len(line.rstrip())
        pos = infile.tell()
    if key != "":
        addKey(index, key, locs)
    
    return index




def readFastaIndex(filename, keyfunc=firstword):
    infile = util.openStream(filename)
    index = {}
    
    for line in infile:
        words = line.split("\t")
        key = words[0]
        locs = []
        
        for i in xrange(1, len(words), 2):
            locs.append((int(words[i]), int(words[i+1])))
        index[keyfunc(key)] = locs
    
    return index
        

def writeFastaIndex(filename, index, order=None):
    out = util.openStream(filename, "w")
    
    if order == None:
        order = index.keys()
        order.sort()
    
    for key in order:
        out.write(key)
        for ind, pos in index[key]:
            out.write("\t%d\t%d" % (ind, pos))
        out.write("\n")
            
#
# end is inclusive
#
def indexSeq(fastaFile, index, key, start=1, end=None):
    infile = util.openStream(fastaFile)
    locs = index[key]
    
    # convert start,end to from 1-based to 0-based
    start -= 1
    end -= 1
    
    
    # find nearest starting index
    i = 0
    while i < len(locs) and locs[i][0] < start: i += 1
    
    # determine desired sequence length
    if end != None:
        seqlen = end - start + 1
    else:
        seqlen = 1e1000 
    
    infile.seek(locs[i][1] + start - locs[i][0])
    seq = ""
    for line in util.SafeReadIter(infile):
        if line[0] == ">":
            break
        else:
            seq += line.rstrip()
            if len(seq) > seqlen:
                seq = seq[:seqlen]
                break
    return seq
"""



"""
def getCodonPos(seq):
    pos = []
    i = 0
    for j in seq:
        if j != "-":
            pos.append(i)
            if i < 2:
                i += 1
            else:
                i = 0
    return pos



# get Ks, Ka stats

def countCodonSub(seq1, seq2):
    assert len(seq1) == len(seq2), "sequences are not aligned"

    # synonymous and non-synonymous subsitutions
    ks = 0
    ka = 0
    
    pos1 = getCodonPos(seq1)
    pos2 = getCodonPos(seq2)
    for i in range(len(pos1)-3):
        # only process aligned codons
        if pos1[i:i+3] == [0,1,2] and \
           pos2[i:i+3] == [0,1,2]:
            # determine if codons are synonymous and count mismatches
            if dna2aa(seq1[i:i+3]) == dna2aa(seq2[i:i+3]):
                for j in range(i,i+3):
                    if seq1[j] != seq2[j]:
                        ks += 1
            else:
                for j in range(i,i+3):
                    if seq1[j] != seq2[j]:
                        ka += 1
    
    return {"ks": ks, "ka": ka}
"""

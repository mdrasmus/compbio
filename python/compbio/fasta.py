"""

    Module for working with FASTA files

"""

# python imports
import sys
import os
from itertools import izip

# rasmus imports
from rasmus import util

# seqlib imports
from . import seqlib
from seqlib import SeqDict



def removestar(value):
    return value.replace("*", "")

def firstword(key):
    """Use the first word as the key"""
    return key.split()[0]


class FastaDict (SeqDict):
    """Store a FASTA file as a dictionary-like object
       
       FastaDict works exactly like a python dict except keys are guaranteed to
       in the same order as they appear in the file.  
    """
       

    def __init__(self, *args, ** keywords):
        SeqDict.__init__(self)
        
        self.index = FastaIndex()
        
        if len(args) > 0:
            self.read(* args, **keywords)
    
    
    def read(self, filename, keyfunc=firstword, valuefunc = lambda x: x, 
              errors=True, useIndex=False):
        """Read sequences from a Fasta file"""
        
        if isinstance(filename, str) and useIndex and has_fasta_index(filename):
            newkeys = self.index.read(filename)
            
            # store None's for when indexing should be used
            for key in newkeys:
                if key not in self:
                    self.names.append(key)
                dict.__setitem__(self, key, None)
        else:
            for key, seq in iter_fasta(filename, keyfunc, valuefunc):
                self.add(key, seq, errors)
    
    
    def write(self, filename=sys.stdout, names=None, width=80):
        """Write sequences in Fasta format"""
        
        out = util.open_stream(filename, "w")
        
        if names is None:
            names = self.names
        
        for key in names:
            print >>out, ">" + key
            util.printwrap(self[key], width, out=out)
    
    
    def __getitem__(self, key):
        """Get a sequence by key"""
        
        val = SeqDict.__getitem__(self, key) 
        
        if val is None:
            # if val == None, then we are using fasta indexing
            val = self.index.get(key)
                
            # cache value
            self[key] = val
            return val
        else:
            return val
    
    
    def getseq(self, key, start=1, end=None, strand=1):
        """Get a sequence (or subsequence) by key"""
        
        val = SeqDict.__getitem__(self, key) 
        
        if val is None:
            # if val == None, then we are using fasta indexing
            return self.index.get(key, start, end, strand)

        else:
            start = util.clamp(start, 1, None)
            end = util.clamp(end, 1, None)
            val = val[start-1:end]
        
            # reverse complement if needed
            if strand == -1:
                val = _revcomp(val)
            
            return val


#=============================================================================
# Convenience functions for input/output
#
    

def read_fasta(filename, keyfunc=firstword, valuefunc = lambda x: x, 
              errors=True, useIndex=True):   
    """Read a FASTA file into a sequence dictionary"""
    
    fa = FastaDict()
    fa.read(filename, keyfunc, valuefunc, errors, useIndex=useIndex)
    return fa
readFasta = read_fasta


def write_fasta(filename, seqs, order = None, width=None):
    """Write a FASTA dictionary into a file"""
    
    out = util.open_stream(filename, "w")
    
    if type(seqs) == list:
        names = map(str, range(len(seqs)))
        writeFastaOrdered(out, names, seqs, width)
    else:
        seqs.write(filename, order, width)
writeFasta = write_fasta


def iter_fasta(filename, keyfunc=firstword, valuefunc = lambda x: x):
    """Iterate through the sequences of a FASTA file"""
    key = ""
    value = ""
    
    for line in util.open_stream(filename):
        if len(line) > 0 and line[0] == ">":
            if key != "":
                yield (key, valuefunc(value))
            key = keyfunc(line[1:].rstrip())
            value = ""
        else:
            assert key != ""
            value += line.rstrip()
    if key != "":
        yield (key, valuefunc(value))


_comp = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N", 
         "a":"t", "c":"g", "g":"c", "t":"a", "n":"n",
         "R":"Y", "Y":"R", "S":"W", "W":"S", "K":"M", "M":"K",
         "r":"y", "y":"r", "s":"w", "w":"s", "k":"m", "m":"k",
         "B":"V", "V":"B", "D":"H", "H":"D",
         "b":"v", "v":"b", "d":"h", "h":"d"}

def _revcomp(seq):
    """Reverse complement a sequence"""
    
    seq2 = []
    for i in xrange(len(seq)-1, -1, -1):
        seq2.append(_comp[seq[i]])
    return "".join(seq2)



#=============================================================================
# Rasmus FASTA Indexing
#

def makeFastaIndex(filename):
    """I also have a faster C program called formatfa"""
    
    infile = util.open_stream(filename)
    
    index = {}
    
    for line in util.SafeReadIter(infile):
        if line.startswith(">"):
            index[line[1:].rstrip()] = infile.tell()
    
    return index


def has_fasta_index(fasta_file):
    """Check to see if fasta_file has an index"""

    return os.path.exists(fasta_file + ".index")


def guessFastaWidth(fastaFile):
    fafile = util.open_stream(fastaFile, "rb")
    
    numlines = 5
    lineno = 0
    width = -1
    width2 = -1
    maxwidth = 0
    
    for line in fafile:
        if len(line) != 0  and line[0] != ">":
            lineno += 1        
            width3 = len(line.rstrip())
            maxwidth = max(maxwidth, width3)
            
            if width == -1:
                # first line
                width = width3

            elif width3 > width:
                # widths cannot get bigger
                return -1
            
            elif width3 == width:
                return width
            
            elif width2 == -1:
                # this should be last line in sequence
                width2 = width3
                return width
            else:
                # width got smaller twice
                return -1
        else:
            # previous sequence had only one line
            # rest widths for next sequence
            if width2 != -1:
                width2 = -1
            else:
                width = -1
    
    return maxwidth
    


class FastaIndex:
    def __init__(self, *filenames):
        self.filelookup = {}
        self.index = {}
        
        for fn in filenames:
            self.read(fn)
        
    
    def read(self, filename):
        # open fasta
        infile = util.open_stream(filename, "rb")
        
        # estimate column width
        self.width = guessFastaWidth(filename)
        if self.width == -1:
            raise Exception("lines do not have consistent width")
        
        # read index
        keys = []
        for key, start, end in util.DelimReader(filename + ".index", delim="\t"):
            keys.append(key)
            self.index[key] = (int(start), int(end))
            self.filelookup[key] = infile
        
        # return keys read
        return keys
    
    
    def get(self, key, start=1, end=None, strand=1):
        """Get a sequence by key
           coordinates are 1-based and end is inclusive"""
    
        assert start > 0, Exception("must specify coordinates one-based")
        assert key in self.index, Exception("key '%s' not in index" % key)
        
        if end != None and end < start:
            return ""
        
        # must translate from one-based to zero-based
        # must account for newlines
        filestart, fileend = self.index[key]
        start -= 1
        seek = filestart + start + (start // self.width)
        
        # if seek is past sequence then return empty sequence
        if seek >= fileend:
            return ""
        
        # seek to beginning
        infile = self.filelookup[key]
        infile.seek(seek)
        
        # read until end of sequence
        seq = []
        if end == None:
            lenNeeded = util.INF
        else:
            lenNeeded = end - start
        
        len2 = 0
        while len2 < lenNeeded:
            line = infile.readline()
            if line.startswith(">") or len(line) == 0:
                break
            seq.append(line.rstrip())
            len2 += len(seq[-1])
            if len2 > lenNeeded:
                seq[-1] = seq[-1][:-int(len2 - lenNeeded)]
                break
        seq = "".join(seq)
        
        # reverse complement if needed
        if strand == -1:
            seq = _revcomp(seq)
        
        return seq


    



#=============================================================================
# FASTA BLAST Indexing
#

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


def hasBlastIndex(fastaFile):
    """Check to see if fastaFile has a formatdb fasta index"""

    return os.path.exists(fastaFile + ".psd") and \
           os.path.exists(fastaFile + ".psi")



#=============================================================================
# Simple FASTA's as lists
# (not used very often)
#

def readFastaOrdered(filename, keyfunc=firstword, valuefunc=lambda x:x):
    """Read a FASTA file into a 'keys' and 'values' lists"""
    
    infile = util.open_stream(filename)
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


def writeFastaOrdered(filename, names, seqs, width=None):
    """Write a FASTA in array style to a file"""
    
    out = util.open_stream(filename, "w")
    
    for name, seq in izip(names, seqs):
        print >>out, ">%s" % name
        util.printwrap(seq, width, out=out)



def array2dict(names, seqs):
    """Convert array style FASTA to FastDict"""
    
    fa = FastaDict()
    for name, seq in zip(names, seqs):
        fa.add(name, seq)
    return fa


def dict2array(fa, order = None):
    """Convert FastaDict to array style FASTA"""
    
    names = []
    seqs = []
    
    if order == None:
        order = fa.keys()
    
    for name in order:
        names.append(name)
        seqs.append(fa[name])
    return (names, seqs)




#=============================================================================
# Special alignment format
# (rarely used. only used it for MIKE's alignments once)
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

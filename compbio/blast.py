# python imports
import os

# rasmus imports
from rasmus import util

# compbio imports
from . import fasta




# NCBI BLASTALL 2.2.10 -m 8 tab-delimited output
# Fields: 
# 0. Query id, 
# 1. Subject id, 
# 2. % identity, 
# 3. alignment length, 
# 4. mismatches, 
# 5. gap openings, 
# 6. q. start, 
# 7. q. end, 
# 8. s. start, 
# 9. s. end, 
# 10. e-value, 
# 11. bit score
#

# TODO: add checks to see if BLAST crashed


class BlastReader:
    """A parser for Blast results in -m8 format"""

    def __init__(self, filename, moreFunc = lambda: False):
        self.infile = util.open_stream(filename)
        self.moreFunc = moreFunc
    
    def __iter__(self):
        return self
    
    def read(self):
        line = self.infile.readline()
        
        # skip comments and blanks but not EOF
        while len(line) > 0 and line[0] in "#\n":
            line = self.infile.readline()
        
        return line.rstrip().split("\t")
    
    def next(self):
        line = self.read()
        while len(line) < 12:
            # if this stream runs dry, try to get another stream
            self.infile = self.moreFunc()
            if not self.infile:
                raise StopIteration
            line = self.read()
        return line


class BlastListReader (BlastReader):
    """A parser for Blast results sorted in multiple files"""

    def __init__(self, filenames):        
        self.infiles = filenames[1:]
        BlastReader.__init__(self, filenames[0], self.nextFile)
    
    def nextFile(self):
        if len(self.infiles) > 0:
            infile = self.infiles[0]
            self.infiles = self.infiles[1:]
            return util.open_stream(infile)
        else:
            return False



#
# funcions for parsing a Blast Hit
#


def query(hit):
    return hit[0]

def subject(hit):
    return hit[1]

def percentIdentity(hit):
    return float(hit[2])

def alignLength(hit):
    return int(hit[3])

def mismatches(hit):
    return int(hit[4])

def gapOpenings(hit):
    return int(hit[5])

def queryStart(hit):
    return int(hit[6])
    
def queryEnd(hit):
    return int(hit[7])

def queryLength(hit):
    return int(hit[7]) - int(hit[6])

def subjectStart(hit):
    return int(hit[8])

def subjectEnd(hit):
    return int(hit[9])

def subjectLength(hit):
    return int(hit[9]) - int(hit[8])

def evalue(hit):
    return float(hit[10])

def bitscore(hit):
    return float(hit[11])



def flipHit(hit):
    """Returns a new hit where query and subject are flipped"""
    return [hit[1],    # 0. Query id, 
            hit[0],    # 1. Subject id, 
            hit[2],    # 2. % identity, 
            hit[3],    # 3. alignment length, 
            hit[4],    # 4. mismatches, 
            hit[5],    # 5. gap openings, 
            hit[8],    # 6. q. start, 
            hit[9],    # 7. q. end, 
            hit[6],    # 8. s. start, 
            hit[7],    # 9. s. end, 
            hit[10],   # 10. e-value, 
            hit[11]]   # 11. bit score




def iterBestHitPerTarget(reader, scorefunc=bitscore):
    """Filters blast hits such that only the best hit between a query and a 
       target is kept.
       
       reader    -- BlastReader
       out       -- file stream
       scorefunc -- function of 1 argument that accepts a hit and returns a
                    score.  "best" score is considered max score.
    """
    lasthit = None
    last = []
    topscore = 0
    
    for line in reader:
        hit = (query(line), subject(line))
        if hit != lasthit and len(last) != 0:
            yield last
            last = line
            topscore = scorefunc(line)
        else:
             score = scorefunc(line)
             if score > topscore:
                topscore = score
                last = line
        lasthit = hit

def filterBestHitPerTarget(reader, out, scorefunc=bitscore):
    """Filters blast hits such that only the best hit between a query and a 
       target is kept.
       
       reader    -- BlastReader
       out       -- file stream
       scorefunc -- function of 1 argument that accepts a hit and returns a
                    score.  "best" score is considered max score.
    """
    lasthit = None
    last = []
    topscore = 0
    
    for line in reader:
        hit = (query(line), subject(line))
        if hit != lasthit and len(last) != 0:
            print >>out, "\t".join(last)
            last = line
            topscore = scorefunc(line)
        else:
             score = scorefunc(line)
             if score > topscore:
                topscore = score
                last = line
        lasthit = hit



def blastp(databaseFile, queryFile, options = "", split=100, resume = None):
    return blast("blastp", databaseFile, queryFile, options = options, 
                 split=split, resume=resume)

def blast(prog, databaseFile, queryFile, options = "", split=100, resume = None):
    """Executes blastp in several smaller batches"""

    if not split:
        # do blasting in one call
        pipe = os.popen("blastall -p %s -d %s -i %s -m 8 %s" % \
            (prog, databaseFile, queryFile, options))
        return BlastReader(pipe)
        
    else:
        # NOTE: split query file into about 100 sequences each
        # this is a work around for ncbi blastall 2.2.10 problem with outputing
        # in -m 8 mode.  error was   "BioseqFindFunc: couldn't uncache"
               
        seqs = fasta.read_fasta(queryFile)
        closure = {
            "index": 0,
            "oldtmp": None,
            "time": 0.0
            }
        
        if resume:
            try:
                closure["index"] = seqs.keys().index(resume)
                util.log("resuming with query '%s' (%d of %d)" % (
                    (resume, closure["index"], len(seqs.keys()))))
            except ValueError:
                raise Exception("Could not resume from last query sequence '%s'" % resume)
        
        
        def processFunc():
            # remove old query tempfile if one exists
            if closure["oldtmp"] != None:
                os.remove(closure["oldtmp"])
                elapse = util.toc()
                closure["time"] += elapse
                
                util.log("blasted %d of %d sequences (%.1f%%), elapse %.0f m, left %.0f m" % (
                    closure["index"], len(seqs.keys()), 
                    100 * float(closure["index"]) / len(seqs.keys()),
                    closure["time"] / 60.0, 
                    elapse / split * (len(seqs.keys()) - closure["index"]) / 60.0))
                
            util.tic()
            
            # find new subset of query sequences
            i = closure["index"]
            names = seqs.keys()[i:i+split]
            
            # if no more sequences then quit
            if len(names) == 0:
                return False
            
            # start blast
            tmpfile = util.tempfile(".", "blastp", ".fasta")
            seqs.write(tmpfile, names = names)
            pipe = os.popen("blastall -p %s -d %s -i %s -m 8 -e .1 %s" % \
                (prog, databaseFile, tmpfile, options))
            
            # update variables
            closure["oldtmp"] = tmpfile
            closure["index"] = i + split
            
            return pipe
         
        filename = processFunc()
        if filename:
            return BlastReader(filename, processFunc)
        else:
            return BlastReader(os.popen("less"))



def bl2seq(seq1, seq2, program="blastp", options="", name1="seq1", name2="seq2"):
    """
    Performs Blast between two sequences 'seq1' and 'seq2'.
    Returns a single Blast hit line or None if no hits found.
    """

    # create temp files for sequences
    file1 = util.tempfile(".", "blastp", ".fasta")
    file2 = util.tempfile(".", "blastp", ".fasta")

    fasta.write_fasta_ordered(file1, [name1], [seq1])
    fasta.write_fasta_ordered(file2, [name2], [seq2])

    
    # execute blast    
    pipe = os.popen("bl2seq -p %s -i %s -j %s -D 1 %s" % \
        (program, file1, file2, options))
    
    # parse hit
    hit = None
    for line in pipe:
        if line[0] == "#":
            continue
        hit = line.rstrip().split("\t")
        break
    
    # remove temp files
    os.remove(file1)
    os.remove(file2)
    
    return hit


def bestBidir(hits, scorefunc=bitscore):
    "find best bidirectional hits"
    
    best = util.Dict(default=[None, 0, None])
    
    for hit in hits:
        gene1 = query(hit)
        gene2 = subject(hit)
        score = scorefunc(hit)
        if score > best[gene1][1]:
            best[gene1] = [gene2, score, hit]
        if score > best[gene2][1]:
            best[gene2] = [gene1, score, hit]
    
    mark = set()
    hits2 = []
    for gene1, (gene2, score, hit) in best.iteritems():
        if best[gene2][0] == gene1 and gene1 not in mark:
            mark.add(gene1)
            mark.add(gene2)
            hits2.append(hit)
    
    return hits2




#=============================================================================
# These databases have not been used much 
# Fri Aug 19 13:24:50 EDT 2011


class BlastDb:
    """A database interface to Blast results"""
    
    def __init__(self, filename):
        import bsddb
        
        self.db = bsddb.db.DB()
        self.db.set_flags(bsddb.db.DB_DUP)
        
        self.db.open(filename, bsddb.db.DB_HASH, bsddb.db.DB_CREATE)
        self.cur = self.db.cursor()  
    
    def __getitem__(self, key):
        hits = {}
        val = self.cur.get(key, None, bsddb.db.DB_SET)
        
        while val:
            tokens = val[1].split("\t")
            hits[tokens[1]] = tokens
            val = self.cur.get(key, None, bsddb.db.DB_NEXT_DUP)
        return hits
    
    def __setitem__(self, key, hits):
        self.db.delete(key)
        for hit in hits.values():
            self.addHit(key, hit)
    
    def __delitem__(self, key):
        self.db.delete(key)
    
    def keys(self):
        return util.sort(util.unique(self.db.keys()))
    
    def has_key(self, key):
        return self.db.has_key(key)
    
    def close(self):
        return self.db.close()
    
    def __len__(self):
        return len(self.db.keys())
    
    
    def addHits(self, reader):
        for hit in reader:
            self.addHit(query(hit), hit)
            self.addHit(subject(hit), flipHit(hit))
    
    def addHit(self, key, hit):
        self.db.put(key, "\t".join(hit))



class BlastDb2:
    """A database interface to Blast results"""
    
    def __init__(self, filename):
        import sqlite
    
        self.con = sqlite.connect(filename)
        self.cur = self.con.cursor()
        
        # create hits table if it does not exist
        self.cur.execute("SELECT name FROM sqlite_master WHERE type = 'table'")
        tables = self.cur.fetchall()
        
        if ("hits",) not in tables:
            self.cur.execute("""
                CREATE TABLE hits (
                    query TEXT,
                    subject TEXT,
                    percid NUMBER, 
                    alignlen INT, 
                    mismatches INT,
                    gaps INT,
                    qstart INT,
                    qend INT,
                    sstart INT,
                    send INT,
                    evalue NUMBER,
                    bitscore NUMBER
                )
                """)
    
    def __del__(self):
        self.close()
    
    def __getitem__(self, key):
        self.cur.execute("SELECT * FROM hits WHERE query='%s' OR subject='%s'" % 
                         (key, key))
        
        hits = {}
        for hit in self.cur:
            if hit[0] == key:
                other = hit[1]
            else:
                other = hit[0]
            hits[other] = hit
    
        return hits
    
    """def __setitem__(self, key, hits):
        self.db.delete(key)
        for hit in hits.values():
            self.addHit(key, hit)
    
    def __delitem__(self, key):
        self.db.delete(key)
    
    def keys(self):
        return util.sort(util.uniq(self.db.keys()))
    
    def has_key(self, key):
        return self.db.has_key(key)
    """
    
    def close(self):
        self.con.commit()
        return self.con.close()
    
    """
    def __len__(self):
        return len(self.db.keys())
    
    """
    
    def addHits(self, reader):
        for hit in reader:
            self.addHit(hit)
    
    
    def addHit(self, hit):
        self.cur.execute("""
            INSERT INTO hits VALUES(
                '%s', '%s', %s, 
                %s, %s, %s, 
                %s, %s, %s, %s,
                %s, %s)
            """ % tuple(map(str, hit)))

    

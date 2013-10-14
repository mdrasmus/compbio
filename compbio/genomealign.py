

from rasmus import util
from rasmus import tablelib


from . import fasta, alignlib




class GenomeAlign (object):
    def __init__(self, master_file=None, seq2species=lambda x: x):
        self.lookup = util.Dict(default=[])
        self.seq2species = seq2species
        
        if master_file != None:
            self.read(master_file)
    
    
    def read(self, master_file):
        for row in tablelib.iter_table(master_file):
            self.lookup[(row['species'], row['chromosome'])].append(row)


    def get(self, species, chrom, start, end):
        records = []
    
        for record in self.lookup[(species, chrom)]:
            if util.overlap(start, end, record["start"], record["end"]):
                records.append(record)
            
        return records


    def get_files(self, species, chrom, start, end):
        return [x['filename'] for x in self.get(species, chrom, start, end)]
    
    
    def get_aligns(self, species, chrom, start, end, 
                  mainspecies=lambda keys: keys[0],
                  collapse=False):
        """By default assumes main species is 1st sequence"""
        
        # get records for this region
        records = self.get(species, chrom, start, end)
        records.sort(key=lambda x: x["start"])
        
        # read alignments
        alns = []
        for record in records:
            aln = fasta.read_fasta(record["filename"])
            
            # collapse alignment
            if collapse:
                ind = util.findneq("-", aln[mainspecies(aln.keys())])
                
                for key, seq in aln.iteritems():
                    if len(seq) != 0:
                        aln[key] = "".join(util.mget(seq, ind))
            
            l2a = alignlib.local2align(aln[mainspecies(aln.keys())])
            
            # trim front
            if start > record["start"]:
                trimstart = l2a[start - record["start"]]
            else:
                trimstart = 0
            
            # trim end
            if end < record["end"]:
                trimend = l2a[-(record["end"]-end)]
            else:
                trimend = aln.alignlen()
            
            # perform trim
            for key, seq in aln.iteritems():
                aln[key] = seq[trimstart:trimend]
                
            alns.append(aln)
        
        return alns
        
    
    def get_align(self, species, chrom, start, end,
                 mainspecies=lambda keys: keys[0],
                 collapse=False):
        # NOTE: assume one2one alignment
        
        cataln = fasta.FastaDict()
        alns = self.get_aligns(species, chrom, start, end,
                               mainspecies=mainspecies, collapse=collapse)
        
        if len(alns) == 0:
            return cataln
            
        lens = []
        
        for i, aln in enumerate(alns):
            lens.append(max(map(len, aln.values())))
            
            for key, seq in aln.iteritems():
                sp = self.seq2species(key)

                # start new species                
                if sp not in cataln:
                    cataln[sp] = []

                # catch up with blank sequences                    
                if len(cataln[sp]) < i:
                    for j in xrange(len(cataln[sp]), i):
                        cataln[sp].append("-" * lens[j])
                
                if len(seq) == 0:
                    cataln[sp].append("-" * lens[-1])
                else:
                    cataln[sp].append(seq)
        
        for sp in cataln:
            if len(cataln[sp]) < len(alns):
                # catch up with blank sequences
                for j in xrange(len(cataln[sp]), len(alns)):
                    cataln[sp].append("-" * lens[j])
        
        for key, seq in cataln.iteritems():
            cataln[key] = "".join(cataln[key])
        
        return cataln
    

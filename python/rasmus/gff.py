
"""
GTF file format
http://genes.cse.wustl.edu/GTF22.html

<seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]

The following feature types are required: "CDS", "start_codon", "stop_codon".
The features "5UTR", "3UTR", "inter", "inter_CNS", "intron_CNS" and "exon" are
optional. All other features will be ignored. The types must have the correct
capitalization shown here.

<start> <end> 
Integer start and end coordinates of the feature relative to the
beginning of the sequence named in <seqname>.  <start> must be less than or equal
to <end>. Sequence numbering starts at 1. Values of <start> and <end> that extend
outside the reference sequence are technically acceptable, but they are
discouraged.

<score> 
The score field indicates a degree of confidence in the feature's
existence and coordinates. The value of this field has no global scale but may
have relative significance when the <source> field indicates the prediction
program used to create this annotation. It may be a floating point number or
integer, and not necessary and may be replaced with a dot.

<frame> 
0 indicates that the feature begins with a whole codon at the 5' most
base. 1 means that there is one extra base (the third base of a codon) before the
first whole codon and 2 means that there are two extra bases (the second and
third bases of the codon) before the first codon. Note that for reverse strand
features, the 5' most base is the <end> coordinate.
"""

import util
import algorithms


TEST_GTF = \
"""
140\tTwinscan\tinter\t5141\t8522\t.\t-\t.\tgene_id ""; transcript_id "";
140\tTwinscan\tinter_CNS\t8523\t9711\t.\t-\t.\tgene_id ""; transcript_id "";
"""


class Region:
    """
    A generic GFF region
    """
    
    def __init__(self, region=None):
        if isinstance(region, Region):
            # copy information from other region
            self.seqname = region.seqname
            self.source  = region.source
            self.feature = region.feature
            self.start   = region.start
            self.end     = region.end
            self.score   = region.score
            self.strand  = region.strand
            self.frame   = region.frame
            self.attrs   = dict(region.attrs)
            self.comment = comment
            
        elif isinstance(region, str):
            self.read(region)
            
        else:
            raise Exception("Cannot parse argument of type '%s'" %
                            str(type(region)))
        
        # initializes species to None
        self.species = None
    
    
    def set(self,
            seqname="", 
            source="",
            feature="",
            start=0,
            end=0,
            score=None,
            strand=1,
            frame=None,
            attrs=None,
            comment=""):
        self.seqname = seqname
        self.source  = source
        self.feature = feature
        self.start   = start
        self.end     = end
        self.score   = score
        self.strand  = strand
        self.frame   = frame

        if attrs == None:
            self.attrs = {}
        else:
            self.attrs = attrs

        self.comment = comment
    
    
    def read(self, line):
        # parse comment
        if "#" in line:
            pos = line.index("#")
            self.comment = line[pos+1:]
            line = line[:pos]
        else:
            self.comment = ""
        
        # split into columns
        tokens = line.split("\t")
        assert len(tokens) == 9, Exception("line does not have 9 columns")
        
        # parse fields
        self.seqname = tokens[0]
        self.source  = tokens[1]
        self.feature = tokens[2]
        self.start   = int(tokens[3])
        self.end     = int(tokens[4])
        
        # parse score
        if tokens[5] == ".":
            self.score = None
        else:
            self.score = float(tokens[5])
        
        # parse strand
        if tokens[6] == "+":
            self.strand = 1
        elif tokens[6] == "-":
            self.strand = -1
        else:
            self.strand = None

        # parse frame
        if tokens[7] == ".":
            self.frame = None
        else:
            self.frame = int(tokens[7])

        # parse attributes
        self.attrs = self.parseAttrs(tokens[8])
    
    
    def __str__(self):
        if self.score == None:
            score = "."
        else:
            score = str(self.score)
        
        if self.strand == None:
            strand = "."
        elif self.strand == 1:
            strand = "+"
        else:
            strand = "-"
        
        if self.frame == None:
            frame = "."
        else:
            frame = str(self.frame)
        
        attr = self.formatAttrs(self.attrs)
        
        if self.comment != "":
            comment = " #%s" % self.comment
        else:
            comment = ""
        
        
        return "%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s%s" % \
            (self.seqname, 
             self.source, 
             self.feature,
             self.start,
             self.end,
             score,
             strand,
             frame,
             attr,
             comment)
    
    def formatAttrs(self, attrs):
        # unparsed attribute
        assert attrs.keys() == [None]
        return attrs[None]
    
    def parseAttrs(self, text):    
        return {None: text}
    
    
    def __repr__(self):
        return self.__str__()



class GtfRegion (Region):
    """
    GTF Region
    
    Parses and formats attributes field according to the GTF specification
    """

    def __init__(self, region=None):
        Region.__init__(self, region)
        

    def formatAttrs(self, attrs):
        lst = []
        for key, val in attrs.items():
            lst.append('%s "%s";' % (key, str(val)))
        return " ".join(lst)
    
    
    def parseAttrs(self, text):
        """Parses an attribute field into a dict of key/value pairs"""

        tokens = text.split(";")
        attrs = {}

        for attr in tokens[:-1]:
            attr = attr.strip()

            pos = attr.index(" ")

            key = attr[:pos]
            val = attr[pos+1:].split("\"")[1]

            attrs[key] = val

        return attrs


class Gff3Region (Region):
    """
    GFF3 Region
    
    Parses and formats attributes field according to the GFF3 specification
    """
    
    def __init__(self, region=None):
        Region.__init__(self, region)


    def formatAttrs(self, attrs):
        lst = []
        for key, val in attrs.items():
            lst.append('%s=%s;' % (key, str(val)))
        return "".join(lst)

    def parseAttrs(self, text):
        """Parses an attribute field into a dict of key/value pairs"""

        tokens = text.split(";")
        attrs = {}

        for attr in tokens:
            attr = attr.strip()
            if len(attr) == 0:
                continue

            pos = attr.index("=")

            key = attr[:pos]
            val = attr[pos+1:]

            attrs[key] = val

        return attrs



class RegionIter:
    """An iterator that walks down a sorted list of regions"""

    def __init__(self, regions, start, end, index=None):
        self.regions  = regions
        self.start  = start
        self.end    = end
        self.nregions = len(regions)
        
        if index != None:
            self.index = index
        else:
            self.index = 0
            
            # find starting index by binary search
            low, top = algorithms.binsearch(regions, start-1, 
                                            lambda a,b: cmp(a.start, b))
            
            if top != None:
                self.index = top
            else:
                self.index = self.nregions
    
    def __iter__(self):
        return self
    
    
    def next(self):
        if (self.index < self.nregions) and \
           (self.regions[self.index].start < self.end):
            gene = self.regions[self.index]
        else:
            raise StopIteration
        
        self.index += 1
        return gene



def overlap(region1, region2):
    return region1.seqname == region2.seqname and \
           region1.species == region2.species and \
           util.overlap(region1.start, region1.end,
                        region2.start, region2.end)

def overlaps(region1, regions):
    return [x for x in regions if overlap(region1, x)]



def readGff(filename, format=Region, 
            lineFilter=lambda x: True,
            regionFilter=lambda x: True):
    """
    Read all regions in a GFF file
    """
    
    infile = iterGff(filename,
                     format, 
                     lineFilter,
                     regionFilter)
    
    regions = []
    
    for region in infile:
        regions.append(region)
    
    return regions


def writeGff(filename, regions):
    out = util.openStream(filename, "w")
    
    for region in regions:
        print >>out, region


def iterGff(filename, format=Region, 
            lineFilter=lambda x: True,
            regionFilter=lambda x: True):
    """
    Iterate over the regions in a GFF file
    """
    
    infile = util.openStream(filename)
    
    for line in infile:
        line = line.rstrip()
        
        # only continue processing if line is not comment and passes filter
        if len(line) == 0 or line[0] == "#" or not lineFilter(line):
            continue
        
        # parse region        
        region = format(line)
        
        # only return region if region passes filter
        if regionFilter(region):
            yield region


if __name__ == "__main__":
    from rasmus.common import *
    
    pc(readGff(strStream(TEST_GTF)))


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



GFF3 File Format
http://song.sourceforge.net/gff3.shtml

"""

import sys

import util
import algorithms


TEST_GTF = \
"""
140\tTwinscan\tinter\t5141\t8522\t.\t-\t.\tgene_id ""; transcript_id "";
140\tTwinscan\tinter_CNS\t8523\t9711\t.\t-\t.\tgene_id ""; transcript_id "";
"""


TEST_GFF3 = \
"""
chr2\tTwinscan\tmRNA\t5141\t8522\t.\t-\t.\tID=gene1;
chr2\tTwinscan\texon\t8523\t9711\t.\t-\t.\tID=exon1; Parent=gene1;
chr2\tTwinscan\texon\t8523\t9711\t.\t-\t.\tID=exon2; Parent=gene1;
chr2\tTwinscan\texon\t8523\t9711\t.\t-\t.\tID=exon3; Parent=gene1;
"""

if __name__ == "__main__":
    import re
    TEST_GFF3_2 = \
re.sub(" +", "\t", """
##gff-version   3
##sequence-region   ctg123 1 1497228       
ctg123	. gene            1000  9000  .  +  .  ID=gene00001;Name=EDEN

ctg123	. TF_binding_site 1000  1012  .  +  .  ID=tfbs00001;Parent=gene00001

ctg123	. mRNA            1050  9000  .  +  .  ID=mRNA00001;Parent=gene00001;Name=EDEN.1
ctg123	. mRNA            1050  9000  .  +  .  ID=mRNA00002;Parent=gene00001;Name=EDEN.2
ctg123	. mRNA            1300  9000  .  +  .  ID=mRNA00003;Parent=gene00001;Name=EDEN.3

ctg123	. exon            1300  1500  .  +  .  ID=exon00001;Parent=mRNA00003
ctg123	. exon            1050  1500  .  +  .  ID=exon00002;Parent=mRNA00001,mRNA00002
ctg123	. exon            3000  3902  .  +  .  ID=exon00003;Parent=mRNA00001,mRNA00003
ctg123	. exon            5000  5500  .  +  .  ID=exon00004;Parent=mRNA00001,mRNA00002,mRNA00003
ctg123	. exon            7000  9000  .  +  .  ID=exon00005;Parent=mRNA00001,mRNA00002,mRNA00003

ctg123	. CDS             1201  1500  .  +  0  ID=cds000011;Parent=mRNA00001;Name=edenprotein.1
ctg123	. CDS             3000  3902  .  +  0  ID=cds000012;Parent=mRNA00001;Name=edenprotein.1
ctg123	. CDS             5000  5500  .  +  0  ID=cds000013;Parent=mRNA00001;Name=edenprotein.1
ctg123	. CDS             7000  7600  .  +  0  ID=cds000014;Parent=mRNA00001;Name=edenprotein.1

ctg123	. CDS             1201  1500  .  +  0  ID=cds000021;Parent=mRNA00002;Name=edenprotein.2
ctg123	. CDS             5000  5500  .  +  0  ID=cds000022;Parent=mRNA00002;Name=edenprotein.2
ctg123	. CDS      7000  7600  .  +  0  ID=cds000023;Parent=mRNA00002;Name=edenprotein.2

ctg123	. CDS             3301  3902  .  +  0  ID=cds000031;Parent=mRNA00003;Name=edenprotein.3
ctg123	. CDS    5000   5500   . +  2  ID=cds000032;Parent=mRNA00003;Name=edenprotein.3
ctg123	. CDS    7000  7600  .  +  2  ID=cds000033;Parent=mRNA00003;Name=edenprotein.3

ctg123	. CDS      3391  3902  .  +  0  ID=cds000041;Parent=mRNA00003;Name=edenprotein.4
ctg123	. CDS         5000  5500  .  +  2  ID=cds000042;Parent=mRNA00003;Name=edenprotein.4
Ctg123	. CDS      7000  7600    .  +  2  ID=cds000043;Parent=mRNA00003;Name=edenprotein.4
""")


#-------------------------------------------------------------------------------
# Region Classes
#-------------------------------------------------------------------------------

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
        
        elif region is None:
            # set defaults
            self.seqname = ""
            self.source  = ""
            self.feature = ""
            self.start   = 0
            self.end     = 0
            self.score   = None
            self.strand  = None
            self.frame   = None
            self.attrs   = {}
            self.comment = ""
        
        else:
            raise Exception("Cannot parse argument of type '%s'" %
                            str(type(region)))
        
        # initializes species to None
        self.species = None
    
    
    def read(self, line):
        # parse comment
        pos = line.find("#")
        if pos > -1:
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
    
    
    def write(self, out=sys.stdout):
        out.write(str(self))
        out.write("\n")
        
    
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


class RegionHierarchy (Region):
    """
    Generic Region class with support for parent and child regions
    """
    
    def __init__(self, region=None):
        Region.__init__(self, region)
        self.parents = []
        self.children = []
    
    
    def addChild(self, child):
        assert self not in child.parents
        assert child not in self.children
        child.parents.append(self)
        self.children.append(child)
    
    def removeChild(self, child):
        child.parents.remove(self)
        self.children.remove(child)
    
    def write(self, out=sys.stdout):
        out.write(str(self))
        out.write("\n")
        
        for child in self.children:
            child.write(out)
        
    
    @staticmethod
    def buildHierarchy(self, regions):
        """
        Produces a hierachy from a list of regions
        Returns list of regions with no parents (roots).
        
        This base class does not implement anything
        """
        
        return []



class GtfRegion (RegionHierarchy):
    """
    GTF Region
    
    Parses and formats attributes field according to the GTF specification
    """

    def __init__(self, region=None):
        RegionHierarchy.__init__(self, region)
        

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
    
    @staticmethod
    def buildHierarchy(regions):
        """
        Produces a hierachy from a list of regions
        Returns list of regions with no parents (roots).
        
        NOT written yet
        """
        
        # TODO: create hierarchy from gene_id, transcript_id attrs
        
        pass



class Gff3Region (RegionHierarchy):
    """
    GFF3 Region
    
    Parses and formats attributes field according to the GFF3 specification
    
    Also includes support for Parents and Children.
    """
    
    def __init__(self, region=None):
        RegionHierarchy.__init__(self, region)


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
    
    @staticmethod
    def buildHierarchy(regions):
        """
        Produces a hierachy from a list of regions
        Returns list of regions with no parents (roots).
        
        Assumes ID and Parent attributes are present.
        """
        
        # make a list of regions in case regions is not a list
        if not isinstance(regions, list):
            regions = list(regions)
        
        # make lookup by id
        roots = set()
        lookup = {}
        for region in regions:
            if "ID" in region.attrs:
                lookup[region.attrs["ID"]] = region
            roots.add(region)
        
        # build hierarchy
        for region in regions:
            if "Parent" in region.attrs:
                parents = region.attrs["Parent"].split(",")
                for parent in parents:
                    lookup[parent].addChild(region)
                roots.remove(region)
        
        # create roots list (regions in same order as they were passed)
        regions2 = [x for x in regions if x in roots]
        
        return regions2
        


#-------------------------------------------------------------------------------
# Region functions
#-------------------------------------------------------------------------------


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
    """Tests whether two regions overlap"""
    
    return region1.seqname == region2.seqname and \
           region1.species == region2.species and \
           util.overlap(region1.start, region1.end,
                        region2.start, region2.end)

def overlaps(region1, regions):
    """Find the regions in list 'regions' that overlap region1"""
    
    return [x for x in regions if overlap(region1, x)]


def regionLookup(regions, key="ID"):
    """
    Returns a dict lookup of regions based on a key (default: ID)
    """
    
    lookup = {}

    for region in regions:
        if key in region.attrs:
            rkey = region.attrs[key]
            assert rkey not in lookup, Exception("duplicate key")
            lookup[rkey] = region
            

    return lookup




#-------------------------------------------------------------------------------
# Region Input/Output
#-------------------------------------------------------------------------------

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
    """
    Write regions to a file stream
    
    filename - a filename or file stream
    regions  - a list of Region objects
    """
    
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
    lineno = 0
    
    for line in infile:
        lineno += 1
        line = line.rstrip()
        
        # only continue processing if line is not comment and passes filter
        if len(line) == 0 or line[0] == "#" or not lineFilter(line):
            continue
        
        # parse region
        try:
            region = format(line)
        except Exception, e:
            raise Exception("%s\nError on line %d: %s" % (e,lineno, line))
            
        
        # only return region if region passes filter
        if regionFilter(region):
            yield region



#
# testing
#
if __name__ == "__main__":
    from rasmus.common import *
    
    regions = readGff(strStream(TEST_GFF3_2), format=Gff3Region)
    
    regions2 = Gff3Region.buildHierarchy(regions)
    
    print regions2
    regions2[0].write()
    
    if 0:
        pc(readGff(strStream(TEST_GTF)))




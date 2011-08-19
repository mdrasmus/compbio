
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

    this format says a feature can have multiple parents

"""


# python imports
import sys

# rasmus imports
from rasmus import util

# compbio imports
from compbio import regionlib


#=============================================================================
# Generic GFF fileformat
#



class Gff (object):
    """Most generic GFF format.  Do not assume any format for attributes field"""
    
    def __init__(self):
        self.nondata = set(["comment", "source", "score", "frame", "species"])
    
    
    def format_data(self, region, ignore=set()):
        # unparsed attribute
        return region.data.get(None, "")

    def parse_data(self, text, ignore=set()):
        return {None: text}


    def read_region(self, line, region=None):
        if region == None:
            region = regionlib.Region()

        # parse comment
        pos = line.find("#")
        if pos > -1:
            region.data["comment"] = line[pos+1:]
            line = line[:pos]


        # split into columns
        tokens = line.split("\t")
        assert len(tokens) == 9, Exception("line does not have 9 columns")

        # parse fields
        region.seqname = tokens[0]
        region.feature = tokens[2]
        region.start   = int(tokens[3])
        region.end     = int(tokens[4])

        # parse strand
        strand = tokens[6]
        if strand == "+" or strand == "1":
            region.strand = 1
        elif strand == "-" or strand == "-1":
            region.strand = -1
        else:
            region.strand = 0

        # parse source
        if tokens[1] != ".":
            region.data["source"]  = tokens[1]

        # parse score
        if tokens[5] != ".":
            region.data["score"] = float(tokens[5])

        # parse frame
        if tokens[7] != ".":
            region.data["frame"] = int(tokens[7])

        # parse attributes
        region.data.update(self.parse_data(tokens[8]))
        
        # parse species
        region.species = region.data.get("species", "")

        return region


    def write_region(self, region, out=sys.stdout):
        score = str(region.data.get("score", "."))
        source = str(region.data.get("source", "."))

        if region.strand == 0:
            strand = "."
        elif region.strand == 1:
            strand = "+"
        else:
            strand = "-"

        frame = str(region.data.get("frame", "."))

        attr = self.format_data(region)

        if "comment" in region.data:
            comment = " #%s" % region.data["comment"]
        else:
            comment = ""


        out.write("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s%s\n" % \
                  (region.seqname, 
                   source, 
                   region.feature,
                   region.start,
                   region.end,
                   score,
                   strand,
                   frame,
                   attr,
                   comment))
    
    def build_hierarchy(self, regions):
        """
        Produces a hierachy from a list of regions
        Returns list of regions with no parents (roots).
        
        This base class function does nothing.  See GFF3
        """
        
        # do nothing
        return []

# singleton
GFF = Gff()

#=============================================================================
# GTF fileformat
#

class Gtf (Gff):

    def format_data(self, region):
        lst = []
        
        if region.species != "":
            lst.append('species "%s";' % region.species)
        
        for key, val in region.data.iteritems():
            if key not in self.nondata:
                lst.append('%s "%s";' % (key, str(val)))
        return " ".join(lst)


    def parse_data(self, text):
        """Parses an attribute field into a dict of key/value pairs"""

        tokens = text.split(";")
        data = {}

        for attr in tokens[:-1]:
            attr = attr.strip()

            pos = attr.find(" ")
            if pos == -1:
                continue

            key = attr[:pos]
            val = attr[pos+1:].split("\"")[1]

            data[key] = val

        return data
    
    def build_hierarchy(self, regions):
        """GTF has its own heirarchy system
           It is currently not implemented"""
        
        return []

GTF = Gtf()

#=============================================================================
# GFF3 fileformat
#

class Gff3 (Gff):

    def format_data(self, region):
        lst = []
        
        if region.species != "":
            lst.append("species=%s;" % region.species)
        
        for key, val in region.data.iteritems():
            if key not in self.nondata:
                lst.append('%s=%s;' % (key, str(val)))
        return "".join(lst)


    def parse_data(self, text):
        """Parses an attribute field into a dict of key/value pairs"""

        tokens = text.split(";")
        data = {}

        for attr in tokens:
            attr = attr.strip()
            if len(attr) == 0:
                continue

            pos = attr.index("=")

            key = attr[:pos]
            val = attr[pos+1:]

            data[key] = val

        return data
    
    
    def build_hierarchy(self, regions):
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
            if "ID" in region.data:
                lookup[region.data["ID"]] = region
            roots.add(region)
        
        # build hierarchy
        for region in regions:
            if "Parent" in region.data:
                parents = region.data["Parent"].split(",")
                for parent in parents:
                    lookup[parent].add_child(region)
                roots.remove(region)
        
        # create roots list (regions in same order as they were passed)
        regions2 = [x for x in regions if x in roots]

        return regions2

GFF3 = Gff3()


#=============================================================================
# Gff Input/Output
#

def read_gff(filename, format=GFF3, 
            lineFilter=lambda x: True,
            regionFilter=lambda x: True):
    """
    Read all regions in a GFF file
    """
    
    infile = iterGff(filename,
                     format, 
                     lineFilter,
                     regionFilter)
    
    return list(infile)
readGff = read_gff


def write_gff(filename, regions, format=GFF3):
    """
    Write regions to a file stream
    
    filename - a filename or file stream
    regions  - a list of Region objects
    """
    
    out = util.open_stream(filename, "w")
    
    for region in regions:
        format.write_region(region, out=out)
writeGff = write_gff


def iter_gff(filename, format=GFF3, 
             line_filter=lambda x: True,
             region_filter=lambda x: True,
             # backcompat
             lineFilter=None,
             regionFilter=None):
    """
    Iterate over the regions in a GFF file
    """

    if lineFilter is not None:
        line_filter = lineFilter
    if regionFilter is not None:
        region_filter = regionFilter

    
    infile = util.open_stream(filename)
    lineno = 0
    
    for line in infile:
        lineno += 1
        line = line.rstrip("\n")
        
        # only continue processing if line is not comment and passes filter
        if len(line) == 0 or line[0] == "#" or not line_filter(line):
            continue
        
        # parse region
        try:
            region = format.read_region(line)
        except Exception, e:
            raise Exception("%s\nError on line %d: %s" % (e,lineno, line))
            
        
        # only return region if region passes filter
        if region_filter(region):
            yield region
iterGff = iter_gff


#
# testing
#
if __name__ == "__main__":
    from rasmus.common import *
    import re

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

    regions = read_gff(strStream(TEST_GFF3_2), format=GFF3)
    regions2 = GFF3.build_hierarchy(regions)
    
    print regions2
    print regions2[0]
    
    pc(read_gff(strStream(TEST_GTF)))

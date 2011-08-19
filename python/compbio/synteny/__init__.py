
import copy

from rasmus import util
from compbio import regionlib



class SyntenyBlock (object):
    def __init__(self, region1, region2, name="", data={}):

        self.name = ""
        self.region1 = region1 # total region span by block in species1
        self.region2 = region2 # total region span by block in species2

        # direction is parallel (1) or anti-parallel (-1)
        self.dir = region1.strand * region2.strand

        # ordered list of ortholog pairs
        self.orths = []

        # extra data
        self.data = copy.copy(data)
    

    def get_direction(self):
        """Returns 1 if region1.strand == region2.strand and
                  -1 if region1.strand != region2.strand
        """
        return self.dir
    
    
    def add_orth(self, orth, direction=None):
        """Initialize synteny block with an ortholog pair"""
        self.orths.append(orth)

        if direction is not None:
            self.dir = direction



    def recalc_regions(self, regiondb):
        """Recalculate the regions of a synteny block from its
           orthologs pairs"""

        self.region1.start = regiondb.get_region(self.orths[0][0][0]).start
        self.region1.end = regiondb.get_region(self.orths[-1][0][-1]).end

        if self.dir == 1 or self.dir == 0:
            self.region2.start = regiondb.get_region(self.orths[0][1][0]).start
            self.region2.end = regiondb.get_region(self.orths[-1][1][-1]).end
        else:
            self.region2.start = regiondb.get_region(self.orths[-1][1][0]).start
            self.region2.end = regiondb.get_region(self.orths[0][1][-1]).end



def make_orth(db, genes1, genes2):
    """
    Returns gene names as a tuple of two sorted lists (by start pos)

    The ortholog format is ((gene1a, gene1b, gene1c), (gene2a, gene2b))
    """
    
    genes1 = list(genes1)
    genes2 = list(genes2)
    genes1.sort(key=lambda x: db.get_region(x).start)
    genes2.sort(key=lambda x: db.get_region(x).start)
    return (genes1, genes2)


def score_block_bbh_num(block):
    """Score a block by the number of BBH it contains"""

    return len(block_bbh_hits(block))


def score_block_bbh_sum(block):
    """Score a block by its sum BBH"""

    return sum(hit[2] for hit in block_bbh_hits(block))


def block_bbh_hits(block):
    """Score a block by the number of BBH it contains"""

    # find all unidirectional best hits
    best = util.Dict(default=[-util.INF, None, None])

    for hit in block.data["hits"]:
        a, b, val = hit[:3]
        a = a.data["ID"]
        b = b.data["ID"]

        if val > best[a][0]:
            best[a] = (val, b, hit)
        if val > best[b][0]:
            best[b] = (val, a, hit)

    # count bi-directional best hits
    hits2 = []
    for a, (val, b, hit) in best.iteritems():
        if best[b][1] == a and a < b:
            hits2.append(hit)

    return hits2



#=============================================================================
# input/output

def write_synteny_blocks(out, blocks, extra=lambda x: ()):
    """Write a list of synteny blocks to file"""

    # write blocks
    for block in blocks:
        out.write("\t".join(map(str, (
                            block.region1.species,
                            block.region1.seqname,
                            block.region1.start,
                            block.region1.end,
                            block.region2.species,
                            block.region2.seqname,
                            block.region2.start,
                            block.region2.end,
                            block.dir,
                            block.name) +
                            extra(block))) + "\n")



def read_synteny_blocks(filename, feature="synteny",
                        extra=lambda r, cols:None):
    infile = util.open_stream(filename)
    blocks = []
    
    for line in infile:
        tokens = line.split("\t")
        species1, chrom1, start1, end1, \
        species2, chrom2, start2, end2, direction = tokens[:9]
        if len(tokens) > 9:
            name = tokens[9]
        else:
            name = ""
        
        blocks.append(SyntenyBlock(
            regionlib.Region(species1, chrom1, feature, 
                             int(start1), int(end1), 1),
            regionlib.Region(species2, chrom2, feature, 
                             int(start2), int(end2), int(direction)),
            name=name))

        extra(blocks[-1], tokens[10:])
        
    return blocks





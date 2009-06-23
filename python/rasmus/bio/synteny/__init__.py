
from rasmus import regionlib, util




class SyntenyBlock (object):
    def __init__(self, region1, region2):
        
        self.region1 = region1 # total region span by block in species1
        self.region2 = region2 # total region span by block in species2

        # direction is parallel (1) or anti-parallel (-1)
        self.dir = region1.strand * region2.strand

        # ordered list of ortholog pair
        self.orths = []
    

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




#=============================================================================
# input/output

def write_synteny_blocks(out, blocks):
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
                            block.dir))) + "\n")



def read_synteny_blocks(filename, feature="synteny"):
    infile = util.open_stream(filename)
    blocks = []
    
    for line in infile:
        species1, chrom1, start1, end1, \
        species2, chrom2, start2, end2, direction = line.split("\t")
        
        blocks.append(SyntenyBlock(
            regionlib.Region(species1, chrom1, feature, 
                             int(start1), int(end1), 1),
            regionlib.Region(species2, chrom2, feature, 
                             int(start2), int(end2), int(direction))))
    return blocks





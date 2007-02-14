# python libs
import copy

# rasmus lib
from rasmus import util
from rasmus import algorithms



class Region (object):
    """
    A generic sequence region
    """
    
    def __init__(self, species="",
                       seqname="",
                       feature="",
                       start=0,
                       end=0,
                       strand=0,
                       data={}):
        if isinstance(species, Region):
            region = species
            
            # copy information from other region
            self.species = region.species
            self.seqname = region.seqname
            self.feature = region.feature
            self.start   = region.start
            self.end     = region.end
            self.strand  = region.strand
            self.data   = copy.copy(region.data)
            self.parents = copy.copy(region.parents)
            self.children = copy.copy(region.children)
        
        else:
            self.species = species
            self.seqname = seqname
            self.feature = feature
            self.start   = start
            self.end     = end
            self.strand  = strand
            self.data   = copy.copy(data)
            self.parents = []
            self.children = []

    def __iter__(self):
        return iter(self.children)
    
    
    def __getitem__(self, key):
        return self.children[key]
    
    
    def __setitem__(self, key,  val):
        self.children[key] = val
    
    def __delitem__(self, key):
        del self.children[key]
    
        
    def length(self):
        return self.end - self.start + 1

    def addChild(self, child):
        assert self not in child.parents
        assert child not in self.children
        child.parents.append(self)
        self.children.append(child)
    
    def removeChild(self, child):
        child.parents.remove(self)
        self.children.remove(child)

    
    def __repr__(self):
        return self.__str__()


    def __str__(self):
        if self.strand == 0:
            strand = "."
        elif self.strand == 1:
            strand = "+"
        else:
            strand = "-"
        
        data = str(self.data)
        
        
        text = ["[%s\t%s\t%s\t%d\t%d\t%s\t%s]" % \
            (self.species,
             self.seqname, 
             self.feature,
             self.start,
             self.end,
             strand,
             data)]
        
        for child in self.children:
            text.append(str(child))
        
        return "\n".join(text)



def findRegionPos(regions, pos):
    """Find the first region that starts after 'pos' in a sorted list of 'regions'"""
    low, top = algorithms.binsearch(regions, pos-1, lambda a,b: cmp(a.start, b))
    return top
    


def findRegion(regions, region):
    """Find a region in a sorted list of 'regions'"""
    low, ind = algorithms.binsearch(regions, region.start-1, 
                                    lambda a,b: cmp(a.start, b))
    if ind == None:
        return None
    
    while ind < len(regions) and regions[ind] != region:
        ind += 1
    
    if ind == len(regions):
        return None
    else:
        return ind
    

# TODO: add reverse iteration
def iterChrom(regions, start, end, index=None):
    """An iterator that walks down a sorted list of regions"""

    nregions = len(regions)
    
    # find index
    if index == None:
        # find starting index by binary search
        index = findRegionPos(regions, start)
        if index == None:
            return
    
    # walk down chromosome
    while  (index < nregions) and \
           (regions[index].start < end):
        yield regions[index]
        index += 1



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
        if key in region.data:
            rkey = region.data[key]
            assert rkey not in lookup, Exception("duplicate key")
            lookup[rkey] = region

    return lookup



class RegionDb (object):
    """Organize regions for easy access"""
    
    # TODO: this is starting to look like my Matching class in genomeutil :)
    
    def __init__(self, regions):
        
        self.species = {}
        self.regions = {}
    
        for region in regions:
            sp = self.species.setdefault(region.species, {})
            chrom = sp.setdefault(region.seqname, [])
            chrom.append(region)
            
            if "ID" in region.data:
                self.regions[region.data["ID"]] = region
            
    
        for sp, chroms in self.species.iteritems():
            for chrom, regs in chroms.iteritems():
                regs.sort(key=lambda x: x.start)


class EndPoint:
    def __init__(self, region, boundary):
        self.region = region
        self.boundary = boundary



# python libs
import copy

# rasmus lib
from rasmus import util



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

    def add_child(self, child):
        assert self not in child.parents
        assert child not in self.children
        child.parents.append(self)
        self.children.append(child)

    
    def remove_child(self, child):
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



def find_region_pos(regions, pos):
    """Find the first region that starts after 'pos' in a sorted list of 'regions'"""
    low, top = util.binsearch(regions, pos-1, lambda a,b: cmp(a.start, b))
    return top
    


def find_region(regions, region):
    """Find a region in a sorted list of 'regions'"""
    low, ind = util.binsearch(regions, region.start-1, 
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
def iter_chrom(regions, start, end, index=None):
    """An iterator that walks down a sorted list of regions"""

    nregions = len(regions)
    
    # find index
    if index == None:
        # find starting index by binary search
        index = find_region_pos(regions, start)
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


def groupby_overlaps(regions, bygroup=True):
    """regions must be sorted by chrom and start"""

    species = None
    seqname = None
    start = -util.INF
    end = -util.INF
    group = None
    groupnum = -1
    for reg in regions:
        
        if (reg.species != species or
            reg.seqname != seqname or
            reg.start > end):
            
            # start new group
            species = reg.species
            seqname = reg.seqname
            start = reg.start
            end = reg.end
            groupnum += 1

            if bygroup:
                if group is not None:
                    yield group
                group = [reg]
            else:
                yield (groupnum, reg)

        else:
            # append to current group
            if reg.end > end:
                end = reg.end

            if bygroup:
                group.append(reg)
            else:
                yield (groupnum, reg)

    
        
        


def region_lookup(regions, key="ID"):
    """
    Returns a dict lookup of regions based on a key (default: ID)
    """
    
    lookup = {}

    for region in regions:
        rkey = None
        
        if isinstance(key, basestring):
            if key in region.data:
                rkey = region.data[key]
        else:
            rkey = key(region)

        if rkey is not None:
            assert rkey not in lookup, Exception("duplicate key '%s'" % rkey)
            lookup[rkey] = region
            

    return lookup



class RegionDb (object):
    """Organize regions for easy access"""
    
    # TODO: this is starting to look like my Matching class in genomeutil :)
    
    def __init__(self, regions):
        
        self.sp2chroms = {} # {species -> {chrom -> regions sorted by start}}
        self.regions = {}   # {region_id -> region}
        self.positions = {} # {region_id -> (species, chrom, position)}
        

        # sort regions into chromosomes
        for region in regions:
            sp = self.sp2chroms.setdefault(region.species, {})
            chrom = sp.setdefault(region.seqname, [])
            chrom.append(region)

            # record region id
            if "ID" in region.data:
                self.regions[region.data["ID"]] = region
            
        # sort each chromosome by start position
        for sp, chroms in self.sp2chroms.iteritems():
            for chrom, regs in chroms.iteritems():
                regs.sort(key=lambda x: x.start)

        # make index lookups
        for sp, chroms in self.sp2chroms.iteritems():
            self.positions[sp] = {}
            for chrom, regs in chroms.iteritems():
                for i, reg in enumerate(regs):
                    self.positions[reg.data["ID"]] = (sp, chrom, i)


    def has_species(self, species):
        return species in self.sp2chroms

    def has_chrom(self, species, chrom):
        return chrom in self.sp2chroms.get(species, {})

    def get_species(self):
        """Returns all species in database"""
        return self.sp2chroms.keys()
    
    def get_chroms(self, species):
        """Get list of chromosomes for a given species"""
        return self.sp2chroms[species]

    def get_regions(self, species, chrom):
        if species not in self.sp2chroms:
            return []        
        return self.sp2chroms[species].get(chrom, [])

    def get_all_regions(self):
        return self.regions.values()

    def iter_all_regions(self):
        return self.regions.itervalues()

    def get_region(self, regionid):
        """Returns the region with regionid""" 
        return self.regions[regionid]

    def has_region(self, regionid):
        return regionid in self.regions

    def get_region_pos(self, regionid):
        """Returns position of region along chromosome"""
        return self.positions[regionid][2]

    def get_region_pos_full(self, regionid):
        """Returns (species, chrom, position) of a region"""
        reg = self.regions[regionid]
        return (reg.species, reg.seqname, self.positions[regionid][2])


class EndPoint:
    def __init__(self, region, boundary):
        self.region = region
        self.boundary = boundary



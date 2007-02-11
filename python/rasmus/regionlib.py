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

    
    


def overlap(region1, region2):
    """Tests whether two regions overlap"""
    
    return region1.seqname == region2.seqname and \
           region1.species == region2.species and \
           util.overlap(region1.start, region1.end,
                        region2.start, region2.end)


class EndPoint:
    def __init__(self, region, boundary):
        self.region = region
        self.boundary = boundary



"""

FFF (Frequent Feature Format)

<feature name> <species> <chrom> <len> <starts ...>

FFI (Frequent Feature Index)

<feature name> <species> <chrom> <len> <start> <character pos>

"""

from rasmus import util



class FeatureRef:
    def __init__(self, name, start, filename, pos):
        self.name     = name
        self.start    = start
        self.filename = filename
        self.pos      = pos


class Feature:
    def __init__(self, name, length):
        self.name = name
        self.length = length        


class FeatureIndex:
    def __init__(self):
        self.data = util.Dict(dim=2, default=[])
        self.infiles = {}
        self.features = {}
    
    
    def read(self, * filenames):
        for filename in filenames:
            indexfile = filename.replace(".fff", ".ffi")
            self.infiles[filename] = file(filename)

            for line in file(indexfile):
                feature, species, chrom, length, start, pos = line.rstrip().split("\t")

                self.data[(species, chrom)][feature].append(
                    FeatureRef(feature, 
                               int(start), filename, 
                               int(pos)))

                self.features[feature] = Feature(feature, int(length))
        
        
        # ensure features are sorted
        for features in self.data.itervalues():
            for sites in features.itervalues():
                sites.sort(key=lambda x: x.start)
    
            

    def lookup(self, species, chrom, start, end):
        key = (species, chrom)
        if key not in self.data:
            return []
        else:
            features = self.data[key]
            lookups = []
            
            for feature, sites in features.iteritems():
                first, junk = util.binsearch(sites, start, 
                                             lambda x, coord:
                                             cmp(x.start, coord))
                
                if first == None:
                    first = 0
                
                i = first
                while i < len(sites) and sites[i].start < end:
                    lookups.append(sites[i])
                    i += 1
            
            return lookups
    
    
    def getFeatures(self, species, chrom, start, end):
        lookups = self.lookup(species, chrom, start, end)
        features = util.Dict(1, [])
        
        for l in lookups:
            infile = self.infiles[l.filename]
            infile.seek(l.pos)
            
            
            starts = []
            while True:
                token, char = util.readUntil(infile, " \t\n")
                if char == "\n":
                    break
                
                coord = int(token)
                
                if abs(coord) > end:
                    break
                
                if abs(coord) > start:
                    starts.append(coord)
            
            features[l.name].extend(starts)
        
        return features
                    
            


def readFffIndex(* filenames):
    findex = FeatureIndex()
    findex.read(* filenames)
    return findex


def makeTestFff(species, chrom, length, out):
    import random
    
    for i in range(50):
        out.write("motif%d\t%s\t%s\t%d\t" % (i, species, chrom, length))
        
        pos = 0
        for j in range(100000):
            pos += random.randint(1, 1000)
            out.write("%d " % pos)
        out.write("\n")


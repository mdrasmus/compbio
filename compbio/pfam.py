
from rasmus import util

class Domain (object):
    def __init__(self, *args):
    
        if len(args) == 1:            
            tokens = args[0].rstrip().split()

            self.name   = tokens[0]
            #self.domain = tokens[1]
            self.start  = int(tokens[2])
            self.end    = int(tokens[3])
            self.score  = float(tokens[8])
            self.evalue = float(tokens[9])
        else:
            self.name, self.start, self.end, self.evalue = args

    def __repr__(self):
        return "Domain(%s,\t%d,\t%d,\t%f)\n" % (self.name, self.start, self.end, self.evalue)
    
    def __str__(self):
        return self.__repr__()


def iterPfam(filename):
    infile = util.open_stream(filename)
    
    
    def getQuery(infile):
        for line in infile:
            if line.startswith("Query sequence"):
                name = line.rstrip().replace("Query sequence: ", "")
                return name
    
    def getDomains(infile):
        domains = []
        
        for line in infile:
            if line.startswith("Parsed for domains:"):
                break
                
        infile.next()   # skip header 1
        infile.next()   # skip header 2

        for line in infile:
            if len(line) <= 1 or line[0] in "\t ":
                break
            domains.append(Domain(line))
        
        return domains
    
    while True:
        query = getQuery(infile)
        if query == None:
            break
        
        domains = getDomains(infile)
        
        yield query, domains



def readPfam(filename):
    result = {}
    
    for query, domains in iterPfam(filename):
        result[query] = domains
    
    return result

from xml.sax import make_parser
from xml.sax.handler import feature_namespaces
import xml.sax.handler


def readGo(filename):
    """DEPRECATED"""
    terms = Dict(default=[])
    
    for line in file(filename):
        if "GI:" in line:# or "KEGG:" in line:
            continue
        tokens = line.rstrip().split("\t")
        try:
            terms[tokens[0]].append(tokens[4])
        except:
            print line
    
    return terms


def readCommonNames(filename):
    """DEPRECATED"""
    commonNames = {}

    for line in file(filename):
        tokens = line.rstrip().split("\t")

        if tokens[1] != '-':
            commonNames[tokens[0]] = tokens[1]
    return commonNames



class GoTerm:
    def __init__(self):
        self.accession = ""
        self.name = ""
        self.definition = ""
        self.is_a = []
        self.part_of = []
        self.synonym = []

class AllTerm(GoTerm):
    def __init__(self):
        GoTerm.__init__(self)
        
        self.accession = "all"
        self.name = "all"
        self.defintion = "top-level term" 

class GoHandler(xml.sax.handler.ContentHandler):
    def __init__(self, base):
        self.terms = {}
        self.term = None
        self.elm = ""
        self.base = base
    
    def startElement(self, name, attrs):
        if name == "go:term":
            self.term = GoTerm()
        elif name == "go:is_a":
            ref = attrs["rdf:resource"]
            if ref.startswith(self.base):
                self.term.is_a.append(ref[len(self.base):])
        elif name == "go:part_of":
            ref = attrs["rdf:resource"]
            if ref.startswith(self.base):
                self.term.part_of.append(ref[len(self.base):])
        self.elm = name
    
    def endElement(self, name):
        if name == "go:term":
            self.terms[self.term.accession] = self.term
        self.elm = ""
    
    def characters(self, text):
        if self.elm == "go:accession":
            self.term.accession += text
        elif self.elm == "go:name":
            self.term.name += text
        elif self.elm == "go:definition":
            self.term.definition += text
        elif self.elm == "go:synonym":
            if 'GO' in text:
                self.term.synonym.append(text)
        

class GoDatabase:
    def __init__(self, filename):
        # Create a parser
        parser = make_parser()

        # Tell the parser we are not interested in XML namespaces
        parser.setFeature(feature_namespaces, 0)

        # Create the handler
        dh = GoHandler("http://www.geneontology.org/go#")

        # Tell the parser to use our handler
        parser.setContentHandler(dh)

        # Parse the input
        parser.parse(filename)

        self.terms = dh.terms
        
        # add top level term
        self.terms["all"] = AllTerm()
    
    
    def getAllParents(self, goid, touched=None, count=0, ret=True):
        if touched == None:
            touched = {}
        
        if goid in self.terms:
            term = self.terms[goid]
            parents =  term.is_a + term.part_of
            
            for parent in parents:
                if parent not in touched and parent != "all":
                    touched[parent] = count
                    count += 1
            
            for parent in parents:
                self.getAllParents(parent, touched, count, False)
        
        if ret:
            parents = touched.keys()
            parents.sort(key=lambda x: touched[x])
            return parents

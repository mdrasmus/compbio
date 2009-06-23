#!/usr/bin/env python

from rasmus import util
import sys

# xml support
from xml.sax import make_parser
import xml.sax.handler




class BibHandler(xml.sax.handler.ContentHandler):
    def __init__(self):
        self.papers = []
        self.attr = None
        self.inpaper = False
    
    def startElement(self, name, attrs):
        if name == "paper":
            self.inpaper = True
            self.papers.append({})
        elif self.inpaper:
            self.attr = name
    
    def endElement(self, name):
        if name == "paper":
            inpaper = False
    
    def characters(self, text):
        if self.attr:
            self.papers[-1][self.attr] = text
            self.attr = None



def readBibXml(filename):
    infile = util.open_stream(filename)

    # Create a parser
    parser = make_parser()

    # Tell the parser we are not interested in XML namespaces
    #parser.setFeature(feature_namespaces, 0)

    # Tell the parser to use our handler
    handler = BibHandler()
    parser.setContentHandler(handler)

    # Parse the input
    parser.parse(infile)
    
    return handler.papers

def writeBib(papers):
    mapping = {
        "author": "author",
        "title": "title",
        "date": "year",
        "pub": "journal",
        "pages": "pages",
        "text": "text"
    }
        

    for paper in papers:
        paper.setdefault("type", "article")
        if "key" not in paper:
            continue
    
        print "@%(type)s{%(key)s," % paper
        
        for key in paper:
            if key in mapping:
                print "  %s = \"%s\"," % (mapping[key], paper[key])
        
        print "}"
        print
        



for f in sys.argv[1:]:
    papers = readBibXml(f)
    writeBib(papers)
    

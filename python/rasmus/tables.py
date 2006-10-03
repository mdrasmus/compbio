"""
tables.py

Implements Manolis style tab-delimited table file format.



#types:string int
#default:'none' 0
#expect:
#delim:whitespace,space,tab,',',';'
#author:
#header:1
#
#
#
#
#
#
name num
data1 data2




"""

import sys
import StringIO

from util import *


class Table (list):
    def __init__(self, rows=None, headers=None):
        if headers == None:
            self.headers = []
        else:
            self.headers = headers
        
        self.types = []
        self.filename = None
        
        if rows != None:
            self.extend(rows)
            
            if headers == None:
                self.headers = sorted(self[0].keys())
    
    
    def add(self, **kargs):
        self.append(kargs)
    
    
    def read(self, filename, delim="\t"):
        """Reads a character delimited file and returns a list of dictionaries
            
           notes:
           Lines that start with '#' are treated as comments and are skiped
           Blank lines are skipped.

           If the first comment starts with '#Types:' the following tokens
           are interpreted as the data type of the column and values in that
           column are automatically converted.

           supported datatypes:
           - string
           - int
           - float
           - bool
           - unknown (no conversion is done, left as a string)

        """

        infile = openStream(filename)
        
        # remember filename for later saving
        if isinstance(filename, str):
            self.filename = filename
        
        self.headers = None
        self[:] = []
        
        
        for line in infile:
            line = line.rstrip()

            # skip blank lines
            if len(line) == 0:
                continue

            # handle comments
            if line[0] == "#":
                if line.startswith("#Types:"):
                    types = parseTableTypes(line, delim)
                continue

            tokens = line.split(delim)


            if not self.headers:
                # parse header
                self.headers = tokens
            else:
                # parse data
                row = {}
                for i in xrange(len(tokens)):
                    if len(types) > 0:
                        row[self.headers[i]] = types[i](tokens[i])
                    else:
                        row[self.headers[i]] = tokens[i]
                self.append(row)
        
        self.types = types
    
    
    def write(self, filename, delim="\t"):
        out = openStream(filename, "w")

        # set default header if needed
        if not self.headers:
            self.headers = self[0].keys()

        # write types    
        entry = self[0]
        types = map(lambda x: type(entry[x]), self.headers)
        out.write("#Types:" + formatTableTypes(types, delim) + "\n")

        # write header
        print >>out, delim.join(self.headers)

        # write data
        for row in self:
            print >>out, delim.join(map(lambda x: str(row[x]), self.headers))
    
    
    def save(self):
        if self.filename != None:
            self.write(self.filename)
        else:
            raise Exception("Table has no filename")
    

    def lookup(self, *keys):
        """return a lookup dict based on a column 'key'
           or multiple keys"""
        
        lookup = Dict(dim=len(keys))
        
        for row in self:
            keys2 = mget(row, keys)
            ptr = lookup
            for i in xrange(len(keys2) - 1):
                ptr = lookup[keys2[i]]
            if keys2[-1] in ptr:
                raise Exception("duplicate key '%s'" % str(keys2[-1]))
            ptr[keys2[-1]] = row
        
        lookup.insert = False
        return lookup
    
    
    def getMatrix(self):
        """returns a copy of the table as a 2D list"""
        
        mat = [self.headers]
        
        for row in self:
            mat.append(mget(row, self.headers))
        
        return mat
    
    
    def filter(self, cond):
        tab = Table(headers = self.headers)
        
        for row in self:
            if cond(row):
                tab.append(row)
        
        return tab
    
    
    def sort(self, cmp=cmp, key=None, reverse=False, field=None):
        if field != None:
            key = lambda row: row[field]        
        list.sort(self, cmp=cmp, key=key, reverse=reverse)
        

    def __repr__(self):
        s = StringIO.StringIO("w")
        printcols(self.getMatrix(), spacing=2, out=s)
        return s.getvalue()
    
    
    def __str__(self):
        s = StringIO.StringIO("w")
        self.write(s)
        return s.getvalue()
    
    


def manoli_str2bool(text):
    """Will parse every way manolis stores a boolean as a string"""
    
    text2 = text.lower()
    
    if text2 == "false":
        return False
    elif text2 == "true":
        return True
    else:
        raise Exception("unknown string for bool '%s'" % text)


tableTypesLookup = {
        "string": str,
        "int": int,
        "float": float,
        "bool": manoli_str2bool,
        "unknown": str
    }



def parseTableTypes(line, delim):
    names = line.replace("#Types:", "").split(delim)
    types = []
    
    for name in names:
        if name in tableTypesLookup:
            types.append(tableTypesLookup[name])
        else:
            types.append(eval(name))
    return types


def formatTableTypes(types, delim):
    lookup = revdict(tableTypesLookup)
    names = []
    
    for t in types:
        if t in lookup:
            names.append(lookup[t])
        else:
            names.append(t.__name__)
    return delim.join(names)
    

def readTable(filename, delim="\t"):
    table = Table()
    table.read(filename, delim=delim)
    return table

def writeTable(filename, table, delim="\t"):
    table.write(filename, delim)


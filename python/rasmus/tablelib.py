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
import copy

import util


# table directives
DIR_TYPES    = 0
DIR_HEADERS  = 1
DIR_DEFAULTS = 2
DIR_DELIM    = 3



class Table (list):
    def __init__(self, rows=[], 
                       headers=None,
                       defaults={},
                       types={},
                       filename=None):
        
        # set table info
        self.headers = copy.copy(headers)
        self.defaults = copy.copy(defaults)
        self.types = copy.copy(types)
        self.filename = filename
        self.comments = []
        
        
        # set data
        if len(rows) > 0:
            # list of dicts
            if isinstance(rows[0], dict):
                self.extend(rows)
                
                if self.headers == None:
                    self.headers = sorted(self[0].keys())
        
            # list of lists
            elif isinstance(rows[0], list):
                if self.headers == None:
                    raise Exception("must specify headers with 2D list")
                for row in rows:
                    self.append(dict(zip(self.headers, row)))
            
            # set info
            for key in self.headers:
                # set types
                if key not in self.types:
                    self.types[key] = type(self[0][key])
            
                # set defaults
                if key not in self.defaults:
                    self.defaults[key] = self.types[key]()
                    
    
    
    def new(self):
        """return a new table with the same info"""
        
        tab = Table(headers=self.headers,
                    types=self.types,
                    defaults=self.defaults,
                    filename=self.filename)
        
        tab.comments = self.comments
        
        return tab        
    
    
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

        infile = util.openStream(filename)
        
        # remember filename for later saving
        if isinstance(filename, str):
            self.filename = filename
        
        # clear table
        self.headers = None
        self.types = {}
        self.defaults = {}
        self[:] = []
        types = None
        
        
        for line in infile:
            line = line.rstrip()

            # skip blank lines
            if len(line) == 0:
                continue

            # handle comments
            if line[0] == "#":
                if line.startswith("#Types:") or \
                   line.startswith("#types:"):
                    types = parseTableTypes(line, delim)
                    self.comments.append(DIR_TYPES)
                else:
                    self.comments.append(line)
                continue
            
            # split row            
            tokens = line.split(delim)

            
            if not self.headers:
                # parse header
                self.headers = tokens
                if types:
                    assert len(types) == len(self.headers)
                    self.types = dict(zip(self.headers, types))
                else:
                    self.types = {}.fromkeys(self.headers, str)
                self.defaults = util.mapdict(self.types,
                                             valfunc=lambda x: x())
            else:
                # parse data
                row = {}
                for i in xrange(len(tokens)):
                    key = self.headers[i]
                    
                    if len(tokens[i]) == 0:
                        # default value
                        row[key] = self.defaults[key]
                    else:
                        row[key] = self.types[key](tokens[i])
                    
                # default values for incomplete rows
                for i in xrange(len(tokens), len(self.headers)):
                    row[key] = self.defaults[self.headers[i]]
                
                self.append(row)
    
    
    def write(self, filename, delim="\t"):
        out = util.openStream(filename, "w")
        
        # ensure all info is complete
        for key in self.headers:
            if key not in self.types:
                if len(self) > 0:
                    self.types[key] = type(self[0][key])
                else:
                    self.types[key] = str
            
            if key not in self.defaults:
                self.defaults[key] = self.types[key]()
            
        
        # ensure types are in directives
        if DIR_TYPES not in self.comments:
            self.comments = [DIR_TYPES] + self.comments

        # write comments
        for line in self.comments:
            if isinstance(line, str):
                print >>out, line
            elif line == DIR_TYPES:
                entry = self[0]
                out.write("#types:" +
                          formatTableTypes(util.mget(self.types,
                                                     self.headers),
                                           delim) + "\n")
        

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
        
        lookup = util.Dict(dim=len(keys))
        
        for row in self:
            keys2 = util.mget(row, keys)
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
            mat.append(util.mget(row, self.headers))
        
        return mat
    
    
    def filter(self, cond):
        tab = self.new()
        
        for row in self:
            if cond(row):
                tab.append(row)
        
        return tab
    
    
    def get(self, rows=None, cols=None):
        """get a table with a subset of the rows and columns"""
        
        tab = self.new()
        
        # determine rows and cols
        if rows == None:
            rows = range(len(self))
        
        if cols == None:
            cols = self.headers
        
        # copy over info
        tab.headers = copy.copy(cols)
        tab.types = util.subdict(self.types, cols)
        tab.defaults = util.subdict(self.defaults, cols)
        
        # copy data        
        for i in rows:
            row = {}
            for j in cols:
                row[j] = self[i][j]
            tab.append(row)
        
        return tab
    
    
    def sort(self, cmp=cmp, key=None, reverse=False, field=None):
        if field != None:
            key = lambda row: row[field]        
        list.sort(self, cmp=cmp, key=key, reverse=reverse)


    def __getslice__(self, a, b):
        tab = self.new()
        tab[:] = list.__getslice__(self, a, b)
        return tab
        

    def __repr__(self):
        s = StringIO.StringIO("w")
        util.printcols(self.getMatrix(), spacing=2, out=s)
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
        "str": str,
        "int": int,
        "float": float,
        "bool": manoli_str2bool,
        "unknown": str
    }



def parseTableTypes(line, delim):
    names = line.replace("#types:", "").replace("#Types:", "").split(delim)
    types = []
    
    for name in names:
        if name in tableTypesLookup:
            types.append(tableTypesLookup[name])
        else:
            types.append(eval(name))
    return types


def formatTableTypes(types, delim):
    lookup = util.revdict(tableTypesLookup)
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




if __name__ == "__main__":
    import StringIO
    
    text="""\
#
#types:str	int
name	num
matt	123
alex	456
mike	789
"""

    tab = readTable(StringIO.StringIO(text))

    

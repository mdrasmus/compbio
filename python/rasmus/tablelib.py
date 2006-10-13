"""
tables.py

Implements Manolis style tab-delimited table file format.


##version:1.0
##types:string int
##default:'none' 0
##header:1
#
#
#
##expect:
##author:
##delim:whitespace,space,tab,',',';'
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
DIR_VERSION  = 0
DIR_TYPES    = 1
DIR_HEADERS  = 2
DIR_DEFAULTS = 3
#DIR_DELIM    = 3

# a special null type 
NULL = util.Bundle()


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
        self.delim = "\t"
        self.nheaders = 1
        self.version = "1.0"
        
        
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
                    self.headers = range(len(rows[0]))
                    self.nheaders = 0
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
        tab.delim = self.delim
        tab.nheaders = self.nheaders
        
        return tab
    
    
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
        self[:] = []
        self.headers = None
        self.types = {}
        self.defaults = {}        
        self.comments = []
        self.delim = delim
        self.nheaders = 1
        self.version = "1.0"
        headerskip = 0
        
        
        # temps for reading only
        self.tmptypes = None
        self.tmpdefaults = None
        
        
        for line in infile:
            line = line.rstrip()

            # skip blank lines
            if len(line) == 0:
                continue

            # handle comments
            if line[0] == "#":
                if not self.readDirective(line):
                    self.comments.append(line)
                continue
            
            # split row            
            tokens = line.split(delim)

            
            if not self.headers:
                # parse header
                if self.nheaders > 0:
                    if headerskip >= self.nheaders - 1:
                        self.headers = tokens
                    else:
                        self.comments.append(line)
                        headerskip += 1
                        continue
                else:
                    self.headers = range(len(tokens))
                
                # populate types
                if self.tmptypes:
                    assert len(self.tmptypes) == len(self.headers)
                    self.types = dict(zip(self.headers, self.tmptypes))
                else:
                    self.types = {}.fromkeys(self.headers, str)
                
                # populate defaults
                if self.tmpdefaults:
                    self.defaults = {}
                    for header, default in zip(self.headers, self.tmpdefaults):
                        self.defaults[header] = self.types[header](default)
                else:        
                    self.defaults = util.mapdict(self.types,
                                                 valfunc=lambda x: x())
                
                # if we used this line as a header then go to next line
                if self.nheaders > 0:
                    continue
            
            
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
                key = self.headers[i]
                row[key] = self.defaults[key]

            self.append(row)
    
    
    def write(self, filename, delim="\t"):
        self.filename = filename
    
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
        
        # ensure version is in directives
        if DIR_VERSION not in self.comments:
            self.comments = [DIR_VERSION] + self.comments


        # write comments
        for line in self.comments:
            if isinstance(line, str):
                print >>out, line
            else:
                self.writeDirective(line, out, delim)
        

        # write header
        if self.nheaders > 0:
            print >>out, delim.join(self.headers)

        # write data
        for row in self:
            rowstr = []
            for header in self.headers:
                if header in row:
                    rowstr.append(str(row[header]))
                else:
                    rowstr.append('')
            print >>out, delim.join(rowstr)
    
    
    def save(self):
        if self.filename != None:
            self.write(self.filename)
        else:
            raise Exception("Table has no filename")
    
    
        
    
    def determineDirective(self, line):
        if line.startswith("##version:"):
            return DIR_VERSION
        elif line.startswith("#Types:") or \
             line.startswith("#types:") or \
             line.startswith("##types:"):
            # backwards compatible
            return DIR_TYPES
            
        elif line.startswith("##defaults:"):
            return DIR_DEFAULTS
            
        elif line.startswith("##headers:"):
            return DIR_HEADERS
            
        else:
            return None
    
    
    def readDirective(self, line):
        """Attempt to read a line with a directive"""
        
        directive = self.determineDirective(line)
        
        if directive == None:
            return False
        
        rest = line[line.index(":")+1:]         
        self.comments.append(directive)
        
        if directive == DIR_VERSION:
            self.version = rest
            return True
            
        elif directive == DIR_TYPES:
            self.tmptypes = parseTableTypes(rest, self.delim)
            return True
            
        elif directive == DIR_DEFAULTS:
            self.tmpdefaults = rest.split(self.delim)
            return True
            
        elif directive == DIR_HEADERS:
            self.nheaders = int(rest)
            #if self.nheaders not in [0, 1]:
            #    raise "Only 0 or 1 headers are allowed"
            return True
        
        else:
            return False
    
    
    def writeDirective(self, line, out, delim):
        """Write a directive"""
        
        if line == DIR_VERSION:
            out.write("##version:%s\n" % self.version)
        
        elif line == DIR_TYPES:
            entry = self[0]
            out.write("##types:" +
                      formatTableTypes(util.mget(self.types,
                                                 self.headers),
                                       delim) + "\n")
        elif line == DIR_DEFAULTS:
            out.write("##defaults:" +
                      delim.join(map(str, 
                                util.mget(self.defaults, self.headers))) + "\n")
        
        elif line == DIR_HEADERS:
            out.write("##headers:%d\n" % self.nheaders)
        
        else:
            raise "unknown directive:", line
    

    
    def add(self, **kargs):
        """Add a row to the table
           
           tab.add(col1=val1, col2=val2, col3=val3)
        """
        self.append(kargs)
    
    
    def addCol(self, header, coltype=str, default=NULL, pos=None):
        """Add a column to the table
        
           header  - name of the column
           coltype - type of the values in that column
           default - default value of the column
           pos     - position to insert column (default: right-end)
        """
        
        if header in self.headers:
            raise Exception("header '%s' is already in table" % header)
        
        if pos == None:
            pos = len(self.headers)
        if default == NULL:
            default = coltype()
        
        self.headers.insert(pos, header)
        self.types[header] = coltype
        self.defaults[header] = default
    
    
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
    
    def filter(self, cond):
        """Returns a table with a subset of rows such that cond(row) == True"""
        tab = self.new()
        
        for row in self:
            if cond(row):
                tab.append(row)
        
        return tab
    
    
    def get(self, rows=None, cols=None):
        """Returns a table with a subset of the rows and columns"""
        
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
    
    
    def cget(self, *cols):
        """Returns columns of the table as separate lists"""
        
        ret = []
        
        for col in cols:
            newcol = []
            ret.append(newcol)
            
            for row in self:
                newcol.append(row[col])
        
        return ret
    
    
    
    def getMatrix(self):
        """Returns a copy of the table as a 2D list"""
        
        mat = [self.headers]
        
        for row in self:
            mat.append(util.mget(row, self.headers))
        
        return mat
    
    
    def sort(self, cmp=cmp, key=None, reverse=False, col=None):
        if col != None:
            key = lambda row: row[col]        
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
    
    





#
# Types handling
#

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
    names = line.split(delim)
    types = []
    
    for name in names:
        if name in tableTypesLookup:
            types.append(tableTypesLookup[name])
        else:
            raise "unknown type '%s'" % name
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



#
# convenience functions
#

def readTable(filename, delim="\t"):
    table = Table()
    table.read(filename, delim=delim)
    return table



if __name__ == "__main__":
    import StringIO
    
    text="""\
##types:str	int	int
name	num	num2
matt	123	3
alex	456	
mike	789	1
"""

    tab = readTable(StringIO.StringIO(text))
    tab.addCol("extra", bool, True)
    
    print tab.defaults
    print tab
    print tab.cget('name', 'num')
    

    text="""\
##types:str	int	int
##defaults:none	0	0
##headers:0
#
# hello
#
name	0	1
matt	123	3
alex	456	
mike	789	1
"""

    tab = readTable(StringIO.StringIO(text))    
    
    print tab.defaults
    print tab
    print tab[0][1]


    text="""\
##types:str	int	int
##defaults:none	0	0
##headers:2
#
# hello
#
skip
name	0	1
matt	123	3
alex	456	
mike	789	1
"""

    tab = readTable(StringIO.StringIO(text))    
    
    print tab

    

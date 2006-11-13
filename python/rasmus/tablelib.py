"""
tables.py

Portable Tabular Format (PTF)

Implements Manolis style tab-delimited table file format.


--Example----------------------------------------------------
##version:1.0
##types:string int
##default:'none' 0
##header:1
#
#
# unimplemented for now
##expect:
##author:
##delim:whitespace,space,tab,',',';'
#
#
#
name num
data1 data2
-------------------------------------------------------------

File is tab delimited.

Directives are on a single line and begin with two hashes '##'
No space after colon is allowed.


"""

import sys
import StringIO
import copy

# rasmus libs
import util


# table directives
DIR_VERSION  = 0
DIR_TYPES    = 1
DIR_HEADERS  = 2
DIR_DEFAULTS = 3
#DIR_DELIM    = 3

# a special unique null type 
NULL = util.Bundle()


class TableException (Exception):
    """Exception class for Table"""
    pass


      


class Table (list):
    """Class implementing the Portable Table Format"""

    def __init__(self, rows=[], 
                       headers=None,
                       defaults={},
                       types={},
                       extraHeaders = [],
                       filename=None):
        
        # set table info
        self.headers = copy.copy(headers)
        self.defaults = copy.copy(defaults)
        self.types = copy.copy(types)
        self.extraHeaders = []
        self.filename = filename
        self.comments = []
        self.delim = "\t"
        self.nheaders = 1
        self.version = "1.0"
        
        
        # set data
        if len(rows) > 0:
            # data is a list of dicts
            if isinstance(rows[0], dict):
                for row in rows:
                    self.append(copy.copy(row))
                
                if self.headers == None:
                    self.headers = sorted(self[0].keys())
        
            # data is a list of lists
            elif isinstance(rows[0], list) or \
                 isinstance(rows[0], tuple):
                if self.headers == None:
                    self.headers = range(len(rows[0]))
                    self.nheaders = 0
                for row in rows:
                    self.append(dict(zip(self.headers, row)))
            
            
            # set table info
            for key in self.headers:
                # set types
                if key not in self.types:
                    self.types[key] = type(self[0][key])
            
                # set defaults
                if key not in self.defaults:
                    self.defaults[key] = self.types[key]()
        
        # set extra headers
        if len(extraHeaders) > 0:
            for row in extraHeaders:
                self.extraHeaders.append(copy.copy(row))
                
                    
    
    
    def new(self):
        """return a new table with the same info"""
        
        tab = Table(headers=self.headers,
                    types=self.types,
                    defaults=self.defaults,
                    filename=self.filename,
                    extraHeaders=self.extraHeaders)
        
        tab.comments = copy.copy(self.comments)
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
        extraHeaders = []
        
        
        # temps for reading only
        self.tmptypes = None
        self.tmpdefaults = None
        
        lineno = 0
        
        try:
            for line in infile:
                line = line.rstrip()        
                lineno += 1


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

                # if no headers read yet, use this line as a header
                if not self.headers:
                    # parse header
                    if self.nheaders > 0:
                        # determine if extra headers exist
                        if len(extraHeaders) >= self.nheaders - 1:
                            self.headers = tokens
                            
                            # check that headers are unique
                            check = set()
                            for header in self.headers:
                                if header in check:
                                    raise TableException("Duplicate header '%s'" % header)
                                check.add(header)
                        else:
                            # this line is an extra header
                            extraHeaders.append(tokens)
                            continue
                    else:
                        # default headers are numbers
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
                
        except Exception, e:
            e = TableException("line %d: %s" % (lineno, str(e)))
            raise e
        
        
        # now that we know the headers we can process extra headers
        for i, row in enumerate(extraHeaders):
            if len(row) != len(self.headers):
                raise TableException("wrong number of columns in extra header %d" % i)
            self.extraHeaders.append(dict(zip(self.headers, row)))
        
        
        # clear temps
        del self.tmptypes
        del self.tmpdefaults
    
    
    
    def write(self, filename, delim="\t"):
        """Write a table to a file or stream.
           
           If 'filename' is a string it will be opened as a file.
           If 'filename' is a stream it will be written to directly.
        """
        
        # remember filename for later saving
        if isinstance(filename, str):
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
        
        # write extra headers
        for row in self.extraHeaders:
            print >>out, delim.join(util.mget(row, self.headers))
        
        
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
        """Writes the table to the last used filename for the read() or write()
           function"""
        
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
            if len(self) > 0:
                entry = self[0]
            else:
                entry = [""] * len(self.headers)
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
        """Add a column to the table.  You must populate column data yourself.
        
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
        
        for row in self.extraHeaders:
            row[header] = ''
    
    
    def removeCol(self, *cols):
        """Removes a column from the table"""
        
        for col in cols:
            self.headers.remove(col)
            del self.types[col]
            del self.defaults[col]

            for row in self.extraHeaders:
                del row[col]
            
            for row in self:
                del row[col]
    
    
    def renameCol(self, oldname, newname):
        """Renames a column"""
        
        # change header
        col = self.headers.index(oldname)
        self.headers[col] = newname
        
        # change info
        self.types[newname] = self.types[oldname]
        del self.types[oldname]
        self.defaults[newname] = self.defaults[oldname]
        del self.defaults[oldname]       
        
        # change data
        for row in self:
            row[newname] = row[oldname]
            del row[oldname]
        
        # change extraHeaders
        for row in self.extraHeaders:
            row[newname] = row[oldname]
            del row[oldname]
            
    
    
    def lookup(self, *keys):
        """Returns a lookup dict based on a column 'key'
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
    
    
    def groupby(self, keyfunc):
        groups = {}
        
        for row in self:
            key = keyfunc(row)
            
            if key not in groups:
                groups[key] = Table(headers=self.headers)
            
            groups[key].append(row)
        
        return groups
    
    
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
        
        if len(ret) == 1:
            return ret[0]
        else:    
            return ret
    
    
    
    def getMatrix(self, rlabel=None):
        """Returns mat, rlabels, clabels
        
           where mat is a copy of the table as a 2D list
                 rlabels are the row labels
                 clabels are the column labels
        """
        
        # get labels
        if rlabel != None:
            rlabels = tab.cget(rlabel)
            clabels = copy.copy(self.headers)
            clabels.remove(rlabel)
        else:
            rlabels = range(len(self))
            clabels = copy.copy(self.headers)

        # get data
        mat = []
        for row in self:
            mat.append(util.mget(row, clabels))
        
        return mat, rlabels, clabels
    
    
    def sort(self, cmp=None, key=None, reverse=False, col=None):
        """Sorts the table inplace"""
        
        if col != None:
            key = lambda row: row[col]
        elif cmp == None and key == None:
            # sort by first column
            key = lambda row: row[self.headers[0]]
        
        list.sort(self, cmp=cmp, key=key, reverse=reverse)


    def __getitem__(self, key):
        if isinstance(key, slice):
            # return another table if key is a slice
            tab = self.new()
            tab[:] = list.__getitem__(self, key)
            return tab
        else:
            return list.__getitem__(self, key)
    
    
    def __getslice__(self, a, b):
        # for python version compatibility
        return self.__getitem__(slice(a, b))
    

    def __repr__(self):
        s = StringIO.StringIO("w")
        mat, rlabels, clabels = self.getMatrix()
        mat = [clabels] + mat
        util.printcols(mat, spacing=2, out=s)
        return s.getvalue()
    
    
    def __str__(self):
        s = StringIO.StringIO("w")
        self.write(s)
        return s.getvalue()
    
    





#
# Types handling
#

def manoli_str2bool(text=None):
    """Will parse every way manolis stores a boolean as a string"""
    
    if text == None:
        return False
    
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


tableTypesNames = {
    str: "string",
    int: "int",
    float: "float",
    bool: "bool"
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
    lookup = tableTypesNames
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



def histTable(items, headers=["item", "count", "percent"]):
    h = util.histDict(items)
    tab = Table(headers=headers)
    tot = float(len(items))

    if len(headers) == 2:    
        for key, val in h.items():
            tab.append({headers[0]: key,
                        headers[1]: val})
    
    elif len(headers) == 3:
        for key, val in h.items():
            tab.append({headers[0]: key,
                        headers[1]: val,
                        headers[2]: val / tot})
    
    else:
        raise Exception("Wrong number of headers (2 or 3 only)")
    
    tab.sort(col=headers[1], reverse=True)
    
    return tab


def matrix2table(matrix, rlabels=None, clabels=None, rowheader="labels"):
    if clabels == None:
        clabels = range(len(matrix[0]))
        nheaders = 0
    else:
        nheaders = 1
    
    if rlabels == None:
        tab = Table(headers=clabels)
    else:
        tab = Table(headers=[rowheader] + clabels)
    tab.nheaders = nheaders
    
   
    for i, row in enumerate(mat):
        if rlabels != None:
            row2 = {rowheader: rlabels[i]}
        else:
            row2 = {}
            
        for j in xrange(len(mat[i])):
            row2[clabels[j]] = mat[i][j]
        
        tab.append(row2)
    
    return tab



#
# testing
#

if __name__ == "__main__":
    import StringIO
    

    
    #################################################
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
    
    
    tab.addCol('extra', bool, False)
    for row in tab:
        row['extra'] = True
    
    
    #################################################
    text="""\
##types:str	int	int
##defaults:none	0	0
##headers:3
#
# hello
#
skip1	skip2	skip3
skip4	skip5	skip6
name	0	1
matt	123	3
alex	456	
mike	789	1
"""


    tab = readTable(StringIO.StringIO(text))    
    tab.renameCol('1', 'num')
    tab.removeCol('0')
    
    print tab
    
    print tab

    #################################################
    text="""\
##types:str	int	int
name	num	num2
matt	123	3
alex	456	
mike	789	1
"""

    tab = readTable(StringIO.StringIO(text))
    tab.sort()
    
    print repr(tab)
    print tab.defaults
    print tab
    print tab.cget('name', 'num')


    #################################################
    # catch parse error
    if 0:
        text="""\
##types:str	int	int
name	num	num
matt	123	0
alex	456	
mike	789	1
"""

        tab = readTable(StringIO.StringIO(text))
        tab.sort()

        print repr(tab)
        print tab.defaults
        print tab
        print tab.cget('name', 'num')

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


Table can also handle custom types.  Custom types must do the following

 1. default value: 
   default = mytype()  
   returns default value

 2. convert from string
   val = mytype(string)  
   converts from string to custom type

 3. convert to string
   string = str(val)
   converts val of type 'mytype' to a string
   TODO: I could change this interface...
   I could just use  mytype.__str__(val)

 4. type inference (optional)
   type(val)
   returns instance of 'mytype'
   TODO: I could not require this (only map() really needs it and __init__())


"""

import sys
import StringIO
import copy

# rasmus libs
from rasmus import util


TABLE_VERSION = "1.0"

# table directives
DIR_VERSION  = 0
DIR_TYPES    = 1
DIR_HEADERS  = 2
DIR_DEFAULTS = 3
#DIR_DELIM    = 3

# a special unique null type (more 'null' than None)
NULL = util.Bundle()


class TableException (Exception):
    """Exception class for Table"""
    def __init__(self, errmsg, filename=None, lineno=None):
        msg = ""
        add_space = False
        add_semicolon = False
        
        if filename:
            msg += "%s" % filename
            add_space = True
            add_semicolon = True
            
        if lineno:
            add_semicolon = True        
            if add_space:
                msg += " "
            msg += "line %d" % lineno
        
        if add_semicolon:
            msg += ": "
        
        msg = msg + errmsg
        
        Exception.__init__(self, msg)


#===========================================================================
# Types handling
#

def guessType(text):
    """Guesses the type of a value encoded in a string"""
    
    if text.isdigit():
        return int

    try:
        float(text)
        return float
    except ValueError:
        pass
    
    try:
        manoli_str2bool(text)
        return bool
    except ValueError:
        pass
    
    return str


def manoli_str2bool(text=None):
    """Will parse every way manolis stores a boolean as a string"""
    
    if text == None:
        # default value
        return False
    
    text2 = text.lower()
    
    if text2 == "false":
        return False
    elif text2 == "true":
        return True
    else:
        raise ValueError("unknown string for bool '%s'" % text)
        
            

class TableTypeLookup:
    def __init__(self, names_types=[]):
        self.name2type = {}
        self.type2name = {}
        
        for name, t in names_types:
            self.name2type[name] = t
            self.type2name[t] = name
    
    
    def extend(self, names_types):
        typeLookup = TableTypeLookup()
        typeLookup.name2type = copy.copy(self.name2type)
        typeLookup.type2name = copy.copy(self.type2name)
        
        for name, t in names_types:
            typeLookup.name2type[name] = t
            typeLookup.type2name[t] = name
        
        return typeLookup
        
    
    def parseTableTypes(self, line, delim):
        lookup = self.name2type
        names = line.split(delim)
        types = []

        try:
            for name in names:
                types.append(lookup[name])
        except KeyError:
            raise TableException("unknown type '%s'" % name)
        return types


    def formatTableTypes(self, types, delim):
        lookup = self.type2name
        names = []

        try:
            for t in types:
                names.append(lookup[t])
        except KeyError:
            raise TableException("unknown type '%s'" % t.__name__)
            #names.append(t.__name__)
        return delim.join(names)


class TableType:
    def __init__(self, parser=None, formatter=None):
        if parser != None:
            self.parserFunc = parser
        else:
            self.parserFunc = self.parser
        
        if formatter != None:
            self.formatterFunc = formatter
        else:
            self.formatterFunc = self.formatter
        
        
    def __call__(self, text=NULL):
        if text == NULL:
            return self.parserFunc()
        else:
            return self.parserFunc(text)
    
    def __str__(self, val):
        return self.formatterFunc(val)
    
    
    def parser(self, text=NULL):
        if text == NULL:
            return ""
        else:
            return text
    
    def formatter(self, val):
        return str(val)
        

def stringQuote(text=None):
    if text == None:
        return ""
    
    text2 = []
    
    i = 0
    while i < len(text):
        if text[i] == "\n":
            text2.append("\\n")
        elif text[i] == "\t":
            text2.append("\\t")
        elif text[i] == "\\":
            text2.append("\\\\")
        else:
            text2.append(text[i])
        
        i += 1
    
    return "".join(text2)
    

def stringExpand(text=None):
    if text == None:
        return ""
    
    print 'eee', text
    text2 = []
    
    i = 0
    while i < len(text):
        if text[i] == "\\" and i < len(text) - 1:
            if text[i+1] == "n":
                text2.append("\n")
            elif text[i+1] == "t":
                text2.append("\t")
            elif text[i+1] == "\\":
                text2.append("\\")
            else:
                text2.append(text[i+1])
            i += 1
        else:
            text2.append(text[i])
        
        i += 1
    
    return "".join(text2)


QuotedString = TableType(stringExpand, stringQuote)
TableBool = TableType(manoli_str2bool, str)
TableStr = TableType(str, str)
TableFloat = TableType(float, str)


_defaultTypeLookup = \
    TableTypeLookup([["unknown", str], # backwards compatiable name                     
                     ["str",    str],   # backwards compatiable name
                     ["string", str],
                     ["int",    int],
                     ["float",  float],
                     ["float",  TableFloat],
                     ["bool",   bool],
                     ["bool",   TableBool],
                     ["quoted_string", QuotedString]])

# NOTE: ordering of name-type pairs is important

      
#===========================================================================
# Table class
#

class Table (list):
    """Class implementing the Portable Table Format"""

    def __init__(self, rows=[], 
                       headers=None,
                       defaults={},
                       types={},
                       extraHeaders = [],
                       filename=None,
                       typeLookup=None):
        
        # set table info
        self.headers = copy.copy(headers)
        self.defaults = copy.copy(defaults)
        self.types = copy.copy(types)
        self.extraHeaders = []
        self.filename = filename
        self.comments = []
        self.delim = "\t"
        self.nheaders = 1
        self.version = TABLE_VERSION
        
        if typeLookup == None:
            self.typeLookup = _defaultTypeLookup
        else:
            self.typeLookup = _defaultTypeLookup.extend(typeLookup)
        
        
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
                # guess any types not specified
                if key not in self.types:
                    self.types[key] = type(self[0][key])
            
                # guess any defaults not specified
                if key not in self.defaults:
                    self.defaults[key] = self.types[key]()
        
        # set extra headers
        if len(extraHeaders) > 0:
            for row in extraHeaders:
                self.extraHeaders.append(copy.copy(row))
            
                
                    
    
    
    def new(self, headers=None):
        """
        return a new table with the same info but no data
            
        headers - if specified, only a subset of the headers will be copied
        """
        
        if headers == None:
            headers = self.headers
        
        tab = type(self)(headers=headers)
        
        tab.types = util.subdict(self.types, headers)
        tab.defaults = util.subdict(self.defaults, headers)
        tab.extraHeaders = []
        for row in self.extraHeaders:
            tab.extraHeaders.append(util.subdict(row, headers))
        
        tab.comments = copy.copy(self.comments)
        tab.delim = self.delim
        tab.nheaders = self.nheaders
        
        tab.typeLookup = copy.copy(self.typeLookup)
        
        return tab
    
    
    #===================================================================
    # Input/Output
    #
    
    def read(self, filename, delim="\t", nheaders=1,
                   headers=None, types=None):
        for row in self.readIter(filename, delim=delim, nheaders=nheaders,
                                 headers=headers, types=types):
            self.append(row)
            
    
    def readIter(self, filename, delim="\t", nheaders=1,
                 headers=None, types=None):
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
        self.headers = copy.copy(headers)
        if types == None:
            self.types = {}
        else:
            self.types = copy.copy(types)
        self.defaults = {}
        self.extraHeaders = []        
        self.comments = []
        self.delim = delim
        self.nheaders = nheaders
        self.version = TABLE_VERSION
        
        
        
        # temps for reading only
        self.tmptypes = None
        self.tmpdefaults = None
        self.tmpextraHeaders = []
        
        # line number for error reporting
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

                # split row into tokens
                tokens = line.split(delim)
                
                # if no headers read yet, use this line as a header
                if not self.headers:
                    # parse header
                    if self.parseHeader(tokens):
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
                
                # return completed row
                yield row
                
                
        except Exception, e:
            # report error in parsing input file
            e = TableException(str(e), self.filename, lineno)
            raise e
        
        
        # now that we know the headers we can process extra headers
        for i, row in enumerate(self.tmpextraHeaders):
            if len(row) != len(self.headers):
                raise TableException("wrong number of columns in extra header %d" % i,
                                     self.filename)
            self.extraHeaders.append(dict(zip(self.headers, row)))
        
        
        # clear temps
        del self.tmptypes
        del self.tmpdefaults
        del self.tmpextraHeaders
        
        raise StopIteration
    
    
    def parseHeader(self, tokens):
        if self.nheaders > 0:
            # determine if extra headers exist
            if len(self.tmpextraHeaders) >= self.nheaders - 1:
                self.headers = tokens

                # check that headers are unique
                check = set()
                for header in self.headers:
                    if header in check:
                        raise TableException("Duplicate header '%s'" % header)
                    check.add(header)
            else:
                # this line is an extra header
                self.tmpextraHeaders.append(tokens)
                return True
        else:
            # default headers are numbers
            self.headers = range(len(tokens))
        
        
        # if we have headers then we can initialize other 
        # things (e.g. types, defaults)
        if self.headers:
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

        # return True if this line has been used as a (extra) header
        # return False if this line should be parsed as data
        if self.nheaders > 0:
            return True
        else:
            return False
    
    
    
    def write(self, filename=sys.stdout, delim="\t"):
        """Write a table to a file or stream.
           
           If 'filename' is a string it will be opened as a file.
           If 'filename' is a stream it will be written to directly.
        """
        
        # remember filename for later saving
        if isinstance(filename, str):
            self.filename = filename
    
        out = util.openStream(filename, "w")
        
        self.writeHeader(out, delim=delim)
        
        # tmp variable
        types = self.types
        
        # write data
        for row in self:
            # code is inlined here for speed
            rowstr = []
            for header in self.headers:
                if header in row:
                    rowstr.append(types[header].__str__(row[header]))
                else:
                    rowstr.append('')
            print >>out, delim.join(rowstr)
    
    
    def writeHeader(self, out=sys.stdout, delim="\t"):
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
    
    
    def writeRow(self, out, row, delim="\t"):
        rowstr = []
        types = self.types
        for header in self.headers:
            if header in row:
                rowstr.append(types[header].__str__(row[header]))
            else:
                rowstr.append('')
        out.write(delim.join(rowstr))
        out.write("\n")
    
    
    def save(self):
        """Writes the table to the last used filename for the read() or write()
           function"""
        
        if self.filename != None:
            self.write(self.filename)
        else:
            raise Exception("Table has no filename")
    
    
    #===================================================================
    # Input/Output: Directives
    #
    
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
            self.tmptypes = self.typeLookup.parseTableTypes(rest, self.delim)
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
                      self.typeLookup.formatTableTypes(
                            util.mget(self.types, self.headers),
                            delim) + "\n")
        elif line == DIR_DEFAULTS:
            out.write("##defaults:" +
                      delim.join(map(str, 
                                util.mget(self.defaults, self.headers))) + "\n")
        
        elif line == DIR_HEADERS:
            out.write("##headers:%d\n" % self.nheaders)
        
        else:
            raise "unknown directive:", line
    

    #===================================================================
    # Table manipulation
    #
    
    def add(self, **kargs):
        """Add a row to the table
           
           tab.add(col1=val1, col2=val2, col3=val3)
        """
        self.append(kargs)
    
    
    def addCol(self, header, coltype=None, default=NULL, pos=None, data=None):
        """Add a column to the table.  You must populate column data yourself.
        
           header  - name of the column
           coltype - type of the values in that column
           default - default value of the column
           pos     - position to insert column (default: right-end)
        """
        # ensure header is unique
        if header in self.headers:
            raise Exception("header '%s' is already in table" % header)
        
        # default column position is last column
        if pos == None:
            pos = len(self.headers)
        
        # default coltype is guessed from data
        if coltype == None:
            if data == None:
                raise Exception("must specify data or coltype")
            else:
                coltype = type(data[0])
        
        # default value is inferred from column type
        if default == NULL:
            default = coltype()
        
        # update table info
        self.headers.insert(pos, header)
        self.types[header] = coltype
        self.defaults[header] = default
        
        for row in self.extraHeaders:
            row[header] = ''
        
        # add data
        if data != None:
            for i in xrange(len(self)):
                self[i][header] = data[i]
    
    
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
        
        if col == -1:
            raise Exception("column '%s' is not in table" % oldname)
        
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
            
       
    def getMatrix(self, rowheader="rlabels"):
        """Returns mat, rlabels, clabels
        
           where mat is a copy of the table as a 2D list
                 rlabels are the row labels
                 clabels are the column labels
        """
        
        # get labels
        if rowheader != None and rowheader in self.headers:
            rlabels = self.cget(rowheader)
            clabels = copy.copy(self.headers)
            clabels.remove(rowheader)
        else:
            rlabels = range(len(self))
            clabels = copy.copy(self.headers)

        # get data
        mat = []
        for row in self:
            mat.append(util.mget(row, clabels))
        
        return mat, rlabels, clabels
    
    
    def filter(self, cond):
        """Returns a table with a subset of rows such that cond(row) == True"""
        tab = self.new()
        
        for row in self:
            if cond(row):
                tab.append(row)
        
        return tab
    
    
    def groupby(self, key=None):
        """Groups the row of the table into separate tables based on the 
           function key(row).  Returns a dict where the keys are the values
           retruned from key(row) and the values are tables.
           
           Ex:
           tab = Table([{'name': 'matt', 'major': 'CS'},
                        {'name': 'mike', 'major': 'CS'},
                        {'name': 'alex', 'major': 'bio'}])
           lookup = tab.groupby(lambda x: x['major'])
           
           lookup ==> {'CS': Table([{'name': 'matt', 'major': 'CS'},
                                    {'name': 'mike', 'major': 'CS'}]),
                       'bio': Table([{'name': 'alex', 'major': 'bio'}])}
            
           Can also use a column name such as:
           tab.groupby('major')
            
        """
           
           
        groups = {}
        
        if isinstance(key, str):
            keystr = key
            key = lambda x: x[keystr]
        
        if key == None:
            raise Exception("must specify keyfunc")
        
        
        for row in self:
            key2 = key(row)
            
            # add new table if necessary
            if key2 not in groups:
                groups[key2] = self.new()
            
            groups[key2].append(row)
        
        return groups
    
    
    def lookup(self, *keys, **options):
        """Returns a lookup dict based on a column 'key'
           or multiple keys
           
           extra options:
           default=None
           uselast=False    # allow multiple rows, just use last
        """
        
        options.setdefault("default", None)
        options.setdefault("uselast", False)
        lookup = util.Dict(dim=len(keys), default=options["default"])
        uselast = options["uselast"]
        
        for row in self:
            keys2 = util.mget(row, keys)
            ptr = lookup
            for i in xrange(len(keys2) - 1):
                ptr = lookup[keys2[i]]
            if not uselast and keys2[-1] in ptr:
                raise Exception("duplicate key '%s'" % str(keys2[-1]))
            ptr[keys2[-1]] = row
        
        lookup.insert = False
        return lookup
    
    
    def get(self, rows=None, cols=None):
        """Returns a table with a subset of the rows and columns"""
        
        # determine rows and cols
        if rows == None:
            rows = range(len(self))
        
        if cols == None:
            cols = self.headers
            
        tab = self.new(cols)
        
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
        self.writePretty(s)
        return s.getvalue()
    
    
    def writePretty(self, out=sys.stdout, spacing=2):
        mat2, rlabels, clabels = self.getMatrix(rowheader=None)
        
        # get extra headers
        mat = []
        for row in self.extraHeaders:
            mat.append(util.mget(row, clabels))
        
        # get headers
        mat.append(clabels)
        
        # get data
        mat.extend(mat2)
        
        util.printcols(mat, spacing=spacing, out=out)
    
    def __str__(self):
        s = StringIO.StringIO("w")
        self.write(s)
        return s.getvalue()
    
    



#===========================================================================
# convenience functions
#

def readTable(filename, delim="\t", nheaders=1, typeLookup=None):
    """Read a Table from a file written in PTF"""
    
    table = Table(typeLookup=typeLookup)
    table.read(filename, delim=delim, nheaders=nheaders)
    return table


def iterTable(filename, delim="\t", nheaders=1):
    """Iterate through the rows of a Table from a file written in PTF"""
    
    table = Table()
    return table.readIter(filename, delim=delim, nheaders=nheaders)


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


def joinTables(* args, **kwargs):
    """Join together tables into one table.
       Each argument is a tuple (table_i, key_i, cols_i)
       
       key_i is either a column name or a function that maps a 
       table row to a unique key
    """
    
    if len(args) == 0:
        return Table()
    
    # determine common keys
    tab, key, cols = args[0]
    if isinstance(key, str):
        keys = tab.cget(key)
        lookups = [tab.lookup(key)]        
    else:
        keys = map(key, tab)
        lookup = {}
        for row in tab:
            lookup[key(row)] = row
        lookups = [lookup]
        
    keyset = set(keys)
    

    for tab, key, cols in args[1:]:
        if isinstance(key, str):
            keyset = keyset & set(tab.cget(key))
            lookups.append(tab.lookup(key))            
        else:
            keyset = keyset & set(map(key, tab))
            lookup = {}
            for row in tab:
                lookup[key(row)] = row
            
            lookups.append(lookup)
    
    keys = filter(lambda x: x in keyset, keys)
    
    
    # build new table
    if "headers" not in kwargs:
        headers = util.concat(*util.cget(args, 2))
    else:
        headers = kwargs["headers"]
    tab = Table(headers=headers)
    
    for key in keys:
        row = {}
        for (tab2, key2, cols), lookup in zip(args, lookups):
            row.update(util.subdict(lookup[key], cols))
        tab.append(row)
    
    return tab
               


#===========================================================================
# Matrix functions
#

def matrix2table(mat, rlabels=None, clabels=None, rowheader="rlabels"):
    """
    convert a matrix into a table
    
    use table.getMatrix()  to convert back to a matrix
    
    """
    
    if clabels == None:
        clabels = range(len(mat[0]))
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


def writeMatrix(filename, mat, rlabels=None, clabels=None, rowheader="rlabels"):
    tab = matrix2table(mat,
                       rlabels=rlabels,
                       clabels=clabels,
                       rowheader=rowheader)
    tab.write(filename)


def readMatrix(filename, rowheader="rlabels"):
    tab = readTable(filename)    
    mat, rlabels, clabels = tab.getMatrix(rowheader=rowheader)
    return mat, rlabels, clabels
    


#===========================================================================
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


    #################################################
    # timing
    if 0:
        from rasmus import util
        
        text=["##types:" + "int\t" * 99 + "int",
              "\t".join(map(str, range(100))) ]

        for i in range(10000):
            text.append("1\t" * 99 + "1")
        text = "\n".join(text)
        
        stream = StringIO.StringIO(text)
        
        util.tic("read table")
        tab = readTable(stream)
        util.toc()
    
    
    #################################################
    # specialized types
    if 1:
        text="""\
##types:str	int	strand_type
name	num	strand
matt	123	+
alex	456	-
mike	789	+
john	0	+
"""
        
        
       
        
        class strand_type:
            def __init__(self, text=None):
                if text == None:
                    self.val = True
                else:
                    if text == "+":
                        self.val = True
                    elif text == "-":
                        self.val = False
                    else:
                        raise Exception("cannot parse '%s' as strand_type" % 
                                        str(text))
                
            
            def __str__(self):
                if self.val:
                    return "+"
                else:
                    return "-"
        

        def strand_parser(text=None):
            if text == None:
                return True
            else:
                if text == "+":
                    return True
                elif text == "-":
                    return False
                else:
                    raise Exception("cannot parse '%s' as strand_type" % 
                                    str(text))
        
        def strand_formatter(val):
            if val:
                return "+"
            else:
                return "-"
        
        strand_type = TableType(strand_parser, strand_formatter)
                    

        stream = StringIO.StringIO(text)
        tab = readTable(stream, typeLookup=[["strand_type", strand_type]])
        print tab.types
        print tab
    
    #################################################
    # quoted strings
    if 1:
        text=\
r"""##types:str	bool	quoted_string
name	num	blah
matt	True	hello\tthere
alex	False	hello\nthere
mike	True	hello\\there
john	False	hello\n\\\nthere
"""                    

        stream = StringIO.StringIO(text)
        tab = readTable(stream)
        print tab.types
        print tab
    
    
    #################################################
    # python data structures/code
    if 1:
        def eval2(text=None):
            if text == None:
                return None
            else:
                return eval(text)
        
        python_type = TableType(eval2, str)
    
    
    
        tab = Table(headers=["name", "list"],
                    types={"list": python_type},
                    typeLookup=[["python", python_type]])
                    
        
        tab.append({"name": "matt", "list": [1,2,3]})
        tab.append({"name": "mike", "list": [4,5,6]})
        tab.append({"name": "alex", "list": [7,8,9]})
        
        tab.write()
    
    ##################################################
    # join tables
    if 1:
        tab1 = Table([[0, 1, 2],
                      [1, 3, 4],
                      [2, 5, 6],
                      [3, 7, 8]],
                     headers=['a', 'b', 'c'])
        tab2 = Table([[0, 6, 6],
                      [1, 7, 7],
                      [3, 8, 8]],
                     headers=['a2', 'b2', 'c2'])
        
        tab3 = joinTables((tab1, lambda x: x['a']+1, ['c', 'b']), (tab2, 'a2', ['b2']))
        
        print tab3
    

"""
 file: util.py 
 authors: Matt Rasmussen
 date: 11/30/05

 Common utilities (timers, progress bars, plotting, list/dict manipulation)
"""

# python libs
import copy
import cPickle
import math
import os
import re
import shutil
import sys
import time
import traceback
import urllib2



#
# see bottom of file for other imports
#

INF = 1e1000



def GLOBALS():
    """Maintains a global set of variables"""

    # install globals root if needed
    if not "RASMUS_GLOBALS" in globals():
        globals()["RASMUS_GLOBALS"] = {}
    return globals()["RASMUS_GLOBALS"]

    


class Closure:
    """
    A small class for creating a closure of variables
    handy for nested functions that need to assign to variables in an 
    outer scope

    Example:

    def func1():
        this = Closure(var1 = 0, var2 = "hello")
        def func2():
            this.var1 += 1
        func2()
        print this.var1   
    func1()
    
    will produce:
    1
    
    """

    def __init__(self, **variables):
        for key, val in variables.iteritems():
            setattr(self, key, val)



class Dict (dict):
    """My personal nested Dictionary (with default values)"""
    
    def __init__(self, dim=1, default=None):
        """dim     - number of dimensions of the dictionary
           default - default value of a dictionary item
        """
        self.dim = dim
        self.null = default
        
        # backwards compatiability
        self.data = self
    
    def __getitem__(self, i):
        if not i in self:
            if self.dim > 1:
                self[i] = Dict(self.dim - 1, self.null)
            else:
                self[i] = copy.copy(self.null)
                return self[i]
        return dict.__getitem__(self, i)
    
    def __len__(self):
        if self.dim == 1:
            return dict.__len__(self)
        else:
            nnz = 0
            for i in self:
                nnz += len(self[i])
            return nnz
        
    def has_keys(self, *keys):
        if len(keys) == 0:
            return True
        elif len(keys) == 1:
            return dict.has_key(self, keys[0])
        else:
            return dict.has_key(self, keys[0]) and \
                   self[keys[0]].has_keys(*keys[1:])
    
    def write(self, out = sys.stdout):
        def walk(node, path):
            if node.dim == 1:
                for i in node:
                    print >>out, "  ",
                    for j in path:
                        print str(j) + ", ",
                    print >>out, i, ":", node[i]
            else:
                for i in node:
                    walk(node[i], path + [i])
        
        print >>out, "< DictMatrix "
        walk(self, [])
        print >>out, ">"


class Percent (float):
    digits = 1
    
    def __str__(self):
        return (("%%.%df" % self.digits) % (float(self) * 100))
    
    def __repr__(self):
        return str(self)


#################################################################################
# list and dict functions
#

def equal(* vals):
    """Returns True if all arguments are equal"""
    if len(vals) < 2:
        return True
    a = vals[0]
    for b in vals[1:]:
        if a != b:
            return False
    return True


def remove(lst, *vals):
    """Returns a copy of list 'lst' with values 'vals' removed
    """
    lst2 = []
    delset = makeset(vals)
    for i in lst:
        if i not in delset:
            lst2.append(i)
    return lst2


def sort(lst, compare=cmp):
    """Returns a sorted copy of a list
       
       arguments:
       lst     - a list to sort
       compare - a function for comparing items (default: cmp)
    """
    lst2 = copy.copy(lst)
    lst2.sort(compare)
    return lst2


def reverse(lst):
    """Returns a reversed copy of a list
    """
    lst2 = copy.copy(lst)
    lst2.reverse()
    return lst2


def cget(mat, *i):
    """Returns the column(s) '*i' of a 2D list 'mat'
    
       notes:
       If one column is given, the column is returned as a list.
       If multiple columns are given, a submatrix is returned.
    """
    if len(i) == 1:
        return [row[i[0]] for row in mat]
    else:
        return [[row[j] for j in i] for row in mat]


def mget(data, keys):
    """Returns a list 'lst2' such that lst2[i] = lst[ind[i]]
       
       Or in otherwords, get the subsequence 'lst'
       Note: same as sublist.  This name is being introduced, because sublist
       sounds like it can't be used on dicts; but it can.
       
    """
    lst = []
    for key in keys:
        lst.append(data[key])
    return lst


def concat(* lists):
    """Concatenates several lists into one
    """
    
    lst = []
    for l in lists:
        lst.extend(l)
    return lst


def sublist(lst, ind):
    """Returns a list 'lst2' such that lst2[i] = lst[ind[i]]
       
       Or in otherwords, get the subsequence 'lst'
       same as mget
       
       DEPRECATED
    """
    lst2 = []
    for i in ind:
        lst2.append(lst[i])
    return lst2


def subdict(dic, keys):
    """Returns a new dictionary dic2 such that
       dic2[i] = dic[i] for all i in keys
       
       arguments:
       dic  - a dictionary
       keys - a list of keys
    """
    dic2 = {}
    for key in keys:
        if key in dic:
            dic2[key] = dic[key]
    return dic2


def revdict(dic, allowdups=False):
    """Reverses a dict 'dic' such that the keys become values and the 
       values become keys.
    """
    
    dic2 = {}
    if allowdups:
        for key, val in dic.iteritems():
            dic2[val] = key
    else:
        for key, val in dic.iteritems():
            assert key not in dic2, "duplicate value '%s' in dict" % val
            dic2[val] = key
    
    return dic2


def list2lookup(lst):
    """Creates a dict where each key is lst[i] and value is i
    """
    
    lookup = {}
    for i in range(len(lst)):
        lookup[lst[i]] = i
    return lookup


def list2dict(lst, val=1):
    """Creates a dict with keys from list 'lst' and values 'val'
    """
    
    dic = {}
    for i in lst:
        dic[i] = val
    return dic


def items2dict(items):
    """Creates a dict from a list of key, value pairs
    
    TODO: this is redundant with dict(items).  Should remove
    """
    
    dic = {}
    for key, val in items:
        dic[key] = val
    return dic


def mapdict(dic, keyfunc=lambda x:x, valfunc=lambda x:x):
    """Creates a new dict where keys and values are mapped
    """
    
    dic2 = {}
    for key, val in dic.iteritems():
        dic2[keyfunc(key)] = valfunc(val)
    
    return dic2


# TODO: make lst a *lst and use zip(lst)
#
def groupby(func, lst):
    """Places i and j of 'lst' into the same group if func(i) == func(j).
       
       arguments:
       func - is a function of one argument that maps items to group objects
       lst  - is a list of items
       
       returns:
       a dictionary such that the keys are groups and values are items found in
       that group
    """
    
    div = {}
    for i in lst:
        lst2 = div.setdefault(func(i), [])
        lst2.append(i)
    
    return div


def unique(lst):
    """Returns a copy of 'lst' with only unique entries
    """
    
    return list2dict(lst).keys()


def flatten(lst, depth=INF):
    """Flattens nested lists into one list
    """
    
    flat = []
    
    for i in lst:
        if isinstance(i, list) and depth > 0:
            flat.extend(flatten(i, depth-1))
        else:
            flat.append(i)
    
    return flat


def mapapply(funcs, lst):
    """apply each function in funcs to one element in lst
    """
    
    lst2 = []
    for func, item in zip(funcs, lst):
        lst2.append(func(item))
    return lst2


def map2(func, matrix):
    """Maps a function onto the elements of a matrix
    """
    
    return [map(func, row) for row in matrix]


def min2(matrix):
    """Finds the minimum of a 2D list or matrix
    """
    return min(map(min, matrix))


def max2(matrix):
    """Finds the maximum of a 2D list or matrix
    """
    return max(map(max, matrix))


class range2:
    def __init__(self, width, height):
        self.width = width
        self.height = height
        self.i = -1
        self.j = 0
    
    def __iter__(self):
        return self
    
    def next(self):
        self.i += 1
        if self.i >= self.width:
            self.j += 1
            self.i = 0
        
        if self.j >= self.height:
            raise StopIteration
        
        return self.i, self.j
    

def frange(start, end, step):
    """Returns a range of floats 
    
       arguments:
       start - begining of range
       end   - end of range
       step  - step size
    """
    
    i = 0
    lst = []
    while start < end:
        lst.append(start + i * step)
        i += 1
    return lst


def count(func, lst):
    """Counts the number of times func(x) is True for x in list 'lst'
       
       See also:
        counteq(a, lst)   count items equal to a
        countneq(a, lst)  count items not equal to a
        countle(a, lst)   count items less than or equal to a
        countlt(a, lst)   count items less than a
        countge(a, lst)   count items greater than or equal to a
        countgt(a, lst)   count items greater than a
    """
    n = 0
    for i in lst:
        if func(i):
            n += 1
    return n

def counteq(a, lst): return count(eqfunc(a), lst)
def countneq(a, lst): return count(neqfunc(a), lst)
def countle(a, lst): return count(lefunc(a), lst)
def countlt(a, lst): return count(ltfunc(a), lst)
def countge(a, lst): return count(gefunc(a), lst)
def countgt(a, lst): return count(gtfunc(a), lst)


def find(func, lst):
    """Returns the indices of func(x) == True for x in list 'lst'"""
    pos = []
    for i in xrange(len(lst)):
        if func(lst[i]):
            pos.append(i)
    return pos

def findeq(a, lst): return find(eqfunc(a), lst)
def findneq(a, lst): return find(neqfunc(a), lst)
def findle(a, lst): return find(lefunc(a), lst)
def findlt(a, lst): return find(ltfunc(a), lst)
def findge(a, lst): return find(gefunc(a), lst)
def findgt(a, lst): return find(gtfunc(a), lst)



def islands(lst):
    """
      takes a list, or a string, and gather islands of identical elements.
      it returns a dictionary counting where
      counting = {element: [(start,end), (start,end), ...],
                  element: [(start,end), (start,end), ...],
                  ...}
      counting.keys() is the list of unique elements of the input list
      counting[element] is the list of all islands of occurence of element
      counting[element][i] = (start,end)
       is such that list[start-1:end] only contains element
       
       note: this comes from Manolis's code
    """
    
    if not lst: return {}

    counting = {}

    i,current_char, current_start = 0,lst[0], 0
    
    while i < len(lst):

        if current_char == lst[i]:
            i = i+1
        else:
            if not counting.has_key(current_char): counting[current_char] = []
            counting[current_char].append((current_start, i))
            current_char = lst[i]
            current_start = i

    if not counting.has_key(current_char): counting[current_char] = []
    counting[current_char].append((current_start, i))

    return counting


def overlap(a, b, x, y, inc=True):
    """Returns True if range [a,b] overlaps [x,y]
    """
    if inc:
        return (y >= a) and (x <= b)
    else:
        return (y > a) and (x < b)


def argmax(lst):
    if len(lst) == 0:
        return -1
    top = 0
    for i in xrange(1, len(lst)):
        if lst[i] > lst[top]:
            top = i
    return top

def argmin(lst):
    if len(lst) == 0:
        return -1
    low = 0
    for i in xrange(1, len(lst)):
        if lst[i] < lst[low]:
            low = i
    return low

def maxfunc(func, lst):
    top = -INF
    for i in lst:
        val = func(i)
        if val > top:
            top = val
    return top

def minfunc(func, lst):
    low = -INF
    for i in lst:
        val = func(i)
        if val < low:
            low = val
    return low

def argmaxfunc(func, lst, index=True):
    if len(lst) == 0:
        return -1
    top = 0
    topval = func(lst[top])
    for i in xrange(1,len(lst)):
        val = func(lst[i])
        if val > topval:
            topval = val
            top = i
    if index:
        return top
    else:
        return lst[top]
    
def argminfunc(func, lst, index=True):
    if len(lst) == 0:
        return -1
    low = 0
    lowval = func(lst[low])
    for i in xrange(1, len(lst)):
        val = func(lst[i])
        if val < lowval:
            lowval = val
            low = i
    if index:
        return low
    else:
        return lst[low]


def eqfunc(a): return lambda x: x == a
def neqfunc(a): return lambda x: x != a
def ltfunc(a): return lambda x: x < a
def gtfunc(a): return lambda x: x > a
def lefunc(a): return lambda x: x <= a
def gefunc(a): return lambda x: x >= a
def withinfunc(a, b, ainc=True, binc=True):
    if ainc:
        if binc:
            return lambda x: a <= x <= b
        else:
            return lambda x: a <= x < b
    else:
        if binc:
            return lambda x: a < x <= b
        else:
            return lambda x: a < x < b


#
# math functions
#

def sign(num):
    """Returns the sign of a number"""
    return cmp(num,0)

def lg(num):
    """Retruns the log_2 of a number"""
    return math.log(num, 2)

def add(a, b): return a + b
def sub(a, b): return a - b
def mul(a, b): return a * b
def idiv(a, b): return a / b
def div(a, b): return a / float(b)

def safediv(a, b, default=INF):
    try:
        return a / float(b)
    except ZeroDivisionError:
        return default

def safelog(x, base=math.e, default=-INF):
    try:
        return math.log(x)
    except OverflowError:
        return default
        
def invcmp(a, b): return cmp(b, a)

def clamp(x, low, high):
    """Clamps a value 'x' between the values 'low' and 'high'"""
    
    if high != None and x > high:
        return high
    elif low != None and x < low:
        return low
    else:
        return x    

def clampfunc(low, high):
    return lambda x: clamp(x, low, high)


def compose(* funcs):
    """Composes two or more functions into one function
    
       example:
       compose(f,g)(x) <==> f(g(x))
    """
    
    def compose2(f, g):
        return lambda x: f(g(x))
    
    f = funcs[-1]
    for i in xrange(len(funcs)-2, -1, -1):
        f = compose2(funcs[i], f)
    return f

#
# set operations
#
def makeset(lst):
    return list2dict(lst)

def union(set1, set2):
    set = {}
    set.update(set1)
    set.update(set2)
    return set

def intersect(set1, set2):
    set = {}
    for i in set1:
        if i in set2:
            set[i] = 1
    return set

def nonintersect(set1, set2):
    diffset = setSubtract(set1, set2)
    diffset.update(setSubtract(set2, set1))
    return diffset

def setSubtract(set1, set2):
    set = {}
    for i in set1:
        if i not in set2:
            set[i] = 1
    return set


#
# regex
#

def match(pattern, line):
    """  remember: to name tokens use (?P<name>pattern) """
    
    m = re.match(pattern, line)
    
    if m == None:
        return {}
    else:
        return m.groupdict()


###############################################################################
# common Input/Output

def readInts(filename):
    """Read a list of integers from a file (one int per line)
    
       filename may also be a stream
    """
    
    infile = openStream(filename)
    vec = []
    for line in infile:
        vec.append(int(line))
    return vec

def readFloats(filename):
    """Read a list of floats from a file (one float per line)
    
       filename may also be a stream
    """
    infile = openStream(filename)
    vec = []
    for line in infile:
        vec.append(float(line))
    return vec

def readStrings(filename):
    """Read a list of strings from a file (one string per line)
    
       filename may also be a stream
    """
    infile = openStream(filename)
    vec = [line.rstrip() for line in infile]
    return vec

def writeVector(filename, vec):
    """Write a list of anything (ints, floats, strings, etc) to a file.
    
       filename may also be a stream
    """
    out = openStream(filename, "w")
    for i in vec:
        print >>out, i


def openStream(filename, mode = "r"):
    """Returns a file stream depending on the type of 'filename' and 'mode'
    
       The following types for 'filename' are handled:
       
       stream         - returns 'filename' unchanged
       URL string     - opens http pipe
       '-'            - opens stdin or stdout, depending on 'mode'
       other string   - opens file with name 'filename'
       
       mode is standard mode for file(): r,w,a
    """
    
    # if filename has a file interface then return it back unchanged
    if "read" in dir(filename) or \
       "write" in dir(filename):
        return filename
    
    # if filename is a string then open it
    elif isinstance(filename, str):
        # open URLs
        if filename.startswith("http://"):
            return urllib2.urlopen(filename)
        
        # open stdin and stdout
        elif filename == "-":
            if mode == "w":
                return sys.stdout
            elif mode == "r":
                return sys.stdin
            else:
                raise Exception("stream '-' can only be opened with modes r/w")
        
        # open regular file
        else:
            return file(filename, mode)
    
    # cannot handle other types for filename
    else:
        raise Exception("unknown filename type '%s'" % type(filename))





#################################################################################
# printting functions
#

def defaultJustify(val):
    if isinstance(val, int) or \
       isinstance(val, float):
        return "right"
    else:
        return "left"


def defaultFormat(val):
    if isinstance(val, int) and \
       not isinstance(val, bool):
        return int2pretty(val)
    elif isinstance(val, Percent):
        return str(val)
    elif isinstance(val, float):
        return "%.3f" % val
    else:
        return str(val)


def printcols(data, width=None, spacing=1, format=defaultFormat, 
              justify=defaultJustify, out=sys.stdout):
    """Prints a list or matrix in aligned columns
        
       data    - a list or matrix
       width   - maxium number of characters per line (default: 75 for lists)
       spacing - number of spaces between columns (default: 1)
       out     - stream to print to (default: sys.stdout)
    """
    
    if len(data) == 0:
        return
    
    if isinstance(data[0], list) or \
       isinstance(data[0], tuple):
        # matrix printing has default width of unlimited
        if width == None:
            width = 100000
        
        mat = data
    else:
        # list printing has default width 75
        if width == None:
            width = 75
        
        ncols = int(width / (max(map(lambda x: len(str(x)), data))+ spacing))
        mat = list2matrix(data, ncols=ncols, bycols=True)
    
    
    # turn all entries into strings
    matstr = map2(format, mat)
    
    
    # ensure every row has same number of columns
    maxcols = max(map(len, matstr))
    for row in matstr:
        if len(row) < maxcols:
            row.extend([""] * (maxcols - len(row)))
    
    
    # find the maximum width char in each column
    maxwidths = map(max, map2(len, zip(* matstr)))
    
    
    # print out matrix with whitespace padding
    for i in xrange(len(mat)):
        fields = []
        for j in xrange(len(mat[i])):
            just = justify(mat[i][j])
            
            if just == "right":
                fields.append((" " * (maxwidths[j] - len(matstr[i][j]))) + \
                              matstr[i][j] + \
                              (" " * spacing))
            else:
                # do left by default            
                fields.append(matstr[i][j] + 
                              (" " * (maxwidths[j] - len(matstr[i][j]) + spacing)))
        out.write("".join(fields)[:width] + "\n")


def list2matrix(lst, nrows=None, ncols=None, bycols=True):
    """Turn a list into a matrix by wrapping its entries"""
    
    mat = []
    
    if nrows == None and ncols == None:
        nrows = int(math.sqrt(len(lst)))
        ncols = int(math.ceil(len(lst) / nrows))
    elif nrows == None:
        nrows = int(math.ceil(len(lst) / min(ncols, len(lst))))
    else:
        ncols = int(math.ceil(len(lst) / min(nrows, len(lst))))

    for i in xrange(nrows):
        mat.append([])
        for j in xrange(ncols):
            if bycols:
                k = i + j*nrows
            else:
                k = i*ncols + j
            if k < len(lst):
                mat[-1].append(lst[k])
    
    return mat


def printwrap(text, width=80, prefix="", out=sys.stdout):
    """Prints text with wrapping"""
    if width == None:
        out.write(text)
        out.write("\n")
        return
    
    pos = 0
    while pos < len(text):
        out.write(prefix)
        out.write(text[pos:pos+width])
        out.write("\n")
        pos += width


def printDict(dic, keyfunc=lambda x: x, valfunc=lambda x: x,
              num=None, compare=lambda a,b: cmp(a[0],b[0]),
              spacing=4, out=sys.stdout,
              format=defaultFormat, 
              justify=defaultJustify):
    if num == None:
        num = len(dic)
    
    dic = mapdict(dic, keyfunc=keyfunc, valfunc=valfunc)
    items = dic.items()
    items.sort(compare)
    
    printcols(items[:num], spacing=spacing, out=out, format=format, 
              justify=justify)


def printDictByKeys(dic, keyfunc=lambda x: x, valfunc=lambda x: x,
                    num=None, spacing=4, compare=cmp, out=sys.stdout):
    printDict(dic, keyfunc=keyfunc, valfunc=valfunc, 
              num=num, compare=lambda a,b: compare(a[0],b[0]),
              spacing=spacing, out=out)


def printDictByValues(dic, keyfunc=lambda x: x, valfunc=lambda x: x,
                      num=None, spacing=4, compare=cmp, out=sys.stdout):
    printDict(dic, keyfunc=keyfunc, valfunc=valfunc, 
              num=num, compare=lambda a,b: compare(a[1],b[1]),
              spacing=spacing, out=out)


def printHistDict(array, keyfunc=lambda x: x, valfunc=lambda x: x,
                  num=None, compare=lambda a,b: cmp(b[1],a[1]),
              spacing=4, out=sys.stdout):
    hist = histDict(array)
    printDict(hist, keyfunc=keyfunc, valfunc=valfunc, 
              num=num, compare=compare, spacing=spacing, out=out)


def printHist(array, ndivs=20, width=75, spacing=2, out=sys.stdout):            
    data = list(hist(array, ndivs))                                             
                                                                                
    # find max bar                                                              
    maxwidths = map(max, map2(compose(len, str), data))                         
    maxbar = width - sum(maxwidths) - 2 * spacing                               
                                                                                
    # make bars                                                                 
    bars = []                                                                   
    maxcount = max(data[1])                                                     
    for count in data[1]:                                                       
        bars.append("*" * int(count * maxbar / float(maxcount)))                
    data.append(bars)                                                           
                                                                                
    printcols(zip(* data), spacing=spacing, out=out)   


def int2pretty(num):
    """Returns a pretty-printed version of an int"""
    
    string = str(num)
    parts = []
    l = len(string)
    for i in xrange(0, l, 3):
        t = l - i
        s = t - 3
        if s < 0: s = 0
        parts.append(string[s:t])
    parts.reverse()
    return ",".join(parts)


def pretty2int(string):
    """Parses a pretty-printed version of an int into an int"""
    return int(string.replace(",", ""))


def str2bool(val):
    """Correctly converts the strings "True" and "False" to the 
       booleans True and False
    """
    
    if val == "True":
        return True
    elif val == "False":
        return False
    else:
        raise Exception("unknown string for bool '%s'" % val)

                

class DelimReader:
    """Reads delimited files"""

    def __init__(self, filename, delim=None, header=False, useDict=True):
        """Constructor for DelimReader
            
           arguments:
           filename  - filename or stream to read from
           delim     - delimiting character
           header    - a bool specifying if the first row is a header
        """
        
        self.infile = openStream(filename)
        self.header = header
        self.delim = delim
        self.useDict = useDict
        self.headers = []
        
        if self.header:
            self.headers = self.split(self.infile.next())

    def __iter__(self):
        return self
    
    def next(self):
        line = self.infile.next()
        fields = self.split(line)
        
        if self.header and self.useDict:
            row = {}
            for i,j in zip(self.headers, fields):
                row[i] = j
            return row
        else:
            return fields

    def split(self, line):
        return line.rstrip().split(self.delim)


def readDelim(filename, delim=None, header=False):
    """Read an entire delimited file into memory as a 2D list"""
    
    reader = DelimReader(filename, delim, header)
    data = [row for row in reader]
    return data

def writeDelim(filename, data, delim="\t"):
    """Write a 2D list into a file using a delimiter"""
    
    out = openStream(filename, "w")
    for line in data:
        print >>out, delim.join(map(str, line))


class IterFunc:
    def __init__(self, func):
        self.func = func
    
    def __iter__(self):
        return self
    
    def next(self):
        return self.func()
    

def selcolIter(myiter, cols):       
    return IterFunc(lambda: sublist(myiter.next(), cols))

def joinIter(myiter, delim):
    return IterFunc(lambda: delim.join(myiter.next()))

def cutIter(myiter, cols, delim=None):
    return IterFunc(lambda: delim.join(sublist(myiter.next().split(delim), cols)))



def clearFile(filename):
    out = file(filename, "w")
    out.close()

def filesize(filename):
    return os.stat(filename)[6]

def openZip(filename):
    (infile, outfile) = os.popen2("zcat '"+filename+"' ")
    return outfile


class SafeReadIter:
    def __init__(self, infile):
        self.infile = infile
    
    def __iter__(self):
        return self
    
    def next(self):
        line = self.infile.readline()
        if line == "":
            raise StopIteration
        else:
            return line

def linecount(filename):
    count = 0
    for line in openStream(filename):
        count += 1
    return count

def readWord(infile, delims = [" ", "\t", "\n"]):
    word = ""
    
    while True:
        char = infile.read(1)
        if char == "":
            return word
        if char not in delims:
            word += char
            break
    
    while True:
        char = infile.read(1)
        if char == "" or char in delims:
            return word
        word += char


def readUntil(stream, chars):
    token = ""
    while True:
        char = stream.read(1)
        if char in chars or char == "":
            return token, char
        token += char


def readWhile(stream, chars):
    token = ""
    while True:
        char = stream.read(1)
        if char not in chars or char == "":
            return token, char
        token += char


def skipComments(infile):
    for line in infile:
        if line.startswith("#"):
            continue
        yield line


class IndentStream:
    """
    Makes any stream into an indent stream.
    
    Indent stream auto indents every line written to it
    """
    
    def __init__(self, stream):
        self.stream = openStream(stream, "w")
        self.linestart = True
        self.depth = 0
    
    def indent(self, num=2):
        self.depth += num
    
    def dedent(self, num=2):
        self.depth -= num
        if self.depth < 0:
            self.depth = 0
    
    def write(self, text):
        lines = text.split("\n")
        
        for line in lines[:-1]:
            if self.linestart:
                self.stream.write(" "*self.depth)
                self.linestart = True
            self.stream.write(line + "\n")
        
        if len(lines) > 0:
            if text.endswith("\n"):
                self.linestart = True
            else:
                self.stream.write(" "*self.depth + lines[-1])
                self.linestart = False




    
    
###############################################################################
# file/directory stuff
#
def listFiles(path, extension=""):
    """Returns a list of files in 'path' with ending with 'extension'"""
    if path[-1] != "/":
        path += "/"
    files = filter(lambda x: x.endswith(extension), os.listdir(path))
    files.sort()
    files = map(lambda x: path + x, files)
    return files

def tempfile(dir, prefix, ext):
    """Generates a a temp filename 'dir/prefix_XXXXXX.ext'"""
    
    import warnings
    warnings.filterwarnings("ignore", ".*", RuntimeWarning)
    filename = os.tempnam(dir, "____")      
    filename = filename.replace("____", prefix) + ext
    warnings.filterwarnings("default", ".*", RuntimeWarning)
    
    return filename

def deldir(path):
    """Recursively remove a directory"""
    
    # This function is slightly more complicated because of a 
    # strange behavior in AFS, that creates .__afsXXXXX files
    
    dirs = []
    
    def cleandir(arg, path, names):
        for name in names:
            filename = os.path.join(path, name)
            if os.path.isfile(filename):
                os.remove(filename)
        dirs.append(path)
    
    # remove files
    os.path.walk(path, cleandir, "")
    
    # remove directories
    for i in xrange(len(dirs)):
        # AFS work around
        afsFiles = listFiles(dirs[-i])
        for f in afsFiles:
            os.remove(f)
        
        while True:
            try:
                os.rmdir(dirs[-i])
            except:
                continue
            break


def replaceExt(filename, oldext, newext):
    if filename.endswith(oldext):
        return filename[:-len(oldext)] + newext
    else:
        raise Exception("file '%s' does not have extension '%s'" % (filename, oldext))



###############################################################################
# sorting
#
def sortInd(array, compare = cmp):
    """Returns list of indices into 'array' sorted by 'compare'"""
    ind = range(len(array))
    ind.sort(lambda x, y: compare(array[x], array[y]))
    return ind

def sortPerm(array, compare = cmp):
    """Returns a copy of list 'array' sorted by 'compare' and the 
       resulting permutation"""
    perm = sortInd(array, compare)
    sorted = permute(array, perm)
    return sorted, perm

def sortTogether(compare, array, *others):
    ind = sortInd(array, compare)
    arrays = [permute(array, ind)]
    
    for other in others:
        arrays.append(permute(other, ind))
    
    return arrays

def permute(lst, perm):
    """Returns a copy of list 'lst' permuted by 'perm'"""
    sorted = [0] * len(lst)
    for i in range(len(sorted)):
        sorted[i] = lst[perm[i]]
    return sorted
    
def invPerm(perm):
    """Returns the inverse of a permutation 'perm'"""
    inv = [0] * len(perm)
    for i in range(len(perm)):
        inv[perm[i]] = i
    return inv


###############################################################################
# histograms, distributions
#
def oneNorm(vals):
    """Normalize values so that they sum to 1"""
    s = float(sum(vals))
    return map(lambda x: x/s, vals)

def bucketSize(array, ndivs = 20):
    """Determine the bucket size needed to divide the values in array into 
       'ndivs' evenly sized buckets"""
    #return abs(max(array) - min(array) + 1) / float(ndivs)
    return abs(max(array) - min(array)) / float(ndivs)

def hist(array, ndivs = 20, size=None):
    """Create a histogram of 'array' with 'ndivs' buckets"""
    
    if size != None:
        ndivs = int(abs(max(array) - min(array) + 1) / float(size))
    
    h = [0] * ndivs
    x = []
    bucketwidth = bucketSize(array, ndivs)
    low = min(array)
    for i in array:
        h[min(int((i - low) / bucketwidth), ndivs-1)] += 1
    for i in range(ndivs):
        x.append(i * bucketwidth + low)
    return (x, h)

def bucket(array, ndivs = 20, low=None, bucketwidth=None):
    """Group elements of 'array' into 'ndivs' lists"""
    
    h = [[] for i in range(ndivs)]
    x = []    
    if bucketwidth == None:
        bucketwidth = bucketSize(array, ndivs)
    if low == None:
        low = min(array)
    
    for i in array:
        h[min(int((i - low) / bucketwidth), ndivs-1)].append(i)
    for i in range(ndivs):
        x.append(i * bucketwidth + low)
    return (x, h)

def hist2(array1, array2, ndivs1=20, ndivs2=20):
    """Perform a 2D histogram"""
    
    h = [[0] * ndivs1 for i in xrange(ndivs2)]
    labels = []
    bucketwidth = bucketSize(array1, ndivs1)
    bucketheight = bucketSize(array2, ndivs2)
    low1 = min(array1)
    low2 = min(array2)
    
    for j,i in zip(array1, array2):
        h[min(int((i - low2) / bucketheight), ndivs2-1)] \
         [min(int((j - low1) / bucketwidth), ndivs1-1)] += 1
    
    for i in range(ndivs2):
        labels.append([])
        for j in range(ndivs1):        
            labels[-1].append([j * bucketwidth + low1,
                               i * bucketheight + low2])
    return labels, h
    
  

def distrib(array, ndivs = 20, size=None):
    """Find the distribution of 'array' using 'ndivs' buckets"""
    
    h = hist(array, ndivs, size=size)
    area = 0
    if size != None:
        delta = size
    else:
        delta = bucketSize(array, ndivs)
    total = float(sum(h[1]))
    return (h[0], map(lambda x: (x/total)/delta, h[1]))

def plothist(array, ndivs = 20, **options):
    """Plot a histogram of array"""
    h = hist(array, ndivs)
    p = options.setdefault("plot", Gnuplot())
    options.setdefault("style", "boxes")
    
    p.plot(h[0], h[1], **options)
    return p

def plotdistrib(array, ndivs = 20, **options):
    """Plot a distribution of array"""
    d = distrib(array, ndivs)
    p = options.setdefault("plot", Gnuplot())
    options.setdefault("style", "boxes")
    
    p.plot(d[0], d[1], **options)
    return p

def histInt(array):
    """Returns a histogram of integers as a list of counts"""
    hist = [0]  * (max(array) + 1)
    negative = []
    for i in array:
        if (i >= 0):
            hist[i] += 1
        else:
            negative.append(i)
    return hist

def histDict(array):
    """Returns a histogram of any items as a dict.
    
       The keys of the returned dict are elements of 'array' and the values
       are the counts of each element in 'array'.
    """
    hist = {}
    for i in array:
        if i in hist.keys():
            hist[i] += 1
        else:
            hist[i] = 1
    return hist





# import common functions from other files, 
# so that only util needs to be included
from timer import *
from vector import *
from progress import *
from options import *
from plotting import *
from workspace import *






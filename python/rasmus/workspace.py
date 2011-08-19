###############################################################################
# pickling
#
# TODO: use python sets

import sys
import cPickle


from rasmus import util


# globals
_varset = {}


def binsave(filename, var):
    cPickle.dump(var, util.open_stream(filename, "w"), 2)

def binload(filename):
    return cPickle.load(util.open_stream(filename))

def varset():
    globals _varset
    return _varset

def addvar(* varnames):
    for name in varnames:
        varset()[name] = 1

def delvar(* varnames):
    for name in varnames:
        del varset()[name]

def getvars(table):
    set = util.subdict(table, varset())
    return set

def setvars(table, dct):
    for name in dct:
        table[name] = dct[name]

def showvars(table=None, width=70, out=sys.stdout):
    names = varset().keys()
    names.sort()
    
    if not table:
        util.printcols(names, width, out)
    else:
        maxname = max(map(len, names))
        
        
        for name in names:
            out.write(name + " "*(maxname-len(name)) + "  ")
            out.write("type: " + type(table[name]).__name__)
            
            if "__len__" in dir(table[name]):
                out.write(", size: %d\n" % len(table[name]))
            else:
                out.write("\n")
            


def saveall(table, filename = "all.pickle"):
    binsave(filename, getvars(table))
    
    # display variables saved
    keys = varset().keys()
    keys.sort()
    for key in keys:
        log("%s: saved '%s'" % (filename, key))

def loadall(table, filename = "all.pickle"):
    set = binload(filename)
    setvars(table, set)
    
    # also add new var names to varset
    set2 = varset()
    keys = set.keys()
    keys.sort()
    for key in keys:
        set2[key] = 1
        log("%s: loaded '%s'" % (filename, key))
        
def clearall():
    varset().clear()

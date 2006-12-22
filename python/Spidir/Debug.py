#-------------------------------------------------------------------------------
# debugging variables and functions
#-------------------------------------------------------------------------------

# python libs
import sys
import math

# rasmus libs
from rasmus import util
from rasmus import stats
from rasmus import treelib



# globals
DEBUG = sys.stdout

# constants
DEBUG_NONE = 0
DEBUG_LOW = 1
DEBUG_MED = 2
DEBUG_HIGH = 3


DEBUG_LEVEL = DEBUG_NONE



def setDebug(level=DEBUG_NONE):
    global DEBUG_LEVEL
    DEBUG_LEVEL = level

def isDebug(level):
    return DEBUG_LEVEL >= level


def debug(* text, **args):
    args.setdefault("level", DEBUG_NONE)

    # check debug level
    if DEBUG_LEVEL < args["level"]:
        return

    output = " ".join(map(str, text))
    
    if "nonl" in args:
        DEBUG.write(output)
    else:
        print >>DEBUG, output



def setDebugStream(stream):
    globals()["DEBUG"] = stream


def drawTreeLogl(tree, out=None, events={}, baserate=1.0):
    labels = {}
    
    if out == None:
        out = DEBUG
    
    if "baserate" in tree.data:
        baserate = tree.data["baserate"]
    
    
    for node in tree.nodes.values():
        notes = ""
        if "extra" in node.data:
            notes += "E"
        if "unfold" in node.data:
            notes += "U"
        
        if "logl" in node.data:
            if isinstance(node.data["logl"], float):
                labels[node.name] = "[%s]\n%.3f (%.3f) %s" % \
                    (node.name, node.dist, node.data["logl"], notes)
                #logl += node.data["logl"]
            else:
                labels[node.name] = "[%s]\n%.3f (%s) %s" % \
                    (node.name, node.dist, str(node.data["logl"]), notes)

        else:
            labels[node.name] = "[%s]\n%.3f (*) %s" % \
                (node.name, node.dist, notes)
        
        if "params" in node.data:
            fracs = map(stats.mean, zip(* node.data["fracs"]))
            mean = sum(util.vmul(util.cget(node.data["params"], 0), fracs))
            sdev = sum(util.vmul(util.cget(node.data["params"], 1), fracs))
            
            mean *= baserate
            sdev *= baserate
            
            labels[node.name] += "\n%.3f %.3f" % (mean, sdev)
        
        
        #if "error" in node.data:
        #    labels[node.name] += "\nerr %.4f" % node.data["error"]
        
        if node in events:
            labels[node.name] += " %s" % events[node]
        
    if "logl" in tree.data:
        debug("logl:      %f" % tree.data["logl"])
        debug("eventlogl: %f" % tree.data["eventlogl"])
        debug("errorlogl: %f" % tree.data["errorlogl"])
    debug("baserate:  %f" % baserate)
    debug("treelen:   %f" % sum(x.dist for x in tree.nodes.values()))
    if "error" in tree.data:
        debug("error:     %f" % tree.data["error"])
    
    treelib.drawTree(tree, minlen=20, maxlen=100, labels=labels, spacing=4, 
                        labelOffset=-3, out=out)



class SindirError (Exception):
    def __init__(self, msg):
        Exception.__init__(self)
        self.msg = msg
    def __str__(self):
        return str(self.msg)



def printVisitedTrees(visited):
    if len(visited) == 0:
        return
    nleaves = len(visited.values()[0][1].leaves())
    
    debug("\n\nmost likily trees out of %d visited (%5.1f total): " % \
          (len(visited), numPossibleTrees(nleaves)))
    
    mat = [[key, logl, 
           tree.data["error"], 
           tree.data["baserate"],
           count]
           for key, (logl, tree, count) in visited.iteritems()]
    mat.sort(key=lambda x: x[1], reverse=True)
    
    util.printcols([["TREE", "LOGL", "ERROR", "BASERATE", "COUNT"]] + 
                   mat[:80], spacing=4, out=DEBUG)
    debug()

    mat.sort(key=lambda x: x[2])
    util.printcols([["TREE", "LOGL", "ERROR", "BASERATE", "COUNT"]] + 
                   mat[:80], spacing=4, out=DEBUG)
    debug()


def numPossibleTrees(nleaves):
    n = 1.0
    
    for i in range(3, 2*nleaves-5+1, 2):
        n *= i
    
    return (2*nleaves - 3) * n


def log(x):
    """Safe logarithm function"""
    
    if x <= 0:
        return -util.INF
    else:
        return math.log(x)

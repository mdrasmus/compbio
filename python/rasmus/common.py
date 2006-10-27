"""
 common.py

 sets up a common environment for interactive scipting.  Used by my version
 of pyshell
"""


# rasmus direct imports
from rasmus.util import *
from rasmus.algorithms import *
from rasmus.fasta import *
from rasmus.matrix import *
from rasmus.stats import *
from rasmus.progress import *
from rasmus.vis.treevis import showTree

# rasmus modules
from rasmus import env, svg, tablelib
from rasmus.tablelib import Table, readTable

# bio tools
from rasmus import muscle, phylip, mrbayes, clustalw, genomeutil, genomeio, blast, alignlib
#from rasmus.phyloutil import viewTree



"""
# common compbio imports
try:
    from compbio.tools import pp
    from compbio import tools
except:
    pass
"""

# python imports
import os, sys, re
from math import *
import StringIO




def readTree(filename):
    """Read a newick tree"""
    tree = Tree()
    tree.readNewick(filename)
    return tree

class ListFiles:
    def __call__(self, args=""):
        return os.system("ls " + args)
    
    def __repr__(self):
        self()
        return ""
ls = ListFiles()


def cd(path = os.environ["HOME"]):
    os.chdir(path)


def strStream(text):
    return StringIO.StringIO(text)

# really quick pretty printing
pc = printcols
pa = alignlib.printAlign
pd = printDict


def pl(lines, out=sys.stdout):
    for line in lines:
        print >>out, line



# try to setup DATAPATH env
RASMUS_COMMON_DATAPATH_LOADED = False
if "DATAPATH" in os.environ:
    env.addEnvPaths("DATAPATH")
    RASMUS_COMMON_DATAPATH_LOADED = True
    
    print "rasmus.common: loaded data paths from DATAPATH"
    for path in env.datapaths:
        print path

    

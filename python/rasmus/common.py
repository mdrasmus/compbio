"""
 common.py

 sets up a common environment for interactive scipting.  Used by my version
 of pyshell
"""

# python libs
import os, sys, re
from math import *
import StringIO
from itertools import izip, imap


# rasmus direct imports
import rasmus
from rasmus.util import *
from rasmus.algorithms import *
from rasmus.matrix import *
from rasmus.stats import *
from rasmus.progress import *
from regionlib import *
#from rasmus.vis.treesvg import showTree

# rasmus modules
from rasmus import env, svg, tablelib, treelib
from rasmus.tablelib import Table, readTable
from rasmus.treelib import *

# bio tools
from rasmus.bio.fasta import *
from rasmus.bio import muscle, phylip, mrbayes, clustalw, genomeutil, blast, alignlib
from rasmus.bio import gff, genomeio








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
ls2 = ListFiles()


def cd2(path = os.environ["HOME"]):
    os.chdir(path)


def strStream(text):
    return StringIO.StringIO(text)

# really quick pretty printing
pc = printcols
pa = alignlib.printAlign
pd = printDict

histtab = tablelib.histTable


def pl(lines, out=sys.stdout):
    for line in lines:
        print >>out, line


def showTree(tree, **options):
    from rasmus.vis import treevis
    vis = treevis.TreeViewer(tree, **options)
    vis.show()
    
    return vis


# try to setup DATAPATH env
RASMUS_COMMON_DATAPATH_LOADED = False
if "DATAPATH" in os.environ:
    env.addEnvPaths("DATAPATH")
    RASMUS_COMMON_DATAPATH_LOADED = True
    
    print "rasmus.common: loaded data paths from DATAPATH"
    for path in env.datapaths:
        print path

    

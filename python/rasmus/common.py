"""
 common.py

 sets up a common environment for interactive scipting.  Used by my version
 of pyshell
"""


# rasmus imports
from util import *
from algorithms import *
from fasta import *
from matrix import *
from stats import *

# rasmus modules
import env, svg, alignlib

# bio tools
import muscle, phylip, mrbayes, clustalw, genomeutil, genomeio, blast

from phyloutil import viewTree

# common compbio imports
from compbio.tools import pp
from compbio import tools

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



# try to setup DATAPATH env
RASMUS_COMMON_DATAPATH_LOADED = False
if "DATAPATH" in os.environ:
    env.addEnvPaths("DATAPATH")
    RASMUS_COMMON_DATAPATH_LOADED = True
    
    print "rasmus.common: loaded data paths from DATAPATH"
    for path in env.datapaths:
        print path

    

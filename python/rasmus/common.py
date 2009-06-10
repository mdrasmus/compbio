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
from rasmus.matrixlib import *
from rasmus.stats import *
from rasmus.progress import *
from rasmus.regionlib import *


# rasmus modules
from rasmus import util, svg, tablelib, treelib
from rasmus.tablelib import Table, read_table, iter_table, histtab
from rasmus.tablelib import showtab, sqltab, sqlget, sqlput, sqlexe
from rasmus.treelib import *

# bio tools
from rasmus.bio.fasta import *
from rasmus.bio import muscle, phylip, mrbayes, blast, alignlib
from rasmus.bio import gff


readTable = read_table


def ipy():
    """start an ipython shell"""
    from IPython.Shell import IPShellEmbed
    ipshell = IPShellEmbed()
    ipshell()



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
pc = util.printcols
pa = alignlib.print_align
pd = util.print_dict



def pl(lines, out=sys.stdout):
    for line in lines:
        print >>out, line


def show_tree(tree, **options):
    from rasmus.vis import treevis
    vis = treevis.TreeViewer(tree, **options)
    vis.show()
    
    return vis

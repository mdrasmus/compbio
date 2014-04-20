"""
 common.py

 sets up a common environment for interactive scipting.  Used by my version
 of pyshell
"""

# python libs
from collections import defaultdict
from itertools import izip, imap, chain, islice, takewhile, dropwhile
from math import *
import os
import re
import StringIO
import sys
from time import sleep

# rasmus direct imports
import rasmus
from rasmus.gnuplot import *
from rasmus.matrixlib import *
from rasmus.plotting import *
from rasmus.rplotting import *
from rasmus.stats import *
from rasmus.treelib import *
from rasmus.util import *

# rasmus modules
from rasmus import util, svg, tablelib, treelib
from rasmus.tablelib import Table, read_table, iter_table, histtab
from rasmus.tablelib import showtab, sqlget, sqlput, sqlexe


# compbio tools
try:
    from compbio.fasta import *
    from compbio.regionlib import *
    from compbio import fasta, alignlib, gff, phylo
except ImportError:
    raise
    #pass


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


def cd2(path=os.environ["HOME"]):
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

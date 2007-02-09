#!/usr/bin/env python
#
# produce a subtree
#

from rasmus import treelib

import sys


# check args
if len(sys.argv) < 2:
    print >>sys.stderr, "usage: subtree.py <newick tree> <leaf1> <leaf2> ..."
    sys.exit(1)

tree = treelib.readTree(sys.argv[1])
keep = set(sys.argv[2:])


leaves = set(tree.leafNames())
delLeaves = leaves - keep


for leaf in delLeaves:
    tree.remove(tree.nodes[leaf])


treelib.removeExposedInternalNodes(tree)
treelib.removeSingleChildren(tree)


tree.write(writeData=lambda x: "")







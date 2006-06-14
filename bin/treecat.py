#!/usr/bin/env python

import sys
from rasmus import algorithms


tree = algorithms.Tree()
tree.makeRoot()

for f in sys.argv[1:]:
    print >>sys.stderr, f
    tree2 = algorithms.Tree()
    tree2.readNewick(f)
    tree.addTree(tree.root, tree2)

tree.write()


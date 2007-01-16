#!/usr/bin/python

import sys
import algorithms
import util

tree = algorithms.Tree()
treefile, labelfile = sys.argv[1:3]
tree.readParentTree(treefile, labelfile)

subtrees = algorithms.smallSubtrees(tree, 1000)

size = 0
for subtree in subtrees:
    part = subtree.leafNames()
    if len(part) > 100:
        size += len(part)
        for gene in part:
            print gene.split(":")[1],
        print

print >>sys.stderr, size

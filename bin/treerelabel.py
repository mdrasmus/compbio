#!/usr/bin/env python


from rasmus import util, algorithms, phylip
import sys

options = [
    ["t:", "tree=", "tree", "AUTO<top tree file>"],
    ["l:", "label=", "label", "AUTO<label file>"]
]

try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    sys.exit(1)


labels = util.readStrings(param["label"][-1])
tree = algorithms.Tree()
tree.readNewick(param["tree"][-1])

phylip.renameTreeWithNames(tree, labels)

tree.writeNewick(sys.stdout)


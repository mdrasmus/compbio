#!/usr/bin/env python


from rasmus import util, treelib
from rasmus.bio import phylip, fasta
import sys

options = [
    ["t:", "tree=", "tree", "<top tree file>"],
    ["l:", "label=", "label", "<label file>"]
]


param = util.parseOptions(sys.argv, options, quit=True)


if file(param["label"][-1]).read()[0] == ">":
    labels = fasta.readFasta(param["label"][-1]).keys()
else:
    labels = util.readStrings(param["label"][-1])
    
tree = treelib.Tree()
tree.readNewick(param["tree"][-1])

phylip.renameTreeWithNames(tree, labels)

tree.writeNewick(sys.stdout)


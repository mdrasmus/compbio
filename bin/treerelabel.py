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
    labels = fasta.read_fasta(param["label"][-1]).keys()
else:
    labels = util.read_strings(param["label"][-1])
    
tree = treelib.read_tree(param["tree"][-1])
phylip.rename_tree_with_names(tree, labels)
tree.write(sys.stdout)


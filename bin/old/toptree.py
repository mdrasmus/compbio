#!/usr/bin/env python

from rasmus import util, algorithms

import sys


options = [
    ["t:", "tree=", "tree", "AUTO<top tree file>"],
    ["p:", "part=", "part", "AUTO<part file>"]
]

try:
    param, rest = util.parseArgs(sys.argv, options)
except:
    sys.exit(1)

# read top tree
tree = algorithms.Tree()
tree.read_newick(param["tree"][-1])

if "part" in param:
    parts = util.read_delim(param["part"][-1])
    
    for i in xrange(len(parts)):
        num = str(i)
        if num in tree.nodes:
            tree2 = algorithms.Tree()
            tree2.make_root()
            for gene in parts[i]:
                tree2.add_child(tree2.root, algorithms.TreeNode(gene))
            tree.replace_tree(tree.nodes[num], tree2)
    
else:
    # read subtrees
    for f in rest:
        tree2 = algorithms.Tree()
        tree2.readNewick(f)

        # find number of tree in name
        num = filter(lambda x: x.isdigit(), f.split("/")[-1])

        if num in tree.nodes:
            tree.replaceTree(tree.nodes[num], tree2)

tree.writeNewick(sys.stdout)


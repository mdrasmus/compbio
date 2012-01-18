

import unittest
from unittest import TestCase

from rasmus.common import *
from rasmus import stats
from rasmus.testing import *

from compbio import phylo


#=============================================================================
# reconciliations

    
class Recon (TestCase):

    def test_enum_recon(self):

        tree = parse_newick("((a,b),(c,d))")
        stree = parse_newick("((a,b),(c,d))")
        gene2species = lambda x: x

        for recon, events in phylo.enum_recon(tree, stree, depth=None,
                                              gene2species=gene2species):
            print "recon"
            print_dict(recon)
            print
            
            print "events"
            print_dict(events)
            print


#=============================================================================
# tree search

class Search (TestCase):

    def test(self):

        tree = parse_newick("((a,b),((c,d),(e,f)))")
        treelib.draw_tree_names(tree, minlen=5)

        a, b = phylo.propose_random_spr(tree)
        print a.name, b.name

        phylo.perform_spr(tree, a, b)
        treelib.draw_tree_names(tree, minlen=5)


        for i in xrange(100):
            top1 = phylo.hash_tree(tree)
            print "forward"
            s = phylo.TreeSearchSpr(tree)
            s.next()
            treelib.draw_tree_names(tree, minlen=5)
            top2 = phylo.hash_tree(tree)

            assert top1 != top2

            print "revert"
            s.revert()
            treelib.draw_tree_names(tree, minlen=5)

            assert phylo.hash_tree(tree) == top1


#=============================================================================


if __name__ == "__main__":
    test_main()

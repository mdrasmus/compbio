
from unittest import TestCase

from rasmus import treelib
from rasmus.treelib import parse_newick

from compbio import phylo


class Recon (TestCase):
    """Gene-tree species-tree reconciliation (recon)"""

    def test_enum_recon(self):
        """Recon enumeration should always produce valid recons"""

        expected_recons = [
            {'a': 'a', 1: 1, 'c': 'c', 3: 3, 'd': 'd', 2: 2, 'b': 'b'},
            {'a': 'a', 1: 1, 'c': 'c', 3: 3, 'd': 'd', 2: 2, 'b': 'b'},
            {'a': 'a', 1: 1, 'c': 'c', 3: 3, 'd': 'd', 2: 2, 'b': 'b'},
            {'a': 'a', 1: 1, 'c': 'c', 3: 3, 'd': 'd', 2: 1, 'b': 'b'},
            {'a': 'a', 1: 1, 'c': 'c', 3: 3, 'd': 'd', 2: 1, 'b': 'b'},
            {'a': 'a', 1: 1, 'c': 'c', 3: 1, 'd': 'd', 2: 1, 'b': 'b'},
            {'a': 'a', 1: 1, 'c': 'c', 3: 3, 'd': 'd', 2: 2, 'b': 'b'},
            {'a': 'a', 1: 1, 'c': 'c', 3: 1, 'd': 'd', 2: 2, 'b': 'b'},
            {'a': 'a', 1: 1, 'c': 'c', 3: 3, 'd': 'd', 2: 2, 'b': 'b'},
            {'a': 'a', 1: 1, 'c': 'c', 3: 1, 'd': 'd', 2: 2, 'b': 'b'},
            {'a': 'a', 1: 1, 'c': 'c', 3: 3, 'd': 'd', 2: 2, 'b': 'b'},
            {'a': 'a', 1: 1, 'c': 'c', 3: 3, 'd': 'd', 2: 1, 'b': 'b'},
            {'a': 'a', 1: 1, 'c': 'c', 3: 3, 'd': 'd', 2: 1, 'b': 'b'},
            {'a': 'a', 1: 1, 'c': 'c', 3: 1, 'd': 'd', 2: 1, 'b': 'b'},
            {'a': 'a', 1: 1, 'c': 'c', 3: 3, 'd': 'd', 2: 2, 'b': 'b'},
            {'a': 'a', 1: 1, 'c': 'c', 3: 1, 'd': 'd', 2: 2, 'b': 'b'},
            {'a': 'a', 1: 1, 'c': 'c', 3: 3, 'd': 'd', 2: 2, 'b': 'b'},
            {'a': 'a', 1: 1, 'c': 'c', 3: 1, 'd': 'd', 2: 2, 'b': 'b'},
        ]

        expected_events = [
            {'a': 'gene', 1: 'spec', 'c': 'gene', 3: 'spec',
             'd': 'gene', 2: 'spec', 'b': 'gene'},
            {'a': 'gene', 1: 'dup', 'c': 'gene', 3: 'spec',
             'd': 'gene', 2: 'spec', 'b': 'gene'},
            {'a': 'gene', 1: 'dup', 'c': 'gene', 3: 'spec',
             'd': 'gene', 2: 'dup', 'b': 'gene'},
            {'a': 'gene', 1: 'dup', 'c': 'gene', 3: 'spec',
             'd': 'gene', 2: 'dup', 'b': 'gene'},
            {'a': 'gene', 1: 'dup', 'c': 'gene', 3: 'dup',
             'd': 'gene', 2: 'dup', 'b': 'gene'},
            {'a': 'gene', 1: 'dup', 'c': 'gene', 3: 'dup',
             'd': 'gene', 2: 'dup', 'b': 'gene'},
            {'a': 'gene', 1: 'dup', 'c': 'gene', 3: 'dup',
             'd': 'gene', 2: 'dup', 'b': 'gene'},
            {'a': 'gene', 1: 'dup', 'c': 'gene', 3: 'dup',
             'd': 'gene', 2: 'dup', 'b': 'gene'},
            {'a': 'gene', 1: 'dup', 'c': 'gene', 3: 'dup',
             'd': 'gene', 2: 'spec', 'b': 'gene'},
            {'a': 'gene', 1: 'dup', 'c': 'gene', 3: 'dup',
             'd': 'gene', 2: 'spec', 'b': 'gene'},
            {'a': 'gene', 1: 'spec', 'c': 'gene', 3: 'spec',
             'd': 'gene', 2: 'dup', 'b': 'gene'},
            {'a': 'gene', 1: 'spec', 'c': 'gene', 3: 'spec',
             'd': 'gene', 2: 'dup', 'b': 'gene'},
            {'a': 'gene', 1: 'spec', 'c': 'gene', 3: 'dup',
             'd': 'gene', 2: 'dup', 'b': 'gene'},
            {'a': 'gene', 1: 'spec', 'c': 'gene', 3: 'dup',
             'd': 'gene', 2: 'dup', 'b': 'gene'},
            {'a': 'gene', 1: 'spec', 'c': 'gene', 3: 'dup',
             'd': 'gene', 2: 'dup', 'b': 'gene'},
            {'a': 'gene', 1: 'spec', 'c': 'gene', 3: 'dup',
             'd': 'gene', 2: 'dup', 'b': 'gene'},
            {'a': 'gene', 1: 'spec', 'c': 'gene', 3: 'dup',
             'd': 'gene', 2: 'spec', 'b': 'gene'},
            {'a': 'gene', 1: 'spec', 'c': 'gene', 3: 'dup',
             'd': 'gene', 2: 'spec', 'b': 'gene'},
        ]

        tree = parse_newick("((a,b),(c,d))")
        stree = parse_newick("((a,b),(c,d))")
        gene2species = lambda x: x

        for i, (recon, events) in enumerate(phylo.enum_recon(
                tree, stree, depth=None, gene2species=gene2species)):
            phylo.assert_recon(tree, stree, recon)
            recon_names = dict((node.name, snode.name)
                               for node, snode in recon.items())
            event_names = dict((node.name, event)
                               for node, event in events.items())

            self.assertEqual(recon_names, expected_recons[i])
            self.assertEqual(event_names, expected_events[i])


class Search (TestCase):
    """Tree search"""

    def test(self):
        """Test a tree search"""

        tree = parse_newick("((a,b),((c,d),(e,f)))")

        a, b = phylo.propose_random_spr(tree)
        phylo.perform_spr(tree, a, b)
        treelib.assert_tree(tree)

        for i in xrange(100):
            top1 = phylo.hash_tree(tree)

            s = phylo.TreeSearchSpr(tree)
            s.next()
            top2 = phylo.hash_tree(tree)

            self.assertNotEqual(top1, top2)

            s.revert()
            self.assertEqual(phylo.hash_tree(tree), top1)


class Splits (TestCase):
    """Tree bi-partitions (splits)"""

    def test(self):

        tree1 = parse_newick("((a,b),c)")
        tree2 = parse_newick("(c,(a,b))")
        assert (phylo.find_splits(tree1, rooted=True) ==
                phylo.find_splits(tree2, rooted=True))

        tree1 = parse_newick("((a,b),(c,d))")
        tree2 = parse_newick("(((c,d),a),b)")
        assert (phylo.find_splits(tree1) ==
                phylo.find_splits(tree2))

        assert phylo.robinson_foulds_error(tree1, tree2) == 0.0

        tree1 = parse_newick("(((a,b),(c,d)),(e,f))")
        tree2 = parse_newick("(((a,c),(b,d)),(e,f))")
        self.assertAlmostEqual(phylo.robinson_foulds_error(tree1, tree2), 2/3.)

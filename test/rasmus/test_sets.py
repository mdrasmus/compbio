
from unittest import TestCase

from rasmus.sets import UnionFind


class SetTests(TestCase):
    def test_union_find(self):
        set1 = UnionFind()
        set2 = UnionFind()
        set3 = UnionFind()

        set1.add(1)
        set1.add(2)

        self.assertEqual(len(set1), 2)
        self.assertTrue(1 in set1)
        self.assertFalse(set1.has(-1))

        set2.add(3)
        set2.add(4)
        set2.add(5)
        self.assertEqual(len(set2), 3)

        set3.add(5)
        set3.add(6)
        set3.add(7)
        self.assertEqual(len(set3), 3)

        self.assertFalse(set1.same(set2))
        set1.union(set2)
        self.assertTrue(set1.same(set2))

        set1.union(set3)
        self.assertTrue(
            set1.members(),
            set([1, 2, 3, 4, 5, 6, 7]))
        self.assertTrue(len(set1), len(set2))

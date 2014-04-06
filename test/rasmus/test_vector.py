import unittest

from rasmus import vector as v


class Test (unittest.TestCase):

    def test_list(self):
        a = [1.0, 2.0, 3.0]
        b = [4.0, 5.0, 6.0]

        self.assertEqual(v.vadd(a, b),
                         [5.0, 7.0, 9.0])
        self.assertEqual(v.vsub(a, b),
                         [-3.0, -3.0, -3.0])
        self.assertEqual(v.vmul(a, b),
                         [4.0, 10.0, 18.0])
        self.assertEqual(v.vdiv(a, b),
                         [0.25, 0.4, 0.5])
        self.assertAlmostEqual(v.vmag(a), 3.74165738677)
        self.assertAlmostEqual(v.vdist(a, b), 5.19615242271)

    def _test_dict(self):
        a = {'x': 1.0, 'y': 2.0, 'z': 3.0}
        b = {'x': 4.0, 'y': 5.0, 'z': 6.0}

        self.assertEqual(v.vadd(a, b),
                         dict(zip('xyz', [5.0, 7.0, 9.0])))
        self.assertEqual(v.vsub(a, b),
                         dict(zip('xyz', [-3.0, -3.0, -3.0])))
        self.assertEqual(v.vmul(a, b),
                         dict(zip('xyz', [4.0, 10.0, 18.0])))
        self.assertEqual(v.vdiv(a, b),
                         dict(zip('xyz', [0.25, 0.4, 0.5])))
        self.assertAlmostEqual(v.vmag(a), 3.74165738677)
        self.assertAlmostEqual(v.vdist(a, b), 5.19615242271)

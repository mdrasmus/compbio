
import sys
import unittest

from rasmus import util


class Test (unittest.TestCase):

    def test_open_stream1(self):
        """open_stream shouldn't close existing stream"""

        infile = util.open_stream(sys.stdin)

        # ensure attribute access
        infile.read

        # make sure file doesn't close
        infile.close()
        assert not sys.stdin.closed

    def test_open_stream2(self):
        """open_stream should close file"""

        # make sure regular files close
        infile = util.open_stream(__file__)
        infile.close()
        assert infile.closed

    def test_buckets(self):
        """Test bucket_bin"""

        bin = util.bucket_bin(50, 20, 0, 5)
        assert bin == 10

    def test_pretty_int(self):

        self.assertEqual(util.int2pretty(0), '0')
        self.assertEqual(util.int2pretty(1), '1')
        self.assertEqual(util.int2pretty(10), '10')
        self.assertEqual(util.int2pretty(100), '100')
        self.assertEqual(util.int2pretty(1000), '1,000')
        self.assertEqual(util.int2pretty(10000), '10,000')
        self.assertEqual(util.int2pretty(100000), '100,000')
        self.assertEqual(util.int2pretty(1000000), '1,000,000')

        self.assertEqual(util.int2pretty(-1), '-1')
        self.assertEqual(util.int2pretty(-10), '-10')
        self.assertEqual(util.int2pretty(-100), '-100')
        self.assertEqual(util.int2pretty(-1000), '-1,000')
        self.assertEqual(util.int2pretty(-10000), '-10,000')
        self.assertEqual(util.int2pretty(-100000), '-100,000')
        self.assertEqual(util.int2pretty(-1000000), '-1,000,000')

        self.assertEqual(util.pretty2int('0'), 0)
        self.assertEqual(util.pretty2int('1'), 1)
        self.assertEqual(util.pretty2int('10'), 10)
        self.assertEqual(util.pretty2int('100'), 100)
        self.assertEqual(util.pretty2int('1,000'), 1000)
        self.assertEqual(util.pretty2int('10,000'), 10000)
        self.assertEqual(util.pretty2int('100,000'), 100000)
        self.assertEqual(util.pretty2int('1,000,000'), 1000000)

        self.assertEqual(util.pretty2int('-1'), -1)
        self.assertEqual(util.pretty2int('-10'), -10)
        self.assertEqual(util.pretty2int('-100'), -100)
        self.assertEqual(util.pretty2int('-1,000'), -1000)
        self.assertEqual(util.pretty2int('-10,000'), -10000)
        self.assertEqual(util.pretty2int('-100,000'), -100000)
        self.assertEqual(util.pretty2int('-1,000,000'), -1000000)

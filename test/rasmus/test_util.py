
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

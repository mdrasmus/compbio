
import unittest

from rasmus.common import *
from rasmus.testing import *

from rasmus import tablelib
from StringIO import StringIO


#=============================================================================


    
class Test (unittest.TestCase):

    def test_open_stream1(self):

        infile = open_stream(sys.stdin)

        # ensure attribute access
        infile.read

        # make sure file doesn't close
        infile.close()
        assert not sys.stdin.closed


    def test_open_stream2(self):
        
        # make sure regular files close
        infile = open_stream(__file__)
        infile.close()
        assert infile.closed
        

    def test_buckets(self):

        bin = bucket_bin(50, 20, 0, 5)
        assert bin == 10


#=============================================================================
if __name__ == "__main__":   
    unittest.main()




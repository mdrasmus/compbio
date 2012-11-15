
import unittest
import timeit

from rasmus.common import *
from rasmus.testing import *


from rasmus import treelib
from StringIO import StringIO


#=============================================================================

fungi = """(((((((scer:7.061760[&&NHX:a=b],spar:7.061760):4.999680,smik:12.061440):5.970600,sbay:18.032040):52.682400,cgla:70.714260):7.220700,scas:77.934960):23.181480,((agos:78.553260,klac:78.553260):10.434960,kwal:88.988220):12.128400):78.883560,(((calb:41.275620,ctro:41.275980):29.632860,(cpar:52.323120,lelo:52.323120):18.585720):31.149540,((cgui:75.615840,dhan:75.615840):14.006880,clus:89.622720):12.435660):77.941620);"""
    
class Test (unittest.TestCase):

    def test_read_tree_speed(self):
        t = timeit.timeit('treelib.parse_newick(fungi)',
                          setup='from rasmus import treelib; ' +
                          'from __main__ import fungi',
                          number=1000)
        print "time", t

        #t = timeit.timeit('treelib.parse_newick2(fungi)',
        #                  setup='from rasmus import treelib;' +
        #                  'from __main__ import fungi',
        #                  number=1000)
        #print "time", t


    def test_tokenize_newick(self):
        fungi_infile = StringIO(fungi)
        print list(treelib.tokenize_newick(fungi_infile))
        print list(treelib.tokenize_newick(fungi))


    def test_read_tree(self):
        tree = treelib.read_tree(StringIO(fungi))

        print tree.nodes
        print tree["scer"].data
        
        newick = tree.get_one_line_newick(writeData=treelib.write_nhx_data)
        print newick
        assert newick == fungi


    def test_iter_trees(self):
        for tree in treelib.iter_trees(StringIO(fungi + fungi + fungi)):
            print tree.get_one_line_newick()
    
    

#=============================================================================
if __name__ == "__main__":   
    test_main()




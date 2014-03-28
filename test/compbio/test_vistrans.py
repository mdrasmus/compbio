
from StringIO import StringIO
import unittest

from compbio import phylo
from compbio.vis import transsvg
from rasmus.testing import make_clean_dir
from rasmus import treelib


stree_newick = """((((sde1:1,sde2:1):2,sdd:3):1,(spy1:2,spy2:2):2):1,see:5);"""


smapfile = StringIO("""\
SDE1*	sde1
SDEG*	sde2
SDD*	sdd
SpyM3_*	spy1
MGAS10750_*	spy2
SEQ_*	see
sde1_*	sde1
sde2_*	sde2
sdd_*	sdd
spy1_*	spy1
spy2_*	spy2
see_*	see""")


#=============================================================================

treefile1 = StringIO("""(
 (
  SEQ_1914:0.000000,
  SEQ_1956:0.000000
 )1.000000:0.112215,
 (
  (
   (
    (
     SDE12394_01485:0.000000,
     SDE12394_01705:0.000000
    )0.800000:0.000000,
    SDEG_0318:0.001920
   )1.000000:0.027307,
   (
    SDD27957_01490:0.000000,
    SDD27957_01740:0.000000
   )0.700000:0.000100
  )1.000000:0.129840,
  (
   (
    SpyM3_0220:0.000000,
    SpyM3_1640:0.000000
   )0.790000:0.000100,
   (
    MGAS10750_Spy1673:0.000000,
    MGAS10750_Spy0253:0.000000
   )1.000000:0.003725
  )1.000000:0.199601
 )1.000000:0.112215
);
""")

breconfile1 = StringIO("""\
1	1	spec
7	sdd	dup
SDE12394_01705	sde1	gene
2	see	dup
4	3	spec
3	2	spec
SDD27957_01490	sdd	gene
MGAS10750_Spy0253	spy2	gene
5	4	spec
8	5	spec
10	spy2	dup
9	spy1	dup
SDE12394_01485	sde1	gene
6	sde1	dup
SEQ_1914	see	gene
MGAS10750_Spy1673	spy2	gene
SpyM3_0220	spy1	gene
SEQ_1956	see	gene
SpyM3_1640	spy1	gene
SDD27957_01740	sdd	gene
SDEG_0318	sde2	gene
""")


treefile2 = StringIO("""(
 SEQ_1075:0.058336,
 (
  (
   (
    (
     SDD27957_10755:0.016483,
     SDEG_2059:0.009287
    )0.980000:0.010526,
    (
     SpyM3_1776:0.005808,
     MGAS10750_Spy1866:0.005792
    )0.990000:0.010198
   )0.710000:0.005186,
   SDE12394_10390:0.002079
  )1.000000:0.305782,
  (
   (
    (
     SDE12394_05780:0.001163,
     SDEG_1033:0.000000
    )0.800000:0.004425,
    SDD27957_05500:0.006104
   )0.520000:0.048296,
   (
    SpyM3_0853:0.003142,
    MGAS10750_Spy1076:0.002729
   )1.000000:0.069646
  )1.000000:0.074916
 )0:0.058336
);
""")

breconfile2 = StringIO("""\
MGAS10750_Spy1076	spy2	gene
8	3	spec
6	5	spec
5	sdd	trans
4	sdd	trans
SEQ_1075	see	gene
7	5	trans
SDD27957_10755	sdd	gene
SDE12394_05780	sde1	gene
1	1	spec
10	5	spec
SpyM3_1776	spy1	gene
SDEG_1033	sde2	gene
SDE12394_10390	4	specloss	sde1	gene
2	2	spec
SDD27957_05500	sdd	gene
3	3	spec
SDEG_2059	sde2	gene
SpyM3_0853	spy1	gene
9	4	spec
MGAS10750_Spy1866	spy2	gene
""")


treefile3 = StringIO("""(
 SEQ_0864:0.094945,
 (
  SDEG_0677:0.117023,
  (
   SDD27957_05160:0.022882,
   (
    SpyM3_1656:0.005147,
    MGAS10750_Spy1732:0.004911
   )1.000000:0.011928
  )1.000000:0.237151
 )0:0.094945
);
""")

breconfile3 = StringIO("""\
SEQ_0864	see	gene
SDEG_0677	3	specloss	4	specloss	sde2	gene
MGAS10750_Spy1732	spy2	gene
4	5	spec
SpyM3_1656	spy1	gene
3	5	trans
2	2	spec
1	1	spec
SDD27957_05160	sdd	gene
""")


class Vis (unittest.TestCase):

    def test1(self):
        outdir = 'test/tmp/test_vistrans/Vis_test1/'
        make_clean_dir(outdir)

        stree = treelib.parse_newick(stree_newick)
        tree = treelib.read_tree(treefile1)
        brecon = phylo.read_brecon(breconfile1, tree, stree)

        transsvg.draw_tree(tree, brecon, stree, filename=outdir + "tree.svg")

    def test2(self):
        outdir = 'test/tmp/test_vistrans/Vis_test2/'
        make_clean_dir(outdir)

        stree = treelib.parse_newick(stree_newick)
        tree = treelib.read_tree(treefile2)
        brecon = phylo.read_brecon(breconfile2, tree, stree)

        transsvg.draw_tree(tree, brecon, stree, filename=outdir + "tree.svg")

    def test3(self):
        outdir = 'test/tmp/test_vistrans/Vis_test3/'
        make_clean_dir(outdir)

        stree = treelib.parse_newick(stree_newick)
        tree = treelib.read_tree(treefile3)
        brecon = phylo.read_brecon(breconfile3, tree, stree)

        phylo.add_implied_spec_nodes_brecon(tree, brecon)
        phylo.write_brecon(open(outdir + 'brecon', 'w'), brecon)

        transsvg.draw_tree(tree, brecon, stree, filename=outdir + "tree.svg")

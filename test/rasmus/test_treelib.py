
from StringIO import StringIO
import timeit
import unittest

from rasmus import treelib
from rasmus.treelib import read_tree
from rasmus.treelib import reorder_tree


fungi = """(((((((scer:7.061760,spar:7.061760):4.999680,smik:12.061440):5.970600,sbay:18.032040):52.682400,cgla:70.714260):7.220700,scas:77.934960):23.181480,((agos:78.553260,klac:78.553260):10.434960,kwal:88.988220):12.128400):78.883560,(((calb:41.275620,ctro:41.275980):29.632860,(cpar:52.323120,lelo:52.323120):18.585720):31.149540,((cgui:75.615840,dhan:75.615840):14.006880,clus:89.622720):12.435660):77.941620);"""  # nopep8

fungi2 = """(((((((scer:7.061760[&&NHX:a=b],spar:7.061760):4.999680,smik:12.061440):5.970600,sbay:18.032040):52.682400,cgla:70.714260):7.220700,scas:77.934960):23.181480,((agos:78.553260,klac:78.553260):10.434960,kwal:88.988220):12.128400):78.883560,(((calb:41.275620,ctro:41.275980):29.632860,(cpar:52.323120,lelo:52.323120):18.585720):31.149540,((cgui:75.615840,dhan:75.615840):14.006880,clus:89.622720):12.435660)xx:77.941620);"""  # nopep8


class TreeTest (unittest.TestCase):

    def test_tokenize_newick(self):
        """Test newick tokenization"""
        text = "((A:10,B:1.2)xx:22,(C,D))aaa;"
        tokens = ['(', '(', 'A', ':', '10', ',', 'B', ':', '1.2', ')',
                  'xx', ':', '22', ',', '(', 'C', ',', 'D', ')', ')',
                  'aaa', ';']

        tokens2 = list(treelib.tokenize_newick(text))
        self.assertEquals(tokens, tokens2)

        text = "((A:10,B:1.2[ aaa bb ])xx:22,(C,D)) aaa [xx] ;"
        tokens = ['(', '(', 'A', ':', '10', ',', 'B', ':', '1.2',
                  '[ aaa bb ]', ')',
                  'xx', ':', '22', ',', '(', 'C', ',', 'D', ')', ')',
                  'aaa', '[xx]', ';']

        tokens2 = list(treelib.tokenize_newick(text))
        self.assertEquals(tokens, tokens2)

        text = "((A:10,B:1.2[ aaa[  bb[ ])xx:22,(C,D))aaa;"
        tokens = ['(', '(', 'A', ':', '10', ',', 'B', ':', '1.2',
                  '[ aaa[  bb[ ]', ')',
                  'xx', ':', '22', ',', '(', 'C', ',', 'D', ')', ')',
                  'aaa', ';']

        tokens2 = list(treelib.tokenize_newick(text))
        self.assertEquals(tokens, tokens2)

    def test_read_tree(self):
        """Test reading tree structure."""
        tree = treelib.read_tree(StringIO(fungi2))
        ptree = dict((node.name, node.parent.name if node.parent else None)
                     for node in tree)
        ptree_expected = {
            1: None, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 6, 8: 2, 9: 8, 10: 'xx',
            11: 10, 12: 10, 'sbay': 5, 14: 13, 'xx': 1, 'scer': 7, 'ctro': 11,
            'scas': 3, 'agos': 9, 'kwal': 8, 'dhan': 14, 'smik': 6, 'cgla': 4,
            'spar': 7, 'calb': 11, 'lelo': 12, 'cpar': 12, 13: 'xx', 'klac': 9,
            'clus': 13, 'cgui': 14}
        self.assertEqual(ptree, ptree_expected)

        newick = tree.get_one_line_newick(writeData=treelib.write_nhx_data)
        self.assertEqual(newick, fungi2)

    def test_read_tree_speed(self):
        """Test speed of tree reading."""
        t = timeit.timeit('treelib.parse_newick(fungi2)',
                          setup='from rasmus import treelib; ' +
                          'fungi2 = %s' % repr(fungi2),
                          number=1000)
        print "time", t

    def test_iter_trees(self):
        """Test iterative parsing of many trees."""
        trees = list(treelib.iter_trees(StringIO(fungi + fungi + fungi)))
        self.assertEqual(len(trees), 3)

    def test_nhx(self):
        """Test parsing of NHX comments."""

        text = """(((ADH2:0.1[&&NHX:S=human:E=1.1.1.1], ADH1:0.11[&&NHX:S=human:E=1.1.1.1]):0.05[&&NHX:S=Primates:E=1.1.1.1:D=Y:B=100], ADHY:0.1[&&NHX:S=nematode:E=1.1.1.1],ADHX:0.12[&&NHX:S=insect:E=1.1.1.1]):0.1[&&NHX:S=Metazoa:E=1.1.1.1:D=N], (ADH4:0.09[&&NHX:S=yeast:E=1.1.1.1],ADH3:0.13[&&NHX:S=yeast:E=1.1.1.1], ADH2:0.12[&&NHX:S=yeast:E=1.1.1.1],ADH1:0.11[&&NHX:S=yeast:E=1.1.1.1]):0.1 [&&NHX:S=Fungi])[&&NHX:E=1.1.1.1:D=N];"""  # nopep8
        tree = read_tree(StringIO(text))

        data = {'ADH3': {'S': 'yeast', 'E': '1.1.1.1'},
                1: {'E': '1.1.1.1', 'D': 'N'},
                2: {'S': 'Metazoa', 'E': '1.1.1.1', 'D': 'N'},
                3: {'S': 'Primates', 'B': '100', 'E': '1.1.1.1', 'D': 'Y'},
                4: {'S': 'Fungi'},
                'ADH2': {'S': 'human', 'E': '1.1.1.1'},
                'ADHY': {'S': 'nematode', 'E': '1.1.1.1'},
                'ADHX': {'S': 'insect', 'E': '1.1.1.1'},
                'ADH1': {'S': 'human', 'E': '1.1.1.1'},
                'ADH1_1': {'S': 'yeast', 'E': '1.1.1.1'},
                'ADH4': {'S': 'yeast', 'E': '1.1.1.1'},
                'ADH2_1': {'S': 'yeast', 'E': '1.1.1.1'}}

        data2 = dict((node.name, node.data) for node in tree)

        for key, val in data.items():
            self.assertEqual(data2[key], val)

    def test_nhx_big(self):
        """Test parsing of big NHX comments."""
        text = """(CFTR_GASAC:0.028272[&&NHX:S=GASAC:O=ENSGACT00000011967.1:T=69293:G=ENSGACG00000009039],((((((((((((((((((CFTR_HUMAN:0.002013[&&NHX:S=HUMAN:O=ENST00000003084.5:T=9606:G=ENSG00000001626],CFTR_PANTR:0.001342[&&NHX:S=PANTR:O=ENSPTRT00000036339.2:T=9598:G=ENSPTRG00000019619]):0.001545,CFTR_PONPY:0.006514[&&NHX:S=PONPY:O=ENSPPYT00000020909.1:T=9600:G=ENSPPYG00000017940]):0.003539,CFTR_MACMU:0.008416[&&NHX:S=MACMU:O=ENSMMUT00000015762.2:T=9544:G=ENSMMUG00000011269]):0.022751,CFTR_TUPGB:0.110613[&&NHX:S=TUPGB:O=ENSTBET00000011046.1:T=37347:G=ENSTBEG00000010974]):0.006474,((CFTR_OTOGA:0.035577[&&NHX:S=OTOGA:O=ENSOGAT00000001759.1:T=30611:G=ENSOGAG00000001756],CFTR_MICMU:0.026588[&&NHX:S=MICMU:O=ENSMICT00000005779.1:T=30608:G=ENSMICG00000005761]):0.010514,CFTR_MYOLU:0.06919[&&NHX:S=MYOLU:O=ENSMLUT00000012267.1:T=59463:G=ENSMLUG00000012244]):0.00395):0.001879,(CFTR_ECHTE:0.065629[&&NHX:S=ECHTE:O=ENSETET00000000538.1:T=9371:G=ENSETEG00000000537],CFTR_LOXAF:0.050347[&&NHX:S=LOXAF:O=ENSLAFT00000005758.1:T=9785:G=ENSLAFG00000005753]):0.016592):0.002471,((CFTR_SORAR:0.056771[&&NHX:S=SORAR:O=ENSSART00000012124.1:T=42254:G=ENSSARG00000012121],CFTR_ERIEU:0.043527[&&NHX:S=ERIEU:O=ENSEEUT00000006570.1:T=9365:G=ENSEEUG00000006484]):0.015585,CFTR_DASNO:0.047157[&&NHX:S=DASNO:O=ENSDNOT00000016544.1:T=9361:G=ENSDNOG00000016541]):0.00431):0.005677,(CFTR_F2_HORSE:0.016035[&&NHX:S=HORSE:O=ENSECAT00000010738.1:T=9796:G=ENSECAG00000009139],((CFTR_CANFA:0.047251[&&NHX:S=CANFA:O=ENSCAFT00000005518.2:T=9615:G=ENSCAFG00000003429],Q9N1D7_FELCA:0.025264[&&NHX:S=FELCA:O=ENSFCAT00000014959.2:T=9685:G=ENSFCAG00000014955]):0.022297,CFTR_BOVIN:0.062409[&&NHX:S=BOVIN:O=ENSBTAT00000053450.1:T=9913:G=ENSBTAG00000006589]):0.00767):0.004191):0.006209,(CFTR_F2_CAVPO:0.136979[&&NHX:S=CAVPO:O=ENSCPOT00000012891.1:T=10141:G=ENSCPOG00000012767],CFTR_SPETR:0.026944[&&NHX:S=SPETR:O=ENSSTOT00000005733.1:T=43179:G=ENSSTOG00000005707]):0.009628):0.007329,(Q29399_RABIT:0.027324[&&NHX:S=RABIT:O=ENSOCUT00000010738.1:T=9986:G=ENSOCUG00000010733],CFTR_OCHPR:0.050953[&&NHX:S=OCHPR:O=ENSOPRT00000014760.1:T=9978:G=ENSOPRG00000014721]):0.017472):0.011797,(Cftr_MOUSE:0.035769[&&NHX:S=MOUSE:O=ENSMUST00000045706.4:T=10090:G=ENSMUSG00000041301],Cftr_RAT:0.049345[&&NHX:S=RAT:O=ENSRNOT00000010981.4:T=10116:G=ENSRNOG00000008284]):0.158692):0.033423,Q2QL94_MONDO:0.08197[&&NHX:S=MONDO:O=ENSMODT00000020031.2:T=13616:G=ENSMODG00000015771]):0.026265,CFTR_ORNAN:0.094961[&&NHX:S=ORNAN:O=ENSOANT00000013974.1:T=9258:G=ENSOANG00000008767]):0.03792,A0M8U4_CHICK:0.119618[&&NHX:S=CHICK:O=ENSGALT00000015182.3:T=9031:G=ENSGALG00000009324]):0.033083,CFTR_XENTR:0.130489[&&NHX:S=XENTR:O=ENSXETT00000047145.1:T=8364:G=ENSXETG00000021796]):0.352249,si_dkey-270i2_F3_BRARE:0.203525[&&NHX:S=BRARE:O=ENSDART00000100729.1:T=7955:G=ENSDARG00000041107]):0.063334,CFTR_ORYLA:0.123603[&&NHX:S=ORYLA:O=ENSORLT00000024332.1:T=8090:G=ENSORLG00000019555]):0.034773,CFTR_TETNG:0.049086[&&NHX:S=TETNG:O=ENSTNIT00000019381.1:T=99883:G=ENSTNIG00000016063]):0.028272)[&&NHX:Loglk=-24078.827174:RatioCons=0.000000;:LoglkSpec=0.000000];"""  # nopep8
        tree = read_tree(StringIO(text))
        expected = {
            29: {},
            30: {},
            'CFTR_MACMU': {'O': 'ENSMMUT00000015762.2',
                           'S': 'MACMU',
                           'T': '9544',
                           'G': 'ENSMMUG00000011269'},
        }

        for name, data in expected.items():
            self.assertEqual(tree[name].data, data)

    def test_write_tree(self):
        """Test tree writing
           Test root data writing
        """

        newick = '''(
 (
  a:1.000000,
  b:2.000000
 )x:3.000000,
 (
  c:4.000000,
  d:5.000000
 )y:6.000000
)rra:0.000000;
'''
        infile = StringIO(newick)
        tree = read_tree(infile)

        out = StringIO()
        tree.write(out, rootData=True)
        self.assertEqual(newick, out.getvalue())

    def test_tree_namefunc(self):
        """Test reading/writing tree with namefunc."""

        count = [0]

        def namefunc(name):
            count[0] += 1
            return 'name%d' % count[0]

        tree = treelib.read_tree(StringIO(fungi2), namefunc=namefunc)
        newick = tree.get_one_line_newick()
        expected_newick = '(((((((name1:7.061760,name2:7.061760):4.999680,name3:12.061440):5.970600,name4:18.032040):52.682400,name5:70.714260):7.220700,name6:77.934960):23.181480,((name7:78.553260,name8:78.553260):10.434960,name9:88.988220):12.128400):78.883560,(((name10:41.275620,name11:41.275980):29.632860,(name12:52.323120,name13:52.323120):18.585720):31.149540,((name14:75.615840,name15:75.615840):14.006880,name16:89.622720):12.435660)xx:77.941620);'  # nopep8

        self.assertEqual(newick, expected_newick)

        def namefunc2(name):
            return 'prefix_' + name

        newick2 = tree.get_one_line_newick(namefunc=namefunc2)
        expected_newick2 = '(((((((prefix_name1:7.061760,prefix_name2:7.061760):4.999680,prefix_name3:12.061440):5.970600,prefix_name4:18.032040):52.682400,prefix_name5:70.714260):7.220700,prefix_name6:77.934960):23.181480,((prefix_name7:78.553260,prefix_name8:78.553260):10.434960,prefix_name9:88.988220):12.128400):78.883560,(((prefix_name10:41.275620,prefix_name11:41.275980):29.632860,(prefix_name12:52.323120,prefix_name13:52.323120):18.585720):31.149540,((prefix_name14:75.615840,prefix_name15:75.615840):14.006880,prefix_name16:89.622720):12.435660)xx:77.941620);'  # nopep8

        self.assertEqual(newick2, expected_newick2)


class Methods(unittest.TestCase):

    def test_reorder(self):
        """Test reordering of tree children."""
        infile = StringIO("((a,b),(c,d));")
        tree = read_tree(infile)

        infile = StringIO("((d,c),(b,a));")
        tree2 = read_tree(infile)

        hashtree1 = tree.get_one_line_newick()
        hashtree2 = tree2.get_one_line_newick()
        self.assertTrue(hashtree1 != hashtree2)

        reorder_tree(tree, tree2)
        hashtree1 = tree.get_one_line_newick()
        hashtree2 = tree2.get_one_line_newick()
        self.assertEqual(hashtree1, hashtree2)


class Draw(unittest.TestCase):
    def test_draw_tree(self):
        """Test tree drawing"""
        text = "((A:10,B:1):5,(C:2,D:3):5);"
        tree = treelib.parse_newick(text)
        out = StringIO()
        treelib.draw_tree(tree, scale=1, spacing=2, out=out,
                          labelOffset=-1, minlen=1)
        drawing = out.getvalue()
        expected = '''\
      /---------  A
 /----+
 |    \  B
++
 |    /-  C
 \----+
      \--  D
'''
        self.assertEqual(drawing, expected)

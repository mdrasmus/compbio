
import StringIO
import unittest

from compbio import arglib
from rasmus import util
from rasmus.common import izip
from rasmus.rplotting import rp
from rasmus.rplotting import rplot_end
from rasmus.rplotting import rplot_start
from rasmus.testing import make_clean_dir


class Arg (unittest.TestCase):

    def test_sample_coal_recomb(self):
        rho = 1.5e-8  # recomb/site/gen
        l = 2000      # length of locus
        k = 10        # number of lineages
        n = 2*10000   # effective popsize
        r = rho * l   # recomb/locus/gen
        nsamples = 10000

        samples = [arglib.sample_coal_recomb(k, n, r)
                   for i in range(nsamples)]
        events = dict(
            (event, count / float(nsamples))
            for event, count in util.hist_dict(util.cget(samples, 0)).items())
        expected = {'coal': 0.88146, 'recomb': 0.11854}

        for key, value in events.items():
            self.assertAlmostEqual(value, expected[key], places=2)

    def test_sample_lineages(self):
        """lineage over time"""

        outdir = 'test/tmp/test_arglib/Arg_test_sample_lineages/'
        make_clean_dir(outdir)

        rho = 1.5e-8  # recomb/site/gen
        l = 5000      # length of locus
        k = 60        # number of lineages
        n = 2*10000   # effective popsize
        r = rho * l   # recomb/locus/gen

        rplot_start(outdir + '/plot.pdf')
        rp.plot([1, 40000], [1, k],  t="n", log="xy")
        times, events = arglib.sample_coal_recomb_times(k, n, r)
        lineages = list(arglib.lineages_over_time(k, events))
        rp.lines(times, lineages)
        rplot_end()

    def test_read_write(self):
        """Read and write an ARG"""

        rho = 1.5e-8   # recomb/site/gen
        l = 10000      # length of locus
        k = 10         # number of lineages
        n = 2*10000    # effective popsize

        arg = arglib.sample_arg(k, n, rho, 0, l)
        # round ages and pos for easy equality
        for node in arg:
            node.age = round(node.age)
            node.pos = round(node.pos)

        stream = StringIO.StringIO()
        arglib.write_arg(stream, arg)
        stream.seek(0)
        arg2 = arglib.read_arg(stream)

        self.assertTrue(arg.equal(arg2))

    def test_local_trees(self):

        rho = 1.5e-8   # recomb/site/gen
        l = 10000      # length of locus
        k = 10         # number of lineages
        n = 2*1e4      # effective popsize

        arg = arglib.sample_arg(k, n, rho, 0, l)
        blocks1 = util.cget(arglib.iter_local_trees(arg, 200, 1200), 0)
        blocks2 = list(arglib.iter_recomb_blocks(arg, 200, 1200))
        self.assertEqual(blocks1, blocks2)

    def test_marginal_leaves(self):

        rho = 1.5e-8   # recomb/site/gen
        l = 10000      # length of locus
        k = 10         # number of lineages
        n = 2*10000    # effective popsize

        arg = arglib.sample_arg(k, n, rho, 0, l)

        for (start, end), tree in arglib.iter_local_trees(arg):
            arglib.remove_single_lineages(tree)
            mid = (start + end) / 2.0
            for node in tree:
                a = set(tree.leaves(node))
                b = set(arglib.get_marginal_leaves(arg, node, mid))
                self.assertEqual(a, b)

    #----------------------------------
    # SPRs

    def test_iter_sprs(self):

        rho = 1.5e-8   # recomb/site/gen
        l = 100000     # length of locus
        k = 6          # number of lineages
        n = 2*10000    # effective popsize

        arg = arglib.sample_arg(k, n, rho, 0, l)

        for a, b in izip(arglib.iter_arg_sprs(arg),
                         arglib.iter_arg_sprs_simple(arg)):
            self.assertEqual(a, b)

    def test_iter_sprs_leaves(self):

        rho = 1.5e-8   # recomb/site/gen
        l = 100000     # length of locus
        k = 40         # number of lineages
        n = 2*10000    # effective popsize

        arg = arglib.sample_arg(k, n, rho, 0, l)

        for a, b in izip(arglib.iter_arg_sprs(arg, use_leaves=True),
                         arglib.iter_arg_sprs_simple(arg, use_leaves=True)):
            a[1][0].sort()
            a[2][0].sort()
            b[1][0].sort()
            b[2][0].sort()
            self.assertEqual(a, b)

    def test_iter_sprs_time(self):

        rho = 1.5e-8   # recomb/site/gen
        l = 100000     # length of locus
        k = 40         # number of lineages
        n = 2*10000    # effective popsize

        arg = arglib.sample_arg(k, n, rho, 0, l)

        x = list(arglib.iter_arg_sprs(arg))
        x = list(arglib.iter_arg_sprs_simple(arg))
        x = list(arglib.iter_arg_sprs(arg, use_leaves=True))
        x = list(arglib.iter_arg_sprs_simple(arg, use_leaves=True))
        x

    def test_iter_sprs_remove_thread(self):

        rho = 1.5e-8   # recomb/site/gen
        l = 100000     # length of locus
        k = 6          # number of lineages
        n = 2*10000    # effective popsize

        arg = arglib.sample_arg(k, n, rho, 0, l)
        remove_chroms = set("n%d" % (k-1))
        keep = [x for x in arg.leaf_names() if x not in remove_chroms]
        arg = arg.copy()
        arglib.subarg_by_leaf_names(arg, keep)

        for a, b in izip(arglib.iter_arg_sprs(arg),
                         arglib.iter_arg_sprs_simple(arg)):
            self.assertEqual(a, b)

    def test_smcify_arg(self):

        rho = 1.5e-8   # recomb/site/gen
        l = 100000     # length of locus
        k = 6          # number of lineages
        n = 2*10000    # effective popsize

        arg = arglib.sample_arg(k, n, rho, 0, l)
        arg = arglib.smcify_arg(arg)

        for pos, (rnode, rtime), (cnode, ctime) in arglib.iter_arg_sprs(arg):
            self.assertNotEqual(rnode, cnode)

    def test_smcify_arg_remove_thread(self):

        rho = 1.5e-8   # recomb/site/gen
        l = 100000      # length of locus
        k = 6         # number of lineages
        n = 2*10000    # effective popsize

        arg = arglib.sample_arg(k, n, rho, 0, l)
        remove_chroms = set("n%d" % (k-1))
        keep = [x for x in arg.leaf_names() if x not in remove_chroms]
        arg = arg.copy()
        arglib.subarg_by_leaf_names(arg, keep)
        arg = arglib.smcify_arg(arg)

    #----------------------------
    # SMC sampling

    def test_sample_arg_smc(self):
        """Sample an ARG using the SMC process"""

        length = 10000     # length of locus
        k = 5              # number of lineages
        n = 1e4            # effective popsize
        rho = 1.5e-8       # recomb/site/gen

        arg = arglib.sample_arg_smc(k, n, rho, 0, length)
        arglib.assert_arg(arg)

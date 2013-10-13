
import StringIO
import unittest

from compbio import arglib
from rasmus.common import ilen, izip
from rasmus.stats import mean
from rasmus import util
from ramsus.util import plot

from rasmus.common import rp


class Arg (unittest.TestCase):

    def test_sample_coal_recomb(self):
        rho = 1.5e-8  # recomb/site/gen
        l = 1000      # length of locus
        k = 10        # number of lineages
        n = 2*10000   # effective popsize
        r = rho * l   # recomb/locus/gen

        event, t = arglib.sample_coal_recomb(k, n, r)
        print event, t

    def test_sample_lineages(self):
        """lineage over time"""

        rho = 1.5e-8  # recomb/site/gen
        l = 5000      # length of locus
        k = 60        # number of lineages
        n = 2*10000   # effective popsize
        r = rho * l   # recomb/locus/gen

        rp.plot([1, 40000], [1, k],  t="n", log="xy")
        times, events = arglib.sample_coal_recomb_times(k, n, r)
        lineages = list(arglib.lineages_over_time(k, events))
        rp.lines(times, lineages)

    def test_read_write(self):
        """Read and write an ARG"""

        rho = 1.5e-8   # recomb/site/gen
        l = 10000      # length of locus
        k = 10         # number of lineages
        n = 2*10000    # effective popsize

        arg = arglib.sample_arg(k, n, rho, 0, l)
        stream = StringIO.StringIO()
        arglib.write_arg(stream, arg)
        stream.seek(0)
        arglib.read_arg(stream)

    def test_local_trees(self):

        rho = 1.5e-8   # recomb/site/gen
        l = 10000      # length of locus
        k = 10         # number of lineages
        n = 2*1e4      # effective popsize

        arg = arglib.sample_arg(k, n, rho, 0, l)

        print arglib.get_recomb_pos(arg)
        blocks1 = util.cget(arglib.iter_local_trees(arg, 200, 1200), 0)
        blocks2 = list(arglib.iter_recomb_blocks(arg, 200, 1200))

        print blocks1
        print blocks2
        assert blocks1 == blocks2

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
                print mid, a, b
                assert a == b

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
            print a, b
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
            print a, b
            self.assertEqual(a, b)

    def test_iter_sprs_time(self):

        rho = 1.5e-8   # recomb/site/gen
        l = 100000      # length of locus
        k = 40         # number of lineages
        n = 2*10000    # effective popsize

        arg = arglib.sample_arg(k, n, rho, 0, l)

        util.tic("arglib.iter_arg_sprs")
        x = list(arglib.iter_arg_sprs(arg))
        util.toc()

        util.tic("arglib.iter_arg_sprs_simple")
        x = list(arglib.iter_arg_sprs_simple(arg))
        util.toc()

        util.tic("arglib.iter_arg_sprs use_leaves=True")
        x = list(arglib.iter_arg_sprs(arg, use_leaves=True))
        util.toc()

        util.tic("arglib.iter_arg_sprs_simple use_leaves=True")
        x = list(arglib.iter_arg_sprs_simple(arg, use_leaves=True))
        util.toc()
        x

    def test_iter_sprs_dsmc(self):

        import arghmm

        rho = 1.5e-8   # recomb/site/gen
        l = 100000      # length of locus
        k = 40         # number of lineages
        n = 2*10000    # effective popsize
        times = arghmm.get_time_points(ntimes=20, maxtime=160000)

        arg = arghmm.sample_arg_dsmc(k, n, rho, 0, l, times=times)

        for a, b in izip(arglib.iter_arg_sprs(arg),
                         arglib.iter_arg_sprs_simple(arg)):
            print a, b
            self.assertEqual(a, b)

    def test_iter_sprs_remove_thread(self):

        rho = 1.5e-8   # recomb/site/gen
        l = 100000      # length of locus
        k = 6         # number of lineages
        n = 2*10000    # effective popsize

        arg = arglib.sample_arg(k, n, rho, 0, l)
        remove_chroms = set("n%d" % (k-1))
        keep = [x for x in arg.leaf_names() if x not in remove_chroms]
        arg = arg.copy()
        arglib.subarg_by_leaf_names(arg, keep)

        for a, b in izip(arglib.iter_arg_sprs(arg),
                         arglib.iter_arg_sprs_simple(arg)):
            print a, b
            self.assertEqual(a, b)

    def test_smcify_arg(self):

        rho = 1.5e-8   # recomb/site/gen
        l = 100000      # length of locus
        k = 6         # number of lineages
        n = 2*10000    # effective popsize

        arg = arglib.sample_arg(k, n, rho, 0, l)
        arg = arglib.smcify_arg(arg)

        for pos, (rnode, rtime), (cnode, ctime) in arglib.iter_arg_sprs(arg):
            print (rnode, rtime), (cnode, rtime), rnode == cnode
            assert rnode != cnode

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

        length = 100000    # length of locus
        k = 2              # number of lineages
        n = 1e4            # effective popsize
        rho = 1.5e-8 * 20  # recomb/site/gen

        arg = arglib.sample_arg_smc(k, n, rho, 0, length)
        arg2 = arglib.smcify_arg(arglib.sample_arg(k, n, rho, 0, length))
        arglib.assert_arg(arg)

        plot([x.age for x in arglib.iter_visible_recombs(arg)],
             main="smc")
        plot([x.age for x in arglib.iter_visible_recombs(arg2)],
             main="coal_recomb")

    def test_sample_arg_smc_cmp(self):
        """Sample an ARG using the SMC process and compare it"""

        k = 10             # number of lineages
        n = 1e4            # effective popsize
        rho = 1.5e-8 * 20  # recomb/site/gen
        length = int(500e3 / 20)  # length of locus

        x = []
        y = []
        for i in range(1, 100):
            print i
            arg = arglib.sample_arg_smc(k, n, i/100. * rho, 0, length)
            arg2 = arglib.smcify_arg(
                arglib.sample_arg(k, n, i/100. * rho, 0, length))
            x.append(ilen(arglib.iter_visible_recombs(arg)))
            y.append(ilen(arglib.iter_visible_recombs(arg2)))

        p = plot(x, y, main="recombs", xlab="smc", ylab="coal_recomb")
        p.plot([0, max(x)], [0, max(x)], style="lines")

    def test_sample_arg_smc_cmp2(self):
        """Sample an ARG using the SMC process and compare it"""

        length = 1000      # length of locus
        k = 4              # number of lineages
        n = 1e4            # effective popsize
        rho = 1.5e-8 * 20  # recomb/site/gen

        x = []
        y = []
        for i in range(100):
            print i
            arg = arglib.sample_arg_smc(k, n, rho, 0, length)
            arg2 = arglib.smcify_arg(arglib.sample_arg(k, n, rho, 0, length))
            x.append(max(x.age for x in arg))
            y.append(max(x.age for x in arg2))

        p = plot(x, y, main="maxage", xlab="smc", ylab="coal_recomb",
                 xlog=10, ylog=10, xmin=1000, ymin=1000)
        p.plot([100, max(x)], [100, max(x)], style="lines")

    def test_sample_arg_smc_cmp3(self):
        """Sample an ARG using the SMC process and compare it"""

        length = 1000      # length of locus
        k = 4              # number of lineages
        n = 1e4            # effective popsize
        rho = 1.5e-8 * 20  # recomb/site/gen

        x = []
        y = []
        for i in range(100):
            print i
            arg = arglib.sample_arg_smc(k, n, rho, 0, length)
            arg2 = arglib.sample_arg(k, n, rho, 0, length)
            x.append(arglib.arglen(arg))
            y.append(arglib.arglen(arg2))

        p = plot(x, y, main="arglen", xlab="smc", ylab="coal_recomb")
        p.plot([0, max(x)], [0, max(x)], style="lines")

    def test_sample_arg_smc_cmp4(self):
        """Sample an ARG using the SMC process and compare it"""

        length = 1000      # length of locus
        k = 4              # number of lineages
        n = 1e4            # effective popsize
        rho = 1.5e-8 * 20  # recomb/site/gen

        x = []
        y = []
        for i in range(100):
            print i
            arg = arglib.smcify_arg(
                arglib.sample_arg_smc(k, n, rho, 0, length))
            arg2 = arglib.sample_arg(k, n, rho, 0, length)
            x.append(mean(x.age for x in arg if x.event == "recomb"))
            y.append(mean(x.age for x in arg2 if x.event == "recomb"))

        p = plot(x, y, main="avg recomb age", xlab="smc", ylab="coal_recomb")
        p.plot([0, max(x)], [0, max(x)], style="lines")


from math import exp
import unittest

from compbio import coal
from compbio import phylo

from rasmus import stats
from rasmus import treelib
from rasmus.gnuplot import Gnuplot
from rasmus.gnuplot import plot
from rasmus.gnuplot import plotfunc
from rasmus.gnuplot import plotdistrib
from rasmus.tablelib import Table
from rasmus.testing import eq_sample_pdf
from rasmus.testing import fequal
from rasmus.testing import fequals
from rasmus.testing import make_clean_dir
from rasmus.util import cumsum
from rasmus.util import distrib
from rasmus.util import frange
from rasmus.util import safelog


# TODO: re-enable broken tests


#=============================================================================
# test coalescence times (normal, censored, bounded)


class Coal (unittest.TestCase):

    def test_prob_coal(self):

        outdir = 'test/tmp/test_coal/Coal_test_prob_coal/'
        make_clean_dir(outdir)

        k = 2
        n = 1000
        p = Gnuplot()
        p.enableOutput(False)
        p.plotfunc(lambda t: coal.prob_coal(t, k, n), 0, 4000, 10,
                   ymin=0)

        # draw single coal samples
        x = [coal.sample_coal(k, n) for i in xrange(200)]
        plotdistrib(x, 40, plot=p)
        p.enableOutput(True)
        p.save(outdir + 'plot.png')

        eq_sample_pdf(x, lambda t: coal.prob_coal(t, k, n), 40)

    def test_prob_coal2(self):

        outdir = 'test/tmp/test_coal/Coal_test_prob_coal2/'
        make_clean_dir(outdir)

        k = 2
        n = 1000
        p = Gnuplot()
        p.enableOutput(False)
        p.plotfunc(lambda t: coal.prob_coal(t, k, n), 0, 4000, 10,
                   ymin=0)

        # draw single coal samples
        x = [coal.sample_coal(k, n) for i in xrange(200)]
        plotdistrib(x, 40, plot=p)
        p.enableOutput(True)
        p.save(outdir + 'plot.png')

        eq_sample_pdf(x, lambda t: coal.prob_coal(t, k, n), 40)

    def test_prob_coal_cond_counts_simple(self):

        # when we condition on b=1, it is the same as the bounded coal
        # PDF.
        # prob_coal_cond_counts is actually a more general version of
        # prob_bounded_coal

        outdir = 'test/tmp/test_coal/Coal_test_prob_coal_counys_simple/'
        make_clean_dir(outdir)

        a = 5
        b = 1
        t = 2000
        n = 1000
        p = Gnuplot()
        p.enableOutput(False)

        for x in frange(0, 2000, 10):
            y = coal.prob_coal_cond_counts_simple(x, a, b, t, n)
            y2 = coal.prob_bounded_coal(x, a, n, t)
            self.assertAlmostEqual(y, y2)

    def test_prob_coal_cond_counts(self):

        # match against a simpler slower implementation
        a = 5
        b = 3
        t = 2000
        n = 1000

        for x in frange(0, 2000, 10):
            y = coal.prob_coal_cond_counts(x, a, b, t, n)
            y2 = coal.prob_coal_cond_counts_simple(x, a, b, t, n)
            self.assertAlmostEqual(y, y2)

    def test_prob_coal_cond_counts2_simple(self):

        # test coalescent pdf when conditioned on future lineage counts

        a = 5
        for b in xrange(2, a):
            t = 500
            n = 1000

            # draw single coal samples using rejection sampling
            x2 = []
            for i in xrange(1000):
                while True:
                    times = coal.sample_coal_times(a, n)
                    if times[a-b-1] < t and (b == 1 or times[a-b]) > t:
                        break
                x2.append(times[0])

            eq_sample_pdf(
                x2, lambda x: coal.prob_coal_cond_counts_simple(
                    x, a, b, t, n), 40)

    def test_prob_coal_cond_counts2(self):

        # test coalescent pdf when conditioned on future lineage counts

        a = 5
        for b in xrange(2, a):
            t = 500
            n = 1000

            # draw single coal samples using rejection sampling
            x2 = []
            for i in xrange(1000):
                while True:
                    times = coal.sample_coal_times(a, n)
                    if times[a-b-1] < t and (b == 1 or times[a-b]) > t:
                        break
                x2.append(times[0])

            eq_sample_pdf(x2, lambda x: coal.prob_coal_cond_counts(
                x, a, b, t, n), 40)

    def test_cdf_coal_cond_counts(self):

        # test coalescent pdf when conditioned on future lineage counts

        outdir = 'test/tmp/test_coal/Coal_test_cdf_coal_cond_counts/'
        make_clean_dir(outdir)

        a = 5
        for b in xrange(2, a):
            t = 500
            n = 1000
            p = Gnuplot()
            p.enableOutput(False)
            p.plotfunc(lambda x: coal.cdf_coal_cond_counts(
                x, a, b, t, n), 0, t, 10)

            # draw single coal samples using rejection sampling
            s = []
            for i in xrange(1000):
                while True:
                    times = coal.sample_coal_times(a, n)
                    if times[a-b-1] < t and (b == 1 or times[a-b] > t):
                        break
                s.append(times[0])

            x2, y2 = stats.cdf(s)
            p.plot(x2, y2, style='lines')
            p.enableOutput(True)
            p.save(outdir + 'plot-%d.png' % b)

            eq_sample_pdf(
                x2, lambda x: coal.prob_coal_cond_counts(x, a, b, t, n), 40)

    def test_sample_coal_cond_counts(self):

        # test coalescent pdf when conditioned on future lineage counts

        a = 5
        for b in xrange(2, a):
            t = 500
            n = 1000

            # draw single coal samples using rejection sampling
            s = [coal.sample_coal_cond_counts(a, b, t, n)
                 for i in xrange(1000)]

            eq_sample_pdf(s, lambda x: coal.prob_coal_cond_counts(
                x, a, b, t, n), 40)

    def test_sample_coal_tree(self):
        n = 1000
        tree = coal.sample_coal_tree(10, n)
        treelib.assert_tree(tree)

    def test_sample_censored_coal(self):
        n = 1000
        tree, lineages = coal.sample_censored_coal_tree(
            10, n, 300, capped=True)
        treelib.assert_tree(tree)

    def test_prob_coal_counts(self):
        self.assertAlmostEqual(
            coal.prob_coal_counts(2, 2, 100, 1000),
            0.904837418036)
        self.assertAlmostEqual(
            coal.prob_coal_counts(5, 2, 100, 1000),
            0.0184034834527)

    def test_prob_mrca(self):
        n = 1000
        k = 50

        x = [coal.sample_coal_times(k, n)[-1] for i in xrange(10000)]
        eq_sample_pdf(x, lambda t: coal.prob_mrca(t, k, n), 40)

    def test_cdf_mrca(self):

        outdir = 'test/tmp/test_coal/Coal_test_cdf_mrca/'
        make_clean_dir(outdir)

        n = 1000
        k = 6
        step = 10
        x = list(frange(0, 5000, step))
        y = [coal.prob_mrca(i, k, n) * step for i in x]
        y2 = cumsum(y)
        y3 = [coal.cdf_mrca(t, k, n) for t in x]

        p = Gnuplot()
        p.enableOutput(False)
        p.plot(x, y2, style="lines")
        p.plot(x, y3, style="lines")
        p.enableOutput(True)
        p.save(outdir + 'plot.png')

        eq_sample_pdf(x, lambda t: coal.cdf_mrca(t, k, n), 40)

    def test_cdf_mrca2(self):
        n = 1000
        k = 6
        step = 10
        x = list(frange(0, 5000, step))
        y = [coal.prob_mrca(i, k, n) * step for i in x]
        y2 = cumsum(y)
        y3 = [coal.cdf_mrca(t, k, n) for t in x]
        y4 = [coal.prob_coal_counts(k, 1, t, n) for t in x]

        fequals(y2, y3, eabs=.01)
        fequals(y3, y4)


class _BoundedCoal (unittest.TestCase):

    def test_plot_bounded_coal(self):
        n = 1000
        k = 6
        T = 500

        # plots should differ
        p = plotfunc(lambda t: coal.prob_bounded_coal(t, k, n, T),
                     0, 1000, 10)
        p.plotfunc(lambda t: coal.prob_coal(t, k, n),
                   0, 1000, 10)

    def test_cdf_coal_bounded(self):
        n = 1000
        k = 4
        t = 500
        step = .1
        x = list(frange(0, 500, step))
        y = [coal.prob_bounded_coal(i, k, n, t) * step for i in x]
        y2 = cumsum(y)
        y3 = [coal.cdf_bounded_coal(i, k, n, t) for i in x]

        p = plot(x, y2, style="lines")
        p.plot(x, y3, style="lines")
        p.plot([0, 500], [1, 1], style="lines")

        fequals(y2, y3, eabs=.01)

    def test_sample_bounded_coal2(self):
        n = 1000
        k = 2
        T = 800

        d = [coal.sample_bounded_coal2(n, T) for i in xrange(2000)]

        p = plotdistrib(d, 40)
        p.plotfunc(lambda t: coal.prob_bounded_coal(t, k, n, T), 0, T, T/200.0)

        eq_sample_pdf(d, lambda t: coal.prob_bounded_coal(t, k, n, T), 40)

    def test_sample_bounded_coal(self):
        n = 1000
        k = 5
        T = 500

        d = [coal.sample_bounded_coal(k, n, T) for i in xrange(2000)]

        p = plotdistrib(d, 40)
        p.plotfunc(lambda t: coal.prob_bounded_coal(t, k, n, T), 0, T, T/200)

        eq_sample_pdf(d, lambda t: coal.prob_bounded_coal(t, k, n, T), 40)

    def test_plot_prob_bounded_coal(self):
        n = 1000
        k = 4
        t = 800
        alltimes = []

        # sample times
        for i in xrange(5000):
            while True:
                times = [0]
                for j in xrange(k, 1, -1):
                    times.append(times[-1] + coal.sample_coal(j, n))
                    if times[-1] >= t:
                        break
                if times[-1] < t:
                    break
            alltimes.append(times)

        p = Gnuplot()
        for i in range(1, 2):
            x, y = distrib([q[i] - q[i-1] for q in alltimes], width=20)
            p.plot(x, y, style="lines", xmax=500)

        x = list(frange(0, 500, 10))
        #for i in range(1, 2): #k):
        y2 = [coal.prob_bounded_coal(j, k, n, t) for j in x]
        p.plot(x, y2, style="lines", xmax=500)

        fequals(y, y2, rel=.05, eabs=.01)

    def test_fast_sample_bounded_coal(self):

        # sample bounded coal times efficiently
        n = 1000
        k = 5
        t = 500
        alltimes = []

        # sample times
        for i in xrange(5000):
            while True:
                times = [0]
                for j in xrange(k, 1, -1):
                    times.append(times[-1] + coal.sample_coal(j, n))
                    if times[-1] >= t:
                        break
                if times[-1] < t:
                    break
            alltimes.append(times)

        p = Gnuplot()
        for i in range(1, k):
            x, y = distrib([q[i] - q[i-1] for q in alltimes], width=30)
            p.plot(x, y, style="lines", xmax=500)
        p.enableOutput(True)
        p.replot()

        # sample times efficently
        alltimes2 = []
        for i in xrange(5000):
            times = [0]
            for j in xrange(k, 1, -1):
                times.append(times[-1] +
                             coal.sample_bounded_coal(j, n, t-times[-1]))
            alltimes2.append(times)

        #p = Gnuplot()
        for i in range(1, k):
            x, y = distrib([q[i] - q[i-1] for q in alltimes2], width=30)
            p.plot(x, y, style="lines", xmax=500)
        p.enableOutput(True)
        p.replot()


class CoalCounts (unittest.TestCase):

    def test_coal_counts(self):
        b = 1
        t = 1000.0
        n = 1000

        for a in xrange(1, 100):
            i = coal.prob_coal_counts(a, b, t, n)
            j = coal.cdf_mrca(t, a, n)
            fequal(i, j)

        for a in xrange(1, 10):
            i = sum(coal.prob_coal_counts(a, b, t, n)
                    for b in xrange(1, a+1))
            fequal(i, 1.0)

    def test_coal_counts2(self):
        b = 3
        t = 1000.0
        n = 1000

        for b in xrange(1, 10):
            for a in xrange(b, 10):
                i = coal.prob_coal_counts(a, b, t, n)
                j = coal.prob_coal_counts_slow(a, b, t, n)
                fequal(i, j)


#=============================================================================
# multicoal

def _test_multicoal_tree(stree, n, nsamples):
    """test multicoal_tree"""
    tops = {}

    for i in xrange(nsamples):
        tree, recon = coal.sample_multicoal_tree(stree, n,
                                                 namefunc=lambda x: x)
        top = phylo.hash_tree(tree)
        tops.setdefault(top, [0, tree, recon])[0] += 1

    tab = Table(headers=["top", "simple_top", "percent", "prob"])
    for top, (num, tree, recon) in tops.items():
        tree2 = tree.copy()
        treelib.remove_single_children(tree2)

        print phylo.hash_tree(tree2)
        print phylo.hash_tree(stree)

        tab.add(top=top,
                simple_top=phylo.hash_tree(tree2),
                percent=num/float(nsamples),
                prob=exp(coal.prob_multicoal_recon_topology(
                    tree, recon, stree, n)))
    tab.sort(col="prob", reverse=True)

    return tab, tops


def _test_bounded_multicoal_tree(stree, n, T, nsamples):
    """test multicoal_tree"""
    tops = {}

    for i in xrange(nsamples):

        # use rejection sampling
        #tree, recon = coal.sample_bounded_multicoal_tree_reject(
        #    stree, n, T, namefunc=lambda x: x)

        # sample tree
        tree, recon = coal.sample_bounded_multicoal_tree(
            stree, n, T, namefunc=lambda x: x)

        top = phylo.hash_tree(tree)
        tops.setdefault(top, [0, tree, recon])[0] += 1

    tab = Table(headers=["top", "simple_top", "percent", "prob"])
    for top, (num, tree, recon) in tops.items():
        tree2 = tree.copy()
        treelib.remove_single_children(tree2)
        tab.add(top=top,
                simple_top=phylo.hash_tree(tree2),
                percent=num/float(nsamples),
                prob=exp(coal.prob_bounded_multicoal_recon_topology(
                    tree, recon, stree, n, T)))
    tab.sort(col="prob", reverse=True)

    return tab, tops


class _MultiCoal (unittest.TestCase):

    def test_1(self):
        """Test multicoal_tree on simple 4 species tree"""

        stree = treelib.parse_newick(
            "((A:1000, B:1000):500, (C:700, D:700):800);")
        n = 500
        nsamples = 1000
        tab, tops = _test_multicoal_tree(stree, n, nsamples)
        print repr(tab[:20].get(cols=["simple_top", "percent", "prob"]))
        a, b = tab[:20].cget("percent", "prob")
        fequals(a, b, eabs=.05)

    def test_flies(self):

        stree = treelib.parse_newick("""(
  (
    (
      (
        (
          (
            dmel:5.32,
            (
              dsec:1.89,
              dsim:1.89
            ):3.43
          ):5.91,
          (
            dere:8.57,
            dyak:8.57
          ):2.66
        ):42.17,
        dana:53.40
      ):2.40,
      (
        dpse:1.37,
        dper:1.37
      ):54.43
    ):6.69,
    dwil:62.49
  ):1.02,
  (
    (
      dmoj:32.74,
      dvir:32.74
    ):4.37,
    dgri:37.11
  ):26.40
);""")
        for node in stree:
            node.dist *= 1e6 * 10
        n = 10e6
        nsamples = 5000
        tab, tops = _test_multicoal_tree(stree, n, nsamples)
        print repr(tab[:20].get(cols=["simple_top", "percent", "prob"]))
        a, b = tab[:20].cget("percent", "prob")
        fequals(a, b, eabs=.05)


class BMC (unittest.TestCase):

    def test_cdf_bmc_simple(self):

        # test cdf mrca BMC
        stree = treelib.parse_newick("((A:1000, B:1000):500, C:1500);")
        n = 1000
        gene_counts = dict.fromkeys(stree.leaf_names(), 1)
        T = 2000

        self.assertAlmostEqual(
            exp(coal.cdf_mrca_bounded_multicoal(gene_counts, T, stree, n)),
            0.27719726132)

    def test_cdf_bmc(self):

        # test cdf mrca BMC
        stree = treelib.parse_newick(
            "((A:1000, B:1000):500, (C:700, D:700):800);")
        n = 1000
        gene_counts = dict.fromkeys(stree.leaf_names(), 1)
        T = 2000

        p = exp(coal.cdf_mrca_bounded_multicoal(gene_counts, T, stree, n))

        nsamples = 5000
        c = 0
        for i in xrange(nsamples):
            tree, recon = coal.sample_multicoal_tree(stree, n)
            if treelib.get_tree_timestamps(tree)[tree.root] < T:
                c += 1
        p2 = c / float(nsamples)

        fequal(p, p2, .05)

    def _test_recon(self):

        # test multicoal_tree on simple 4 species tree
        stree = treelib.parse_newick(
            "((A:1000, B:1000):500, (C:700, D:700):800);")
        n = 500
        T = 2000
        nsamples = 10000
        tab, tops = _test_bounded_multicoal_tree(stree, n, T, nsamples)
        print repr(tab[:20].get(cols=["simple_top", "percent", "prob"]))

        a, b = tab[:20].cget("percent", "prob")
        fequals(a, b, rel=.05, eabs=.005)

    def test_top(self):

        outdir = 'test/tmp/test_coal/BMC_test_top/'
        make_clean_dir(outdir)

        stree = treelib.parse_newick(
            "(((A:200, E:200):800, B:1000):500, (C:700, D:700):800);")
        n = 500
        T = 2000
        nsamples = 4000

        # compare top hist with simpler rejection sampling
        tops = {}
        tops2 = {}

        for i in xrange(nsamples):
            # use rejection sampling
            tree, recon = coal.sample_bounded_multicoal_tree_reject(
                stree, n, T, namefunc=lambda x: x)

            # sample tree
            tree2, recon2 = coal.sample_bounded_multicoal_tree(
                stree, n, T, namefunc=lambda x: x)

            top = phylo.hash_tree(tree)
            top2 = phylo.hash_tree(tree2)

            tops.setdefault(top, [0, tree, recon])[0] += 1
            tops.setdefault(top2, [0, tree2, recon2])

            tops2.setdefault(top2, [0, tree2, recon2])[0] += 1
            tops2.setdefault(top, [0, tree, recon])

        keys = tops.keys()
        x = [safelog(tops[i][0], default=0) for i in keys]
        y = [safelog(tops2[i][0], default=0) for i in keys]

        self.assertTrue(stats.corr(x, y) > .9)

        p = Gnuplot()
        p.enableOutput(False)
        p.plot(x, y)
        p.plot([min(x), max(x)], [min(x), max(x)], style="lines")
        p.enableOutput(True)
        p.save(outdir + 'plot.png')

    '''
    def test_prob_coal(self):

        # test the waiting time PDF of a coal in the BMC
        stree = treelib.parse_newick(
            "(((A:700, B:700):300, E:1000):500, (C:700, D:700):800);")
        n = 500
        T = 1700
        nsamples = 5000

        t = 200
        u = stree.nodes["A"].parent.parent
        ucount = 3
        utime = 1000
        gene_counts = None

        c = 1.0 - exp(coal.prob_no_coal_bmc(
            u, utime, ucount, gene_counts, T, stree, n,
            sroot=None, sleaves=None, stimes=None,
            tree=None, recon=None))

        def pdf(t):
            return c * exp(coal.prob_coal_bmc(
                t, u, utime, ucount, gene_counts, T, stree, n,
                sroot=None, sleaves=None, stimes=None,
                tree=None, recon=None))


        print "sum", integrate(pdf, 0, 500, 1)
        print "coal", c


        p = plotfunc(pdf, 0, 500, 10)
        plotfunc(lambda t: coal.prob_coal(t, ucount, n), 0, 500, 10, plot=p)
        pause()

        draw_tree_names(stree, maxlen=8)
    '''

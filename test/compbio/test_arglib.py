

import sys
import StringIO

from rasmus.common import *
from rasmus.testing import *
from compbio import arglib



#=============================================================================


    
class Arg (unittest.TestCase):

    def test_sample_coal_recomb(self):
        rho = 1.5e-8 # recomb/site/gen
        l = 1000     # length of locus
        k = 10       # number of lineages
        n = 2*10000  # effective popsize
        r = rho * l  # recomb/locus/gen
    
        event, t = arglib.sample_coal_recomb(k, n, r)    
        print event, t


    def test_sample_lineages(self):
        """lineage over time"""

        rho = 1.5e-8 # recomb/site/gen
        l = 5000     # length of locus
        k = 60       # number of lineages
        n = 2*10000  # effective popsize
        r = rho * l  # recomb/locus/gen

        rp.plot([1, 40000], [1, k],  t="n", log="xy")
        times, events = arglib.sample_coal_recomb_times(k, n, r)
        lineages = list(arglib.lineages_over_time(k, events))
        rp.lines(times, lineages)
        pause()

    def test_read_write(self):
        """Read and write an ARG"""
        
        rho = 1.5e-8   # recomb/site/gen
        l = 10000      # length of locus
        k = 10         # number of lineages
        n = 2*10000    # effective popsize
        r = rho * l    # recomb/locus/gen

        arg = arglib.sample_arg(k, n, rho, 0, l)
        stream = StringIO.StringIO()
        arglib.write_arg(stream, arg)
        stream.seek(0)
        arg2 = arglib.read_arg(stream)

    def test_marginal_leaves(self):

        rho = 1.5e-8   # recomb/site/gen
        l = 10000      # length of locus
        k = 10         # number of lineages
        n = 2*10000    # effective popsize
        r = rho * l    # recomb/locus/gen

        arg = arglib.sample_arg(k, n, rho, 0, l)
        
        for (start, end), tree in arglib.iter_tree_tracks(arg):
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
        l = 100000      # length of locus
        k = 6         # number of lineages
        n = 2*10000    # effective popsize
        r = rho * l    # recomb/locus/gen

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
        r = rho * l    # recomb/locus/gen

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
        r = rho * l    # recomb/locus/gen

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


    def test_iter_sprs_dsmc(self):

        import arghmm

        rho = 1.5e-8   # recomb/site/gen
        l = 100000      # length of locus
        k = 40         # number of lineages
        n = 2*10000    # effective popsize
        r = rho * l    # recomb/locus/gen
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
        r = rho * l    # recomb/locus/gen

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
        r = rho * l    # recomb/locus/gen

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
        r = rho * l    # recomb/locus/gen

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
        k = 2             # number of lineages
        n = 1e4           # effective popsize
        rho = 1.5e-8 * 20 # recomb/site/gen

        arg = arglib.sample_arg_smc(k, n, rho, 0, length)
        arg2 = arglib.smcify_arg(arglib.sample_arg(k, n, rho, 0, length))
        arglib.assert_arg(arg)
        
        p = plot([x.age for x in arglib.iter_visible_recombs(arg)],
                 main="smc")
        p2 = plot([x.age for x in arglib.iter_visible_recombs(arg2)],
                  main="coal_recomb")
        
        pause()


    def test_sample_arg_smc_cmp(self):
        """Sample an ARG using the SMC process and compare it"""
        
        k = 10             # number of lineages
        n = 1e4            # effective popsize
        rho = 1.5e-8 * 20  # recomb/site/gen
        length = int(500e3 / 20) # length of locus

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
        
        pause()


    def test_sample_arg_smc_cmp2(self):
        """Sample an ARG using the SMC process and compare it"""
        
        length = 1000    # length of locus
        k = 4             # number of lineages
        n = 1e4           # effective popsize
        rho = 1.5e-8 * 20 # recomb/site/gen

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
        
        pause()


    def test_sample_arg_smc_cmp3(self):
        """Sample an ARG using the SMC process and compare it"""
        
        length = 1000    # length of locus
        k = 4             # number of lineages
        n = 1e4           # effective popsize
        rho = 1.5e-8 * 20 # recomb/site/gen

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
        
        pause()
        

    def test_sample_arg_smc_cmp4(self):
        """Sample an ARG using the SMC process and compare it"""
        
        length = 1000    # length of locus
        k = 4             # number of lineages
        n = 1e4           # effective popsize
        rho = 1.5e-8 * 20 # recomb/site/gen

        x = []
        y = []
        for i in range(100):
            print i
            arg = arglib.smcify_arg(arglib.sample_arg_smc(k, n, rho, 0, length))
            arg2 = arglib.sample_arg(k, n, rho, 0, length)
            x.append(mean(x.age for x in arg if x.event == "recomb"))
            y.append(mean(x.age for x in arg2 if x.event == "recomb"))

        p = plot(x, y, main="avg recomb age", xlab="smc", ylab="coal_recomb")
        p.plot([0, max(x)], [0, max(x)], style="lines")
        
        pause()

    

# lineages over time is deterministic with high k
if 0:
    rho = 1.5e-8   # recomb/site/gen
    l = 5000       # length of locus
    k = 1000       # number of lineages
    n = 2*10000    # effective popsize
    r = rho * l    # recomb/locus/gen

    rplot_start("test/figures/many_lineages.pdf")
    rplot("plot", [1, 40000], [1, k],  t="n", log="xy")
    for i in xrange(20):
        times, events = arglib.sample_coal_recomb_times(k, n, r)
        lineages = list(arglib.lineages_over_time(k, events))
        rplot("lines", times, lineages)

    # change Ne
    n = .8 * 2*10000
    for i in xrange(20):
        times, events = arglib.sample_coal_recomb_times(k, n, r)
        lineages = list(arglib.lineages_over_time(k, events))
        rplot("lines", times, lineages, col="red")

    rplot_end(True)


# test postorder
if 0:
    rho = 1.5e-8   # recomb/site/gen
    l = 5000     # length of locus
    k = 10       # number of lineages
    n = 2*10000  # effective popsize
    r = rho * l  # recomb/locus/gen
    
    times, events = arglib.sample_coal_recomb_times(k, n, r)
    arg = arglib.make_arg_from_times(k, times, events)

    seen = set()
    for node in arg.postorder():
        for child in node.children:
            assert child in seen, node
        seen.add(node)
    print "DONE"


# test sample_arg
if 0:
    rho = 1.5e-8   # recomb/site/gen
    l = 10000       # length of locus
    k = 10         # number of lineages
    n = 2*10000    # effective popsize
    r = rho * l    # recomb/locus/gen

    arg = arglib.sample_arg(k, n, rho, 0, l)
    arglib.write_arg(sys.stdout, arg)


# test read/write
if 0:
    rho = 1.5e-8   # recomb/site/gen
    l = 10000      # length of locus
    k = 10         # number of lineages
    n = 2*10000    # effective popsize
    r = rho * l    # recomb/locus/gen

    from compbio.vis import argvis

    arg = arglib.sample_arg(k, n, rho, 0, l)
    arg.write()

    stream = StringIO.StringIO()
    arglib.write_arg(stream, arg)
    stream.seek(0)
    arg2 = arglib.read_arg(stream)
    argvis.show_arg(arg)
    argvis.show_arg(arg2)
    

# test coal_regions
if 0:
    regions1 = [(0, 10), (10, 20), (30, 31)]
    regions2 = [(0, 13), (15, 25)]
    regions3 = [(12, 40)]
    for start, end, count in arglib.count_region_overlaps(
        regions1, regions2, regions3):
        print start, end, count


# visualize arg
if 0:
    rho = 1.5e-8   # recomb/site/gen
    l = 10000       # length of locus
    k = 100         # number of lineages
    n = 2*10000    # effective popsize
    r = rho * l    # recomb/locus/gen

    #arg = arglib.sample_arg(k, n, rho, 0, l)
    arg = arglib.sample_arg_smc(k, n, rho, 0, l)
    #arglib.show_arg(arg)

    print len([node for node in arglib.remove_single_lineages(
               arg.get_marginal_tree(l+1))
               if len(node.children) == 2])


# extract marginal trees
if 0:
    rho = 1e-8   # recomb/site/gen
    l = 5000     # length of locus
    k = 10       # number of lineages
    n = 2*10000  # effective popsize
    r = rho * l  # recomb/locus/gen
    
    times, events = arglib.sample_coal_recomb_times(k, n, r)
    arg = arglib.make_arg_from_times(k, times, events)
    arg.set_recomb_pos()
    arg.set_ancestral()

    tree = arg.get_tree(.9)
    draw_tree_names(tree, maxlen=5, minlen=5)


# extract marginal trees more directly
if 0:
    rho = 1e-8   # recomb/site/gen
    l = 5000     # length of locus
    k = 10       # number of lineages
    n = 2*10000  # effective popsize
    r = rho * l  # recomb/locus/gen
    
    arg = arglib.sample_arg(k, n, rho, 0, l)
    tree = arg.get_tree(.9)
    draw_tree_names(tree, maxlen=5, minlen=5)
    

# histogram of how many nodes to expect in arg
if 0:
    rho = 1.5e-8   # recomb/site/gen
    l = 5000     # length of locus
    k = 60       # number of lineages
    n = 2*10000  # effective popsize
    r = rho * l  # recomb/locus/gen

    ntimes = []
    for i in xrange(100):
        times, events = arglib.sample_coal_recomb_times(k, n, r)
        ntimes.append(len(times))
    plothist(ntimes, width=1)


# how often do we have recombination in the locus?
if 0:
    rho = 1.5e-8   # recomb/site/gen
    l = 1000       # length of locus
    k = 10         # number of lineages
    n = 2*10000    # effective popsize
    r = rho * l    # recomb/locus/gen

    nrecombs = []
    for i in xrange(10000):
        arg = arglib.sample_arg(k, n, rho, 0, l)
        nrecombs.append(len([x for x in arg if x.event == "recomb"]))

    plothist(nrecombs, 50)


# how long are non-recombining loci?
if 0:
    rho = 1.5e-8   # recomb/site/gen
    l = 10000      # length of locus
    k = 10         # number of lineages
    n = 2*10000    # effective popsize
    r = rho * l    # recomb/locus/gen

    nonrecomblen = []
    nonrecomblenv = []
    num = []
    ks = range(2, 500, 10)
    for k in ks:
        lens = []
        for j in xrange(50):
            arg = arglib.sample_arg(k, n, rho, 0, l)
            blocks = arglib.iter_recomb_blocks(arg)
            lens.extend(r[1] - r[0] for r in blocks)

        nonrecomblen.append(mean(lens))
        num.append(len(lens))
        print k, nonrecomblen[-1]

    p = plot(ks, nonrecomblen, style="lines")



# sample mutations for 10 individual
if 0:
    l = 1000       # length of locus
    k = 100        # number of lineages
    n = 2*10000    # effective popsize
    rho = 1.5e-8   # recomb/site/gen
    r = rho * l    # recomb/locus/gen
    mu = 2.5e-8    # mut/site/gen
    u = mu * l     # mut/locus/gen

    nmut = []
    for i in xrange(300):
        arg = arglib.sample_arg(k, n, rho, 0, l)
    
        mut = arglib.sample_mutations(arg, u)
        nmut.append(len(mut))

    plothist(nmut)
    print mean(nmut)




# how does total branch relate to k
def treelen(k, times):
    tot = 0.0
    t = 0.0
    for i in xrange(len(times)):
        tot += (times[i] - t) * k
        t = times[i]
        k -= 1
    return tot
if 0:
    ks = range(2, 1000, 10)
    n = 2*10000    # effective popsize

    totlen = []
    for k in ks:
        totlen.append(mean(treelen(k, coal.sample_coal_times(k, n))
                           for i in xrange(50)))
    plot(ks, totlen, style="lines")
        
        

# determine fraction of ancestral sequence
if 0:
    rho = 1e-8   # recomb/site/gen
    l = 5000     # length of locus
    k = 10       # number of lineages
    n = 2*10000  # effective popsize
    r = rho * l  # recomb/locus/gen

    nnodes = []
    noanc = []
    for i in xrange(100):
        times, events = arglib.sample_coal_recomb_times(k, n, r)
        arg = arglib.make_arg_from_times(k, times, events)
        arg.set_recomb_pos()
        arg.set_ancestral()        

        nnodes.append(len(times))
        noanc.append(len([node for node in arg
                          if len(node.data["ancestral"]) == 0]))

    plothist(noanc)
    plot(nnodes, noanc)



# test that all ancestral sequence appears in marginal trees and vice versa
if 0:
    rho = 1e-8   # recomb/site/gen
    l = 1000     # length of locus
    k = 100       # number of lineages
    n = 2*10000  # effective popsize
    r = rho * l  # recomb/locus/gen

    noanc = []
    for i in xrange(100):
        arg = arglib.sample_arg(k, n, rho, 0, l)
        anc = set(node.name for node in arg
                  if len(node.data["ancestral"]) > 0)
        rpos = sorted([node.pos for node in arg if node.event == "recomb"])
        mnodes = set(node.name for tree in arglib.iter_marginal_trees(arg)
                     for node in tree)

        arg.prune(False)
        prune = set(arg.nodes.keys())

        assert (anc - mnodes) == (mnodes - anc), (anc - mnodes, mnodes - anc)
        assert mnodes == prune
        
    print "DONE"


# what is the frequency distribution of mutations
if 0:
    l = 100000     # length of locus
    k = 100        # number of lineages
    n = 2*10000    # effective popsize
    rho = 1.5e-8   # recomb/site/gen
    mu = 2.5e-8    # mut/site/gen
    
    f = []
    for j in xrange(1):
        print j
        arg = arglib.sample_arg(k, n, rho, 0, l)
        mut = arglib.sample_arg_mutations(arg, mu)
        for node, parent, pos, t in mut:
            f.append(float(ilen(arglib.get_marginal_leaves(arg, node, pos))) / k)

    plotdistrib(f, 50)

    

# what is the frequency distribution of mutations
if 0:
    l = 100000     # length of locus
    k = 100        # number of lineages
    n = 2*10000    # effective popsize
    rho = 1.5e-8   # recomb/site/gen
    mu = 2.5e-8    # mut/site/gen
    
    f = []
    for j in xrange(1):
        print j
        arg = arglib.sample_arg(k, n, rho, 0, l)
        mut = arglib.sample_arg_mutations(arg, mu)
        for node, parent, pos, t in mut:
            f.append(float(ilen(arglib.get_marginal_leaves(arg, node, pos))) / k)

    plotdistrib(f, 50)

    


def sample_expected_coal_time(k, n):
    return (2.0 * n) / (k*(k-1))

# is there a simple way to polarize?
# TODO:
if 0:
    k = 100
    n = 2*1e4

    x = []
    y = []
    for i in xrange(20):
        times = coal.sample_coal_times(k, n)
        for j, t in enumerate(times):
            x.append(j)
            y.append(t)
    p = plot(x, y, ylog=10, ymin=.1)

    x = []
    y = []
    times = (sample_expected_coal_time(k2, n)
             for k2 in xrange(k, 1, -1))
    t = 0
    for j, t2 in enumerate(times):
        x.append(j)
        t += t2
        y.append(t)
    p.plot(x, y, ylog=10, ymin=.1)
    


def iter_tree_random_postorder(tree, node=None):
    
    if node is None:
        node = tree.root
    children = list(node.children)
    random.shuffle(children)
    stack = [[node, children]]

    while len(stack) > 0:
        node, children = stack[-1]

        if len(children) > 0:
            child = children.pop()
            children2 = list(child.children)
            random.shuffle(children2)
            stack.append([child, children2])
        else:
            yield node
            stack.pop()

if 0:
    tree = treelib.parse_newick("(((A,B),C),(D,E));")

    print
    treelib.draw_tree_names(tree, minlen=5)

    print [x.name for x in iter_tree_random_postorder(tree)]



# how long are split tracks?
if 0:
    rho = 1.5e-8   # recomb/site/gen
    l = 100000     # length of locus
    k = 100        # number of lineages
    n = 2*10000    # effective popsize
    r = rho * l    # recomb/locus/gen
    
    track = []
    for j in xrange(1):
        arg = arglib.sample_arg(k, n, rho, 0, l)
        print "ARG sim"

        blocks = arglib.iter_recomb_blocks(arg)
        split_regions = defaultdict(lambda: [])

        # find regions for splits
        for block, tree in izip(blocks, arglib.iter_marginal_trees(arg)):
            for node in tree:
                if len(node.children) != 2 or node == tree.root: continue
                split = tuple(sorted(tree.leaf_names(node)))
                if len(split) == 1: continue
                split_regions[split].append(block)

        # union regions
        for split, regions in split_regions.iteritems():
            split_regions[split] =  [
                (min(r[0] for r in g), max(r[1] for r in g))
                for g in arglib.groupby_overlaps(regions)]

        for split, regions in split_regions.iteritems():
            for region in regions:
                track.append(region[1] - region[0])

    plothist(track,50)
    print "mean track", mean(track)


# visualize split tracks
if 0:
    rho = 1.5e-8   # recomb/site/gen
    l = 10000      # length of locus
    k = 100        # number of lineages
    n = 2*10000    # effective popsize
    r = rho * l    # recomb/locus/gen
    mu = 2.5e-8    # mut/site/gen
    u = mu * l     # mut/locus/gen
    
    arg = arglib.sample_arg_smc(k, n, rho, 0, l)
    arg.set_ancestral()
    mut = arglib.sample_mutations(arg, u)
    mut.sort(key=lambda x: x[2])
    print "sim", len(mut)

    blocks = arglib.iter_recomb_blocks(arg)
    split_regions = defaultdict(lambda: [])
    split_nodes = defaultdict(lambda: [])

    # find regions for splits
    for block, tree in izip(blocks, arglib.iter_marginal_trees(arg)):
        for node in tree:
            if len(node.children) != 2 or node.children[0] == node.children[1]:
                continue
            split = tuple(sorted(tree.leaf_names(node)))
            if len(split) == 1: continue
            split_regions[split].append(block)
            split_nodes[split].append((block[0], node))
    
    # sort splits
    tracks = split_regions.items()
    tracks.sort(key=lambda x: (x[1][0][0], x[1][-1][1])) # sort by first appearance
    split_lookup = list2lookup(cget(tracks, 0))

    split_mut = set()
    mut_splits = []
    for pos, split in arglib.iter_mutation_splits(arg, mut):
        split_mut.add(split)
        mut_splits.append((pos, split))


    def click(i):
        def func():
            split = tracks[i][0]
            print i, split, split_nodes[tracks[i][0]][0]
            print [
                (min(r[0] for r in g), max(r[1] for r in g))
                for g in arglib.groupby_overlaps(tracks[i][1])]
        return func

    import summon
    from summon.core import *
    win = summon.Window()

    # draw splits
    for i, (split, regions) in enumerate(tracks):
        col = (0, 1, 0) if split in split_mut else (1, 1, 1, .5)
        win.add_group(group(color(*col),
                group(* (lines(r[0], i, r[1], i) for r in regions)),
                hotspot("click", regions[0][0], i-.5, regions[-1][1], i+.5,
                        click(i))))
    
    # draw mutations
    for node, parent, pos, t in mut:
        split = tuple(sorted(x.name for x in 
                             arglib.get_marginal_leaves(arg, node, pos)))
        if len(split) == 1: continue
        win.add_group(arglib.draw_mark(pos, split_lookup[split], col=(0,0,1)))
        
    print len(split_mut), len(split_lookup)


    # infer splits from mut_splits
    tracks2 = []
    cur = {}
    for pos, split in mut_splits:
        for split2 in cur.keys():
            if not arglib.is_split_compatible(split, split2):
                tracks2.append(((cur[split2], pos), split2))
                del cur[split2]
        if split not in cur:
            cur[split] = pos
    for split2 in cur:
        tracks2.append(((cur[split2], l), split2))
    
    # draw infered splits
    #for r, split in tracks2:
    #    col = (1, 0, 0)
    #    y = split_lookup[split]
    #    win.add_group(lines(color(*col), r[0], y+.2, r[1], y+.2))

    



# how many mutations per split track?
if 0:
    rho = 1.5e-8   # recomb/site/gen
    l = 100000      # length of locus
    k = 100       # number of lineages
    n = 2*10000    # effective popsize
    r = rho * l    # recomb/locus/gen
    mu = 2.5e-8    # mut/site/gen
    u = mu * l     # mut/locus/gen

    nmuts = []
    for j in xrange(1):
        arg = arglib.sample_arg(k, n, rho, 0, l)
        muts = arglib.sample_mutations(arg, u)
        print "ARG sim"

        split_muts = defaultdict(lambda: 0)
        
        for node, parent, pos, t in muts:
            split = tuple(sorted(x.name for x in 
                                 arglib.get_marginal_leaves(arg, node, pos)))
            if len(split) == 1: continue
            split_muts[split] += 1

    #rp.hist(split_muts.values(), main="", xlab="# mutations")
    plotdistrib(split_muts.values(), low=0, width=1)
    print mean(split_muts.values())


# visualize conflict graph (remove redundant conflicts)
if 0:
    rho = 1.5e-8   # recomb/site/gen
    mu = 2.5e-8    # mut/site/gen
    l = 100000     # length of locus
    k = 100        # number of lineages
    n = 2*10000    # effective popsize


    nmuts = []
    for j in xrange(1):
        arg = arglib.sample_arg(k, n, rho, 0, l)
        muts = arglib.sample_arg_mutations(arg, mu)
        print "ARG sim"

        leaves = tuple(sorted(arg.leaf_names()))
        splits = []
        for pos, split in arglib.iter_mutation_splits(arg, muts):
            splits.append(split)
        
        conflicts = {}
        
        for split in splits:
            conflicts[split] = sets.UnionFind([split])
        for split1 in splits:
            for split2 in splits:
                if not arglib.is_split_compatible_unpolar(split1, split2, leaves):
                    conflicts[split1].union(conflicts[split2])
        components = set(s.root() for s in conflicts.itervalues())
        print [len(x) for x in components]

        cmat = []
        for split1 in splits:
            cmat.append([])
            for split2 in splits:
                cmat[-1].append(int(not arglib.is_split_compatible(split1, 
                                                                   split2)))

        #from rasmus import cluto
        #partids = cluto.cluster(cmat, 10, "vcluster")
        #tree = cluto.clusterTree(cmat, 10, prog="scluster")
        #perm = cluto.reorderTree(tree, cmat, prog="scluster")
        #perm = sortindex(map(sum,cmat))

        #from summon import matrix
        #v = matrix.show_dmat(cmat, rperm=perm, cperm=perm)
        v = matrix.show_dmat(cmat)




def get_unique_conflicts(splits):
    """Ignore redundant conflicts"""

    n = len(splits)
    conflicts = []
    right = [n] * n # nearest conflict to right
    left = [-1] * n  # nearest conflict to left

    leaves = tuple(sorted(arg.leaf_names()))

    # visit conflict edges in increase itervals
    for k in range(1, n):
        for i in range(n - k):
            j = i + k
            if right[i] < left[j]:
                # redundant, skip
                continue
            if not arglib.is_split_compatible_unpolar(splits[i], splits[j], leaves):
                # conflict
                conflicts.append((i, j))
                right[i] = min(right[i], j)
                left[j] = max(left[j], i)

    return conflicts, left, right



# visualize conflict graph
if 0:
    rho = 1.5e-8   # recomb/site/gen
    mu = 2.5e-8    # mut/site/gen
    l = 100000     # length of locus
    k = 100        # number of lineages
    n = 2*10000    # effective popsize


    nmuts = []
    for j in xrange(1):
        arg = arglib.sample_arg(k, n, rho, 0, l)
        muts = arglib.sample_arg_mutations(arg, mu)
        print "ARG sim"

        splits = []
        for pos, split in arglib.iter_mutation_splits(arg, muts):
            splits.append(split)
        
        conflicts, left, right = get_unique_conflicts(splits)
        imat = []
        for i, j in conflicts:
            imat.append((i, j, 1))
            imat.append((j, i, 1))

        from summon import matrix
        v = matrix.show_imat(imat)




# how far does one need to go along an alignment to find an incompatible mut?
if 0:
    rho = 1.5e-8   # recomb/site/gen
    l = 50000      # length of locus
    k = 1000         # number of lineages
    n = 2*10000    # effective popsize
    r = rho * l    # recomb/locus/gen
    mu = 2.5e-8    # mut/site/gen
    u = mu * l     # mut/locus/gen

    nmuts = []
    for j in xrange(1):
        arg = arglib.sample_arg(k, n, rho, 0, l)
        muts = arglib.sample_mutations(arg, u)        
        muts.sort(key=lambda x: x[2])

        blocks = arglib.iter_recomb_blocks(arg)
        splits = set()
        incompats = []

        for node, parent, pos, t in muts:
            if node.is_leaf(): continue
            split = tuple(sorted(x.name for x in 
                                 arglib.get_marginal_leaves(arg, node, pos)))

            for split2 in splits:
                if not arglib.is_split_compatible(split, split2):
                    incompats.append(pos)
                    splits.clear()
                    break
            splits.add(split)

    print incompats
    lens = []
    for i in xrange(len(incompats) - 1):
        lens.append(incompats[i+1] - incompats[i])
    plothist(lens)
    print mean(lens)


#=============================================================================
# visualization

# visualize arg and mutations
if 0:
    l = 10000      # length of locus
    k = 100         # number of lineages
    n = 2*10000    # effective popsize
    rho = 1.5e-8   # recomb/site/gen
    r = rho * l    # recomb/locus/gen
    mu = 2.5e-8    # mut/site/gen
    u = mu * l     # mut/locus/gen

    arg = arglib.sample_arg(k, n, rho, 0, l)
    mut = arglib.sample_mutations(arg, u)

    def minlog(x, low=10):
        return log(max(x, low))
    #ymap = lambda x: x
    ymap = minlog

    import summon
    from summon.core import *
    win = summon.Window()
    x = 0
    step = 2
    treewidth = ilen(arg.leaves()) + step

    def trans_camera(win, x, y):
        v = win.get_visible()
        win.set_visible(v[0]+x, v[1]+y, v[2]+x, v[3]+y, "exact")

    win.set_binding(input_key("]"), lambda : trans_camera(win, treewidth, 0))
    win.set_binding(input_key("["), lambda : trans_camera(win, -treewidth, 0))

    blocks = list(arglib.iter_recomb_blocks(arg))

    for tree, block in izip(arglib.iter_marginal_trees(arg), blocks):
        pos = block[0]
        print pos
        
        layout = arglib.layout_arg(tree, yfunc=ymap)
        win.add_group(translate(x, 0, color(1,1,1),
                                arglib.draw_tree(tree, layout)))

        # mark responsible recomb node
        for node in tree:
            if pos != 0.0 and node.pos == pos:
                nx, ny = layout[node]
                win.add_group(arglib.draw_mark(x + nx, ny))

        # draw mut
        for node, parent, mpos, t in mut:
            if (node.name in tree and node.name != tree.root.name and
                block[0] < mpos < block[1]):
                nx, ny = layout[tree[node.name]]
                win.add_group(arglib.draw_mark(x + nx, ymap(t), col=(0,0,1)))

                for leaf in arg.leaves(node):
                    nx, ny = layout[tree[leaf.name]]
                    win.add_group(arglib.draw_mark(x + nx, ymap(0),col=(0,1,0)))
                
        x += treewidth

    # home
    win.home("exact")





# visualize arg and mutations
if 0:
    l = 10000      # length of locus
    k = 100         # number of lineages
    n = 2*10000    # effective popsize
    rho = 1.5e-8   # recomb/site/gen
    r = rho * l    # recomb/locus/gen
    mu = 2.5e-8    # mut/site/gen
    u = mu * l     # mut/locus/gen

    #times, events = arglib.sample_coal_recomb_times(k, n, r)
    #arg = arglib.make_arg_from_times(k, times, events)
    #arg.set_recomb_pos(0, l)
    #arg.set_ancestral()
    #arg.prune()
    arg = arglib.sample_arg(k, n, rho, 0, l)
    mut = arglib.sample_mutations(arg, u)

    #reload(arglib)
    #layout = arglib.layout_arg(arg)
    #win = arglib.show_arg(arg, mut=mut)

    cmap = rainbowColorMap(low=0, high=l)

    def minlog(x, low=10):
        return log(max(x, low))
    #ymap = lambda x: x
    ymap = minlog

    import summon
    from summon.core import *
    win = summon.Window()
    x = 0
    step = 2
    treewidth = ilen(arg.leaves()) + step

    def trans_camera(win, x, y):
        v = win.get_visible()
        win.set_visible(v[0]+x, v[1]+y, v[2]+x, v[3]+y, "exact")

    win.set_binding(input_key("]"), lambda : trans_camera(win, treewidth, 0))
    win.set_binding(input_key("["), lambda : trans_camera(win, -treewidth, 0))

    blocks = arglib.iter_recomb_blocks(arg)

    for tree, block in izip(arglib.iter_marginal_trees(arg), blocks):
        pos = block[0]
        print pos
        
        leaves = sorted((x for x in tree.leaves()), key=lambda x: x.name)
        layout = arglib.layout_arg(tree, leaves, yfunc=ymap)
        win.add_group(translate(x, 0, color(1,1,1,.2),
                                arglib.draw_tree(tree, layout)))

        # mark responsible recomb node
        for node in tree:
            if pos != 0.0 and node.pos == pos:
                nx, ny = layout[node]
                win.add_group(arglib.draw_mark(x + nx, ny))

        # draw mut
        for node, parent, mpos, t in mut:
            if (node.name in tree and node.name != tree.root.name and
                block[0] < mpos < block[1]):
                nx, ny = layout[tree[node.name]]
                win.add_group(arglib.draw_mark(x + nx, ymap(t), col=(0,0,1)))
            if node.name in tree and tree[node.name].parents:
                nx, ny = layout[tree[node.name]]
                py = layout[tree[node.name].parents[0]][1]
                start = arg[node.name].data["ancestral"][0][0]
                win.add_group(lines(color(0,1,0), #*cmap.get(start)), 
                                    x+nx, ny, x+nx, py,
                                    color(1,1,1)))
            
                
        x += treewidth

    win.set_visible(* win.get_root().get_bounding() + ("exact",))





#=============================================================================
# understand recombinations and compatiability


# visualize compatiability
if 0:
    rho = 1.5e-8   # recomb/site/gen
    l = 100000      # length of locus
    k = 60         # number of lineages
    n = 2*10000    # effective popsize
    r = rho * l    # recomb/locus/gen
    mu = 2.5e-8    # mut/site/gen
    u = mu * l     # mut/locus/gen
    
    arg = arglib.sample_arg(k, n, rho, 0, l)
    mut = arglib.sample_mutations(arg, u)
    mut.sort(key=lambda x: x[2])
    print "sim", len(mut)

    # remove uninformative muts (singletons)
    splits, mut = zip(* (a for a in ((arglib.get_mutation_split(arg, m), m)
                                  for m in mut)
                         if len(a[0]) > 1))

    # get recomb points
    rpos = arglib.get_recomb_pos(arg)
    part = []
    i = 0
    for node, parent, pos, t in mut:
        while i < len(rpos) and pos > rpos[i]:
            i += 1
        part.append(i)
    
    splits = [arglib.get_mutation_split(arg, m) for m in mut]
    
    compat = []
    for i in xrange(len(splits)):
        print "mut", i
        for j in xrange(i+1, len(mut)):
            rel = arglib.split_relation(splits[i], splits[j])
            if rel == "disjoint":
                compat.append((i, j, 1))
            elif rel in ("parent", "child", "equal"):
                compat.append((i, j, 2))
            elif rel == "conflict":
                pass
            else:
                raise Exception("unknown rel" + rel)

    import summon.matrix
    v = summon.matrix.show_imat(compat, rpart=part, cpart=part, style="quads")



# visualize redundant conflicts
if 0:
    rho = 1.5e-8   # recomb/site/gen
    l = 100000      # length of locus
    k = 60         # number of lineages
    n = 2*10000    # effective popsize
    r = rho * l    # recomb/locus/gen
    mu = 2.5e-8    # mut/site/gen
    u = mu * l     # mut/locus/gen
    
    arg = arglib.sample_arg(k, n, rho, 0, l)
    mut = arglib.sample_mutations(arg, u)
    mut.sort(key=lambda x: x[2])
    print "sim", len(mut)

    # remove uninformative muts (singletons)
    splits, mut = zip(* (a for a in ((arglib.get_mutation_split(arg, m), m)
                                  for m in mut)
                         if len(a[0]) > 1))

    # get recomb points
    rpos = arglib.get_recomb_pos(arg)
    part = []
    i = 0
    for node, parent, pos, t in mut:
        while i < len(rpos) and pos > rpos[i]:
            i += 1
        part.append(i)
    
    splits = [arglib.get_mutation_split(arg, m) for m in mut]
    
    compat = []
    for i in xrange(len(splits)):
        print "mut", i
        compat.append([0] * (i+1))
        for j in xrange(i+1, len(splits)):
            rel = arglib.split_relation(splits[i], splits[j])
            if rel == "disjoint":
                compat[-1].append(1)
            elif rel in ("parent", "child", "equal"):
                compat[-1].append(1)
            elif rel == "conflict":
                compat[-1].append(2)
            else:
                raise Exception("unknown rel" + rel)


    # remove redudant conflicts
    for i in xrange(len(splits)):
        for j in xrange(i+1, len(splits)):
            followl = 0
            followr = 0
            dist = j - i
            for k in xrange(1, dist):
                if compat[i][i+k] == 2:
                    followl = 1 # follow left
                if compat[i+k][j] == 2:
                    followr = 1 # follow right
                if followl + followr == 2:
                    compat[i][j] = 0
                    break


    import summon.matrix
    v = summon.matrix.show_dmat(compat, rpart=part, cpart=part, style="quads")

    

# visualize conflict graph
if 0:
    rho = 1.5e-8   # recomb/site/gen
    l = 100000      # length of locus
    k = 60         # number of lineages
    n = 2*10000    # effective popsize
    r = rho * l    # recomb/locus/gen
    mu = 2.5e-8    # mut/site/gen
    u = mu * l     # mut/locus/gen
    
    arg = arglib.sample_arg(k, n, rho, 0, l)
    mut = arglib.sample_mutations(arg, u)
    mut.sort(key=lambda x: x[2])
    print "sim", len(mut)

    mat = []
    last = None
    lookup = {}
    for tree in arglib.iter_marginal_trees(arg):
        splits = set(arglib.iter_tree_splits(tree))
        for split in splits:
            if split not in lookup:
                lookup[split] = len(lookup)

        if last:
            remove = last - splits
            add = splits - last
            for r in remove:
                for a in add:
                    mat.append((lookup[r], lookup[a], 1))
            
        last = set(splits)
        
    v = summon.matrix.show_imat(mat)
    



# visualize supporting mutations
if 0:
    rho = 1.5e-8   # recomb/site/gen
    l = 100000      # length of locus
    k = 60         # number of lineages
    n = 2*10000    # effective popsize
    r = rho * l    # recomb/locus/gen
    mu = 2.5e-8    # mut/site/gen
    u = mu * l     # mut/locus/gen
    
    arg = arglib.sample_arg(k, n, rho, 0, l)
    mut = arglib.sample_mutations(arg, u)
    mut.sort(key=lambda x: x[2])
    print "sim", len(mut)

    # remove uninformative muts (singletons)
    splits, mut = zip(* (a for a in ((arglib.get_mutation_split(arg, m), m)
                                  for m in mut)
                         if len(a[0]) > 1))


    split_loc = defaultdict(lambda:[])
    for m in mut:
        split = arglib.get_mutation_split(arg, m)
        split_loc[split].append(m[2])

    # get recomb points
    rpos = arglib.get_recomb_pos(arg)
    
    import summon
    from summon.core import *
    win = summon.Window()
    
    for i, r in enumerate(rpos):
        splits1 = set(arglib.iter_tree_splits(arg.get_marginal_tree(r-.01)))
        splits2 = set(arglib.iter_tree_splits(arg.get_marginal_tree(r+.01)))
        remove = splits1 - splits2
        add = splits2 - splits1
        all = list(add | remove)

        x = [loc for split in all for loc in split_loc[split]]

        if len(x) > 0:
            win.add_group(lines(color(1,1,1),
                                min(min(x), r), i, max(max(x), r), i))
            for mx in x:
                win.add_group(arglib.draw_mark(mx, i, col=(0,0,1)))
        
        win.add_group(arglib.draw_mark(r, i, col=(1,0,0)))
        
        

    win.home("exact")


#=============================================================================
if __name__ == "__main__":
    test_main()

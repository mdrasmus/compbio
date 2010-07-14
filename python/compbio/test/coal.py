
from rasmus.common import *
from rasmus import stats
from rasmus.testing import *

from compbio import coal
reload(coal)


#=============================================================================
# test coalescence times (normal, fixed, bounded)

# plot single coal PDF
if 0:
    k = 2
    n = 1000
    p = plotfunc(lambda t: coal.prob_coal(t, k, n), 0, 4000, 10,
                 ymin=0)

    # draw single coal samples
    x = [coal.sample_coal(k, n) for i in xrange(200)]
    plotdistrib(x, 40, plot=p)


# draw sampled coal tree
if 0:
    n = 1000
    tree = coal.sample_coal_tree(10, n)
    draw_tree(tree, scale=.01)


# draw sampled coal tree with fixed time
if 0:
    n = 1000
    tree, lineages = coal.sample_coal_tree_fixed(10, n, 300, capped=True)
    print lineages
    draw_tree(tree, scale=.01)
    #show_tree(tree)

# test prob counts
if 0:
    print coal.prob_coal_counts(2, 2, 48e6, 20000)


# test MRCA
if 0:
    n = 1000
    k = 50

    x = [coal.sample_coal_times(k, n)[-1] for i in xrange(10000)]
    p = plotdistrib(x, width=100)

    p.plotfunc(lambda i: coal.prob_mrca(i, k, n),
               0, max(x), max(x) / 100.0)


# test CDF MRCA
if 0:
    n = 1000
    k = 6
    step = 10
    x = list(frange(0, 5000, step))
    y = [coal.prob_mrca(i, k, n) * step for i in x]
    y2 = cumsum(y)
    y3 = [coal.cdf_mrca(t, k, n) for t in x]
    
    p = plot(x, y2, style="lines")
    p.plot(x, y3, style="lines")


# test CDF MRCA
if 0:
    n = 1000
    k = 6
    step = 10
    x = list(frange(0, 5000, step))
    y = [coal.prob_mrca(i, k, n) * step for i in x]
    y2 = cumsum(y)
    y3 = [coal.cdf_mrca(t, k, n) for t in x]
    y4 = [coal.prob_coal_counts(k, 1, t, n) for t in x]
    
    p = plot(x, y2, style="lines")
    p.plot(x, y3, style="lines")
    p.plot(x, y4, style="lines")


# plot pdf of bounded coal
if 0:
    n = 1000
    k = 6
    T = 500
    
    #x = [coal.sample_coal_bounded(k, n, t)
    #     for i in xrange(500)]
    #plothist(x)

    p = plotfunc(lambda t: coal.prob_coal_bounded(t, k, n, T),
                 0, 1000, 10)
    p.plotfunc(lambda t: coal.prob_coal(t, k, n),
               0, 1000, 10)


# test CDF coal_bounded
if 0:
    n = 1000
    k = 4
    t = 500
    step = .1
    x = list(frange(0, 500, step))
    y = [coal.prob_coal_bounded(i, k, n, t) * step for i in x]
    y2 = cumsum(y)
    y3 = [coal.cdf_coal_bounded(i, k, n, t) for i in x]
    
    p = plot(x, y2, style="lines")
    p.plot(x, y3, style="lines")
    p.plot([0, 500], [1, 1], style="lines")



# test sample bounded coal (k=2)
if 0:
    n = 1000
    k = 2
    T = 800

    tic("sample")
    d = [coal.sample_coal_bounded2(n, T) for i in xrange(2000)]
    toc()

    p = plotdistrib(d, 40)
    p.plotfunc(lambda t: coal.prob_coal_bounded(t, k, n, T), 0, T, t/200)
    

# test sample bounded coal
if 0:
    n = 1000
    k = 5
    T = 500

    tic("sample")
    d = [coal.sample_coal_bounded(k, n, T) for i in xrange(2000)]
    toc()

    p = plotdistrib(d, 40)
    p.plotfunc(lambda t: coal.prob_coal_bounded(t, k, n, T), 0, T, t/200)
    


# plot pdf of bounded coal
if 0:
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
    for i in range(1, 2): #k):
        x, y = distrib([q[i] - q[i-1] for q in alltimes], width=20)
        p.plot(x, y, style="lines", xmax=500)


    x = list(frange(0, 500, 10))
    #for i in range(1, 2): #k):
    y = [coal.prob_coal_bounded(j, k, n, t) for j in x]
    p.plot(x, y, style="lines", xmax=500)

    

# efficently sample bounded coal times
if 0:
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
    time.sleep(1)
    p.enableOutput(True)
    p.replot()


    # sample times efficently
    alltimes2 = []    
    for i in xrange(5000):
        times = [0]
        for j in xrange(k, 1, -1):
            times.append(times[-1] +
                         coal.sample_coal_bounded(j, n, t-times[-1]))
        alltimes2.append(times)

    #p = Gnuplot()
    for i in range(1, k):
        x, y = distrib([q[i] - q[i-1] for q in alltimes2], width=30)
        p.plot(x, y, style="lines", xmax=500)
    time.sleep(1)
    p.enableOutput(True)
    p.replot()

    

# test coal counts
if 0:
    b = 1
    t = 1000.0
    n = 1000

    util.tic("test coal counts")
    for a in xrange(1, 100):
        i = coal.prob_coal_counts(a, b, t, n)
        j = coal.cdf_mrca(t, a, n)
        print i, j
        fequal(i, j)
    toc()

    for a in xrange(1, 10):
        i = sum(coal.prob_coal_counts(a, b, t, n)
                  for b in xrange(1, a+1))
        print a, i
        fequal(i, 1.0)



# test coal counts
if 0:
    b = 3
    t = 1000.0
    n = 1000

    util.tic("test coal counts")
    for b in xrange(1, 10):
        for a in xrange(b, 10):
            i = coal.prob_coal_counts(a, b, t, n)
            j = coal.prob_coal_counts_slow(a, b, t, n)
            print i, j
            fequal(i, j)
    toc()


# test cdf mrca BMC
if 1:
    stree = treelib.parse_newick("((A:1000, B:1000):500, C:1500);")
    n = 1000
    gene_counts = dict.fromkeys(stree.leaf_names(), 1)
    T = 2000
                   
    print exp(coal.cdf_mrca_bounded_multicoal(gene_counts, T, stree, n))


# test cdf mrca BMC
if 0:
    stree = treelib.parse_newick("((A:1000, B:1000):500, (C:700, D:700):800);")
    n = 1000
    gene_counts = dict.fromkeys(stree.leaf_names(), 1)
    T = 2000
                   
    print exp(coal.cdf_mrca_bounded_multicoal(gene_counts, T, stree, n))

    nsamples = 5000
    c = 0
    for i in xrange(nsamples):
        tree, recon = coal.sample_multicoal_tree(stree, n)
        if treelib.get_tree_timestamps(tree)[tree.root] < T:
            c += 1
    print c / float(nsamples)

    
    


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


# plot pdf of bounded coal
if 0:
    n = 1000
    k = 45
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
    v = 1
    t = 1000.0
    n = 1000

    for u in xrange(1, 20):
        a = coal.prob_coal_counts(u, v, t, n)
        b = coal.cdf_mrca(t, u, n)
        print a, b
        fequal(a, b)

    
    for u in xrange(1, 10):
        a = sum(coal.prob_coal_counts(u, v, t, n)
                  for v in xrange(1, u+1))
        print u, a
        fequal(a, 1.0)


    

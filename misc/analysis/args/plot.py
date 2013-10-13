

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

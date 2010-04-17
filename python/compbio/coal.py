"""

   Coalescent methods

"""

from __future__ import division

from itertools import chain, izip

from math import *
import random
from rasmus import treelib, stats, util
from rasmus.symbolic import *



def prob_coal(t, k, n):
    """
    Returns the probability density of observing the first coalesce of 'k'
    individuals in a population size of 'n' at generation 't'    
    """

    # k choose 2
    k2 = k * (k-1) / 2
    k2n = k2 / n    
    return k2n * exp(- k2n * t)


def sample_coal(k, n):
    """
    Returns a smaple coalescent time for 'k' individuals in a population size 'n'
    """

    # k choose 2
    k2 = k * (k-1) / 2
    k2n = k2 / n
    return random.expovariate(k2n)



def sample_coal_tree(k, n):
    """
    Returns a simulated coalescence tree for k leaves from a population n.
    """

    times = [0]
    for j in xrange(k, 1, -1):
        times.append(times[-1] + sample_coal(j, n))

    tree = treelib.Tree()

    # initialize k children
    children = set(treelib.TreeNode(tree.new_name()) for i in xrange(k))
    for child in children:
        child.data["time"] = 0.0

    # perform k - 1 random merges
    for i in xrange(1, len(times)):
        # make new parent and merge children
        parent = treelib.TreeNode(tree.new_name())
        parent.data["time"] = times[i]
        a, b = random.sample(children, 2)
        
        tree.add_child(parent, a)
        tree.add_child(parent, b)

        # adjust children set
        children.remove(a)
        children.remove(b)
        children.add(parent)


    # set branch lengths
    for node in tree:
        if not node.parent:
            node.dist = 0.0
        else:
            node.dist = node.parent.data["time"] - node.data["time"]
        
    tree.root = children.pop()

    return tree
        



def sample_coal_tree_fixed(k, n, t, capped=False):
    """
    Returns a simulated coalescence tree for k leaves from a population n
    with a fixed maximum time t.

    The return value is the tuple (tree, lineages)
    """

    times = [0]
    for j in xrange(k, 1, -1):
        times.append(times[-1] + sample_coal(j, n))
        if times[-1] > t:
            times.pop()
            break

    tree = treelib.Tree()

    # initialize k children
    children = set(treelib.TreeNode(tree.new_name()) for i in xrange(k))
    for child in children:
        tree.add(child)
        child.data["time"] = 0.0

    # perform random merges
    for i in xrange(1, len(times)):
        # make new parent and merge children
        parent = treelib.TreeNode(tree.new_name())
        parent.data["time"] = times[i]
        a, b = random.sample(children, 2)
        
        tree.add_child(parent, a)
        tree.add_child(parent, b)

        # adjust children set
        children.remove(a)
        children.remove(b)
        children.add(parent)


    # set branch lengths
    for node in tree:
        if not node.parent:
            node.dist = t - node.data["time"]
        else:
            node.dist = node.parent.data["time"] - node.data["time"]


    # for convenience cap the tree for easy drawing/manipulation
    if capped:
        tree.make_root()
        for node in children:
            tree.add_child(tree.root, node)
            
    return tree, children



def sample_multicoal_tree(stree, n, leaf_counts=None,
                          namefunc=None):
    """
    Returns a gene tree from a multi-species coalescence process
    n -- population size (int or dict)
         If n is a dict it must map from species name to population size
    """

    # initialize vector for how many genes per extant species
    if leaf_counts is None:
        leaf_counts = dict((l, 1) for l in stree.leaf_names())

    # initialize function for generating new gene names
    if namefunc is None:
        spcounts = dict((l, 1) for l in stree.leaf_names())
        def namefunc(sp):
            name = sp + "_" + str(spcounts[sp])
            spcounts[sp] += 1
            return name

    # initialize population sizes
    if isinstance(n, (int, float)):
        popsizes = dict.fromkeys(stree.nodes.keys(), n)
    elif isinstance(n, dict):
        popsizes = n
    else:
        raise Exception("n must be a int or dict.")

    # init gene counts
    counts = dict((n.name, 0) for n in stree)
    counts.update(leaf_counts)

    # init reconciliation, events
    recon = {}

    subtrees = {}

    # loop through species tree
    for snode in stree.postorder():        
        # simulate population for one branch
        k = counts[snode.name]
        if snode.parent:
            # non basal branch
            subtree, lineages = sample_coal_tree_fixed(k, popsizes[snode.name],
                                                       snode.dist,
                                                       capped=True)
        else:
            # basal branch
            subtree = sample_coal_tree(k, popsizes[snode.name])
            lineages = subtree.root
        subtrees[snode] = (subtree, lineages)
        if snode.parent:
            counts[snode.parent.name] += len(lineages)
        for node in subtree:
            recon[node] = snode


    # stitch subtrees together
    tree = treelib.Tree()

    # add all nodes to total tree
    for subtree, lineages in subtrees.values():
        tree.merge_names(subtree)
        tree.remove(subtree.root)    
    
    for snode in stree:
        if not snode.is_leaf():
            subtree, lineages = subtrees[snode]

            # get lineages from child subtrees
            lineages2 = chain(*[subtrees[child][1]
                                for child in snode.children])

            # ensure leaves are randomly attached
            leaves = subtree.leaves()
            random.shuffle(leaves)

            # stitch leaves of the subtree to children subtree lineages
            for leaf, lineage in izip(leaves, lineages2):
                tree.add_child(leaf, lineage)


    # set root
    tree.root = subtrees[stree.root][0].root
    tree.add(tree.root)

    # name leaves
    for leaf in tree.leaves():
        tree.rename(leaf.name, namefunc(recon[leaf].name))
        
    return tree, recon




def sample_allele_freq(p, n):
    """
    Sample a new allele frequency using starting allele frequency p and
    population size n
    """
        
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0

    if p < 0.05:
        return min(float(stats.poissonvariate(p*n))/n, n)
    if p > 0.95:
        return 1.0 - min(float(stats.poissonvariate((1-p)*n))/n, n)
    
    mu = p * n
    sigma = sqrt(n * p*(1 - p))
    p1 = random.normalvariate(mu, sigma) / n

    if p1 < 0:
        return 0.0
    if p1 > 1:
        return 1.0
    return p1

    

# Legendre polynomial
def legendre_poly(n):

    """ \frac{1}{2^n n!} d^n/dx^n [(x^2 - 1)^n] """

    return simplify(('mult', ('scalar', 1.0 / (2 ** n * stats.factorial(n))),
                    derivate(('power', ('add', ('power', ('var', 'x'),
                                                         ('scalar', 2)),
                                               ('scalar', -1)),
                                       ('scalar', n)),
                             'x', n)))

def legendre(n, r):
    l = simplify(assign_vars(legendre_poly(n), {'x': r}))
    assert l[0] == 'scalar'
    return l[1]


def gegenbauer(i, r):
    return ((i * (i+1)) / float((2*i+1)*(1-r*r)) *
            (legendre(i-1, r) - legendre(i+1, r)))

def gegenbauer2(i, r):
    return ((i * (i+1)) / 2.0 * hypergeo(i+2, 1 - i, 2, (1 - r) / 2.0))


def gegenbauer3(n, a, z):

    tot = 0
    for k in xrange(int(n/2)+1):
        tot += ((-1)**k * stats.gamma(n - k + a) / (
                stats.gamma(a) * stats.factorial(k) * stats.factorial(n - 2*k))
                * ((2*z) ** (n - 2*k)))
    return tot




def prob_fix(p, n, t, k=8, esp=0.001):
    r = 1 - 2*p
    
    prob = p
    for i in xrange(1, k+1):
        term = (.5 * (-1)**i * (legendre(i-1, r) - legendre(i+1, r)) *
                 exp(-t * i * (i+1) / (4 * n)))
        if term != 0.0 and abs(term) < esp:
            return prob + term
        prob += term

    return prob


def hypergeo_fast(a, b, c, z, k=100):
    """Hypergeometric function"""
    terms = [1.0]
    for i in xrange(1, k+1):
        terms.append(float((i+a-1)*(i+b-1)*z)/(i+c-1)/i * terms[i-1])
    return sum(terms)


def hypergeo(a, b, c, z, k=100):
    """Hypergeometric function"""
    terms = [0.0]
    signs = [1.0]
    for i in xrange(1, k+1):
        term = float((i+a-1)*(i+b-1)*z)/(i+c-1)/i
        signs.append(util.sign(term) * signs[-1])
        if term == 0.0:
            break
        terms.append(log(abs(term)) + terms[i-1])
    #util.plot(terms)
    return sum(s*exp(i) for s, i in zip(signs, terms))


def loghypergeo(a, b, c, z, k=100):
    """Hypergeometric function"""
    terms = [0.0]
    signs = [1.0]
    for i in xrange(1, k+1):
        term = float((i+a-1)*(i+b-1)*z)/(i+c-1)/i
        signs.append(util.sign(term) * signs[-1])
        if term == 0.0:
            break
        terms.append(log(abs(term)) + terms[i-1])

    sgn = 1
    tot = -util.INF

    for s, t in zip(signs, terms):
        sgn, tot = stats.logadd_sign(sgn, tot, s, t)
    return sgn, tot


def hypergeo_mult(i, z1, z2, k=100):
    
     h1 = hypergeo(1-i, i+2, 2, z1, k)
     h2 = hypergeo(1-i, i+2, 2, z2, k)
     return h1 * h2
     


def freq_pdf_old(x, p, n, t, k=8):

    if x > 0.5:
        return freq_pdf2(1.0-x, 1.0-p, n, t, k)
    
    q = 1.0 - p
    prob = -util.INF
    sgn = 1
    t4n = t / (4*n)
    
    for i in xrange(1, k+1):
        #term = (p * q * i * (i+1) * (2*i+1) *
        #        hypergeo(1-i,i+2,2,p) * hypergeo(1-i,i+2,2,x) *
        #        exp(-t * i * (i+1) / (4*n)))

        lcoff = log(p * q * i * (i+1) * (2*i+1))
        h1 = hypergeo(1-i,i+2,2,p, i+2)
        h2 = hypergeo(1-i,i+2,2,x, i+2)
        sgn2 = util.sign(h1) * util.sign(h2)

        if sgn2 != 0:
            term = (lcoff + log(abs(h1)) + log(abs(h2)) +
                    (- i * (i+1) * t4n))
            sgn, prob = stats.logadd_sign(sgn, prob, sgn2, term)

    return sgn * exp(prob)


def freq_pdf(x, p, n, t, k=8):

    if x > 0.5:
        return freq_pdf(1.0-x, 1.0-p, n, t, k)
    
    q = 1.0 - p
    prob = -util.INF
    sgn = 1
    t4n = t / (4*n)
    
    for i in xrange(1, k+1):
        #term = (p * q * i * (i+1) * (2*i+1) *
        #        hypergeo(1-i,i+2,2,p) * hypergeo(1-i,i+2,2,x) *
        #        exp(-t * i * (i+1) / (4*n)))

        lcoff = log(p * q * i * (i+1) * (2*i+1))
        s1, h1 = loghypergeo(1-i,i+2,2,p, i+2)
        s2, h2 = loghypergeo(1-i,i+2,2,x, i+2)
        sgn2 = s1 * s2
        term = (lcoff + h1 + h2 - (i * (i+1) * t4n))
        
        sgn, prob = stats.logadd_sign(sgn, prob, sgn2, term)

    return sgn * exp(prob)


def freq_pdf2(x, p, n, t, k=8):
    r = 1 - 2*p
    z = 1 - 2*x

    prob = 0.0
    for i in xrange(1, k+1):
        term = ((2*i + 1) * (i - r*r) / float(i * (i+1)) *
                gegenbauer(i, r) * gegenbauer(i, z) *
                exp(-t * i * (i+1) / (4*n)))
        print term
        prob += term

    return prob


def freq_pdf3(x, p, n, t, k=8):
    q = 1.0 - p
    prob = 0.0
    for i in xrange(1, k+1):
        term = (p * q * i * (i+1) * (2*i+1) *
                hypergeo(1-i,i+2,2,p,40) * hypergeo(1-i,i+2,2,x,40) *
                exp(-t * i * (i+1) / (4*n)))
        prob += term

    return prob


def freq_pdf4(x, p, n, t, k=8):
    q = 1.0 - p
    prob = 0.0
    for i in xrange(1, k+1):
        term = (p * q * i * (i+1) * (2*i+1) *
                hypergeo_mult(i, p, x, 100) *
                exp(-t * i * (i+1) / (4*n)))
        prob += term

    return prob

    
    
if __name__ == "__main__":
    from rasmus.common import plotfunc

    if 0:
        for i in range(5):
            print "P_%d(x) = " % i, legendre_poly(i)
            print


    #========================
    # hypergeo speed

    a, b, c, z, k = 30, 20, 12, .3, 40
    
    util.tic("hypergeo_fast")
    for i in range(100):
        hypergeo_fast(a, b, c, z, k)
    util.toc()


    util.tic("hypergeo")
    for i in range(100):
        hypergeo(a, b, c, z, k)
    util.toc()
    
    util.tic("loghypergeo")
    for i in range(100):
        loghypergeo(a, b, c, z, k)
    util.toc()
    

if 0:
    p0 = .5
    k=30

    p = plotfunc(lambda x: freq_pdf(x, p0, 1000, 100, k=k),
                 .01, .99, .01, style="lines")
    p.plotfunc(lambda x: freq_pdf(x, p0, 1000, 200, k=k),
               .01, .99, .01, style="lines")
    p.plotfunc(lambda x: freq_pdf(x, p0, 1000, 500, k=k),
               .01, .99, .01, style="lines")
    p.plotfunc(lambda x: freq_pdf(x, p0, 1000, 1000, k=k),
               .01, .99, .01, style="lines")
    p.plotfunc(lambda x: freq_pdf(x, p0, 1000, 2000, k=k),
               .01, .99, .01, style="lines")
    p.plotfunc(lambda x: freq_pdf(x, p0, 1000, 3000, k=k),
               .01, .99, .01, style="lines")
    p.enableOutput(True)
    p.replot()

    #p.plotfunc(lambda x: normalPdf(x, (.5, .1135)),
    #           .01, .99, .01, style="lines")



if 0:
    p0 = .1

    p = plotfunc(lambda x: freq_pdf(x, p0, 1000, 100, k=25),
                 .01, .99, .01, style="lines")
    p.plotfunc(lambda x: freq_pdf(x, p0, 1000, 200, k=25),
               .01, .99, .01, style="lines")
    p.plotfunc(lambda x: freq_pdf(x, p0, 1000, 500, k=25),
               .01, .99, .01, style="lines")
    p.plotfunc(lambda x: freq_pdf(x, p0, 1000, 1000, k=25),
               .01, .99, .01, style="lines")
    p.plotfunc(lambda x: freq_pdf(x, p0, 1000, 2000, k=25),
               .01, .99, .01, style="lines")
    p.plotfunc(lambda x: freq_pdf(x, p0, 1000, 3000, k=25),
               .01, .99, .01, style="lines")
    p.enableOutput(True)
    p.replot()

    #p.plotfunc(lambda x: freq_pdf3(x, .5, 1000, 1000/10, k=40),
    #             .01, .99, .01, style="lines")


if 1:
    p0 = .5
    k=30

    p = plotfunc(lambda x: freq_pdf(x, p0, 1000, 30, k=k),
                 .01, .99, .01, style="lines")
    p.enableOutput(True)
    p.replot()


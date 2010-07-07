"""

   Coalescent methods


A note about population size.  In this code all population sizes N or n are
uncorrected.  If you need to compute a coalescent for a diploid species
you must multiply N by 2 before passing it to any of these functions.

"""

#=============================================================================
# imports

from __future__ import division

# python imports
import itertools
from itertools import chain, izip
from math import *
import random

# rasmus imports
from rasmus import treelib, stats, util
from rasmus.symbolic import *

# compbio imports
from . import birthdeath

#=============================================================================
# single coalescent PDFs, CDFs, and sampling functions


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
    Returns a sample coalescent time for 'k' individuals in a population 'n'
    """

    # k choose 2
    k2 = k * (k-1) / 2
    k2n = k2 / n
    return random.expovariate(k2n)


def sample_coal_times(k, n):
    """
    Returns a sampling of (k-1) coalescences for 'k' lineages in a
    population of size 'n'.
    """
    times = [0]
    for j in xrange(k, 1, -1):
        times.append(times[-1] + sample_coal(j, n))
    return times[1:]


def prob_mrca(t, k, n):
    """
    Probability density function of the age 't' of the most recent
    common ancestor (MRCA) of 'k' lineages in a population size 'n'
    """

    s = 0.0
    for i in xrange(1, k):
        lam = (i+1) * i / 2.0 / n
        s += lam * exp(- lam * t) * mrca_const(i, 1, k-1)
    return s


def cdf_mrca(t, k, n):
    """
    Cumulative probability density of the age 't' of the most recent common
    ancestor (MRCA) of 'k' lineages in a population size 'n'
    """
    
    if k == 1:
        return 1.0
    
    s = 0.0
    for i in xrange(1, k+1):
        lam = i * (i-1) / (2.0 * n)
        p = 1.0
        for y in xrange(1, i):
            p *= (y-k) / (k+y)
        s += exp(-lam * t) * (2*i - 1) * p
    return s


def mrca_const(i, a, b):
    """A constant used in calculating MRCA"""

    # i+1 choose 2
    y = (i+1) * i / 2.0
    prod = 1.0

    for j in xrange(a, b+1):
        if j == i:
            continue
        # j+1 choose 2
        x = (j+1) * j / 2.0
        prod *= x / (x - y)
    return prod



def prob_coal_bounded(t, k, n, T):
    """
    Probability density function of seeing a coalescence at 't' from
    'k' lineages in a population of size 'n' with bounding time 'T'
    """

    if t > T:
        return 0.0
    
    if k == 2:
        prob_coal(t, k, n)
    return prob_coal(t, k, n) * cdf_mrca(T-t, k-1, n) / \
           cdf_mrca(T, k, n)


def cdf_coal_bounded(t, k, n, T):
    """
    Cumalative density function of seeing a coalescence at 't' from
    'k' lineages in a population of size 'n' with bounding time 'T'
    """
    i = k - 1

    lam_i = (i+1)*i/2.0 / n
    C = [mrca_const(j, 1, i-1) for j in xrange(1, i)]
    A = lam_i / n / cdf_mrca(T, k, n)
    B = sum(C) / lam_i
    F = [C[j-1] * exp(-(j+1)*j/2.0/n * T) / ((j+1)*j/2.0/n - lam_i)
         for j in xrange(1, i)]
    
    return (lam_i / cdf_mrca(T, k, n) *
            (B * (1-exp(-lam_i * t))
             - sum(F[j-1] * (exp(((j+1)*j/2.0/n - lam_i)*t)-1)
                   for j in xrange(1, i))))

    

def sample_coal_bounded(k, n, T):
    """
    Sample a coalescent time 't' for 'k' lineages and population 'n'
    on the condition that the MRCA is before 'T'
    """

    # special case
    if k == 2:
        return sample_coal_bounded2(n, T)

    # this code solves this equation for t
    #   cdf(t) - p = 0
    # where p ~ U(0, 1)

    import scipy.optimize

    i = k - 1
    p = random.random()

    # compute constants
    lam_i = (i+1)*i/2.0 / n
    C = [mrca_const(j, 1, i-1) for j in xrange(1, i)]
    A = lam_i / cdf_mrca(T, k, n)
    B = sum(C) / lam_i
    F = [C[j-1] * exp(-(j+1)*j/2.0/n * T) / ((j+1)*j/2.0/n - lam_i)
         for j in xrange(1, i)]

    # CDF(t) - p
    def f(t):
        if t <= 0:
            return t - p
        if t >= T:
            return 1.0 - p + (t - T)

        return (A * (B * (1-exp(-lam_i * t))
             - sum(F[j-1] * (exp(((j+1)*j/2.0/n - lam_i)*t)-1)
                   for j in xrange(1, i)))) - p
    
    return scipy.optimize.brentq(f, 0.0, T, disp=False)



def sample_coal_bounded2(n, T):
    """
    Sample a coalescent time 't' for 'k=2' lineages and population 'n'
    on the condition that the MRCA is before 'T'
    """

    # sample from a truncated expontial distribution

    # k choose 2
    lam = 1 / n
    p = exp(-lam * T)
    return - log(random.uniform(p, 1.0)) / lam


def sample_coal_bounded_reject(k, n, T):
    """
    Sample a coalescent time 't' for 'k' lineages and population 'n'
    on the condition that the MRCA is before 'T'

    Uses rejection sampling.  It works but is very inefficient.
    """

    i = k - 1
    consts = [mrca_const(j, 1, i-1) for j in xrange(1, i)]
    x = sum(consts)

    while True:
        while True:
            t = sample_coal(k, n)
            if t < T:
                break

        if i == 1:
            return t

        y = sum(mrca_const(j, 1, i-1) * exp(-((j+1) * j / 2.0 / n) * (T - t))
                for j in xrange(1, i))

        r = 1 - y / x

        if random.random() < r:
            return t


def prob_coal_counts(a, b, t, n):
    """
    The probabiluty of going from 'a' lineages to 'b' lineages in time 't'
    with population size 'n'
    """
    
    s = 0.0
    for k in xrange(b, a+1):
        i = exp(-k*(k-1)*t/2.0/n) * \
            (2*k-1)*(-1)**(k-b) / stats.factorial(b) / \
            stats.factorial(k-b) / (k+b-1) * \
            stats.prod((b+y)*(a-y)/(a+y) for y in xrange(k))
        s += i
    return s

def count_lineages_per_branch(tree, recon, stree, rev_recon=None):
    """
    Returns the count of gene lineages present at each node in the species
    tree 'tree' given a gene tree 'tree' and reconciliation 'recon'
    """

    
    # init reverse reconciliation
    if rev_recon is None:
        rev_recon = {}
        nodes = set(tree.postorder())
        for node, snode in recon.iteritems():
            if node not in nodes:
                raise Exception("node '%s' not in tree" % node.name)
            rev_recon.setdefault(snode, []).append(node)

    # init lineage counts
    lineages = {}
    for snode in stree:
        if snode.is_leaf():
            lineages[snode] = [len([x for x in rev_recon[snode]
                                   if x.is_leaf()]), 0]
        else:
            lineages[snode] = [0, 0]
    
    # iterate through species tree branches
    for snode in stree.postorder():
        if snode.parent:
            # non root branch
            a = lineages[snode][0]

            # subtract number of coals in branch
            b = a - len([x for x in rev_recon.get(snode, [])
                         if not x.is_leaf()])
            lineages[snode][1] = b
            lineages[snode.parent][0] += b
        else:
            lineages[snode][1] = 1

    return lineages



def prob_coal_recon_topology(tree, recon, stree, n):
    """
    Returns the log probability of a reconciled gene tree ('tree', 'recon')
    from the coalescent model given a species tree 'stree' and
    population sizes 'n'
    """
    
    popsizes = init_popsizes(stree, n)

    lineages = count_lineages_per_branch(tree, recon, stree)

    # log probability
    lnp = 0.0


    # iterate through species tree branches
    for snode in stree.postorder():
        if snode.parent:
            # non root branch
            a, b = lineages[snode]
            
            lnp += log(prob_coal_counts(a, b, snode.dist,
                                        popsizes[snode.name]))
            lnp -= log(num_labeled_histories(a, b))            
        else:
            a = lineages[snode][0]
            lnp -= log(num_labeled_histories(a, 1))

    
    # correct for topologies H(T)
    # find connected subtrees that are in the same species branch
    subtrees = []
    subtree_root = {}
    for node in tree.preorder():
        if node.parent and recon[node] == recon[node.parent]:
            subtree_root[node] = subtree_root[node.parent]
        else:
            subtrees.append(node)
            subtree_root[node] = node

    # find leaves through recursion
    def walk(node, subtree, leaves):
        if node.is_leaf():
            leaves.append(node)
        elif (subtree_root[node.children[0]] != subtree and
              subtree_root[node.children[1]] != subtree):
            leaves.append(node)
        else:
            for child in node.children:
                walk(child, subtree, leaves)

    # apply correction for each subtree
    for subtree in subtrees:
        leaves = []
        for child in subtree.children:
            walk(subtree, subtree, leaves)
        if len(leaves) > 2:
            lnp += log(birthdeath.num_topology_histories(subtree, leaves))

    return lnp



def prob_coal_recon_topology2(tree, recon, stree, n):
    """
    Returns the log probability of a reconciled gene tree ('tree', 'recon')
    from the coalescent model given a species tree 'stree' and
    population sizes 'n'
    """
    
    popsizes = init_popsizes(stree, n)

    # log probability
    lnp = 0.0 

    # init reverse reconciliation
    rev_recon = {}
    nodes = set(tree.postorder())
    for node, snode in recon.iteritems():
        if node not in nodes:
            raise Exception("node '%s' not in tree" % node.name)
        rev_recon.setdefault(snode, []).append(node)

    # init lineage counts
    lineages = {}
    for snode in stree:
        if snode.is_leaf():
            lineages[snode] = len([x for x in rev_recon[snode]
                                   if x.is_leaf()])
        else:
            lineages[snode] = 0

    # iterate through species tree branches
    for snode in stree.postorder():
        if snode.parent:
            # non root branch
            u = lineages[snode]

            # subtract number of coals in branch
            v = u - len([x for x in rev_recon.get(snode, [])
                         if not x.is_leaf()])            
            lineages[snode.parent] += v

            lnp += log(prob_coal_counts(u, v, snode.dist,
                                        popsizes[snode.name]))
            lnp -= log(num_labeled_histories(u, v))            
        else:
            u = lineages[snode]
            lnp -= log(num_labeled_histories(u, 1))

    
    # correct for topologies H(T)
    # find connected subtrees that are in the same species branch
    subtrees = []
    subtree_root = {}
    for node in tree.preorder():
        if node.parent and recon[node] == recon[node.parent]:
            subtree_root[node] = subtree_root[node.parent]
        else:
            subtrees.append(node)
            subtree_root[node] = node

    # find leaves through recursion
    def walk(node, subtree, leaves):
        if node.is_leaf():
            leaves.append(node)
        elif (subtree_root[node.children[0]] != subtree and
              subtree_root[node.children[1]] != subtree):
            leaves.append(node)
        else:
            for child in node.children:
                walk(child, subtree, leaves)

    # apply correction for each subtree
    for subtree in subtrees:
        leaves = []
        for child in subtree.children:
            walk(subtree, subtree, leaves)
        if len(leaves) > 2:
            lnp += log(birthdeath.num_topology_histories(subtree, leaves))

    return lnp


def cdf_mrca_bounded_tree(gene_counts, T, stree, sroot, n,
                          tree=None, recon=None):
    """
    What is the log probability that multispecies coalescent in species
    tree 'stree' with population sizes 'n' and extant gene counts 'gene_counts'
    will have a MRCA that occurs in branch 'sroot' before time 'T'.

    As a convenience, you can pass None for gene_counts and give a reconciled
    gene tree instead ('tree', 'recon').
    """

    if sroot is None:
        sroot = stree.root

    # init gene counts
    if gene_counts is None:
        gene_counts = dict.fromkeys(sroot.leaf_names(), 0)
        for leaf in tree.leaves():
            gene_counts[recon[leaf.name]] += 1

    popsizes = init_popsizes(stree, n)

    # get time to MRCA above sroot
    stimes = treelib.get_tree_timestamps(stree, sroot)
    root_time = T - stimes[sroot]

    # use dynamic programming to calc prob of lineage counts
    prob_counts = {}
    def walk(node):
        if node.is_leaf():
            M = gene_counts[node.name]
            prob_counts[node] = [0.0] * (M+1)
            prob_counts[node][M] = 1.0
            return M
        
        else:
            assert len(node.children) == 2
            c1 = node.children[0]
            c2 = node.children[1]
            t1 = c1.dist
            t2 = c2.dist
            M1 = walk(c1)
            M2 = walk(c2)
            M = M1 + M2 # max lineage counts in this snode
            n1 = popsizes[c1.name]
            n2 = popsizes[c2.name]

            prob_counts[node] = [0, 0]
            for k in xrange(2, M+1):
                prob_counts[node].append(sum(
                    sum(prob_coal_counts(i, m, t1, n1) *
                        prob_counts[c1][i]
                        for i in xrange(m, M1+1)) * 
                    sum(prob_coal_counts(i, k-m, t2, n2) *
                        prob_counts[c2][i]
                        for i in xrange(k-m, M2+1))
                    for m in xrange(1, k)))
                #print prob_counts[node]

            assert abs(sum(prob_counts[node]) - 1.0) < .001
            return M
    M = walk(sroot)

    #util.print_dict(prob_counts, key=lambda x: x.name)
    
    # count final sum
    p = sum(cdf_mrca(root_time, k, popsizes[sroot.name]) *
            prob_counts[sroot][k]
            for k in xrange(2, M+1))

    return log(p)


def num_labeled_histories(nleaves, nroots):
    n = 1.0
    for i in xrange(nroots + 1, nleaves+1):
        n *= i * (i - 1) / 2.0
    return n


def prob_coal_bounded_recon_topology(tree, recon, stree, n, T):
    """
    Returns the log probability of a reconciled gene tree ('tree', 'recon')
    from the coalescent model given a species tree 'stree' and
    population sizes 'n' and stopping time 'T'
    """

    popsizes = init_popsizes(stree, n)

    p = prob_coal_recon_topology(tree, recon, stree, n)
    times = treelib.get_tree_timestamps(tree)
    lineages = count_lineages_per_branch(tree, recon, stree)
    k_root = lineages[tree.root][0]
    T_root = T - times[tree.root]
    return log(cdf_mrca(T_root, k_root, popsizes[recon[tree.root].name])) + p \
           - cdf_mrca_bounded_tree(None, T, stree, stree.root, n,
                                   tree=tree, recon=recon)
    


#=============================================================================
# sampling coalescent trees
#
#  - normal kingman coalescent
#  - censored coalescent
#  - bounded coalescent (conditioned on completion before a fixed time)
#


def sample_coal_tree(k, n):
    """
    Returns a simulated coalescent tree for 'k' leaves from a population 'n'.
    """
    times = [0]
    for j in xrange(k, 1, -1):
        times.append(times[-1] + sample_coal(j, n))
    return make_tree_from_times(times)[0]


def sample_coal_tree_bounded(k, n, T, capped=False):
    """
    Returns a simulated coalescent tree for 'k' leaves from a populations 'n'
    with fixed maximum time 't'.  The simulation is conditioned on returning
    a tree that completely coaleces before time 'T'.

    capped -- if True an artificial root to the tree.  Used primarily by
              other methods.
    """    
    times = [0]
    for j in xrange(k, 1, -1):
        times.append(times[-1] + sample_coal_bounded(j, n, T - times[-1]))
    return make_tree_from_times(times, t=T, capped=capped)[0]


def sample_coal_tree_bounded_reject(k, n, T, capped=False):
    """
    Returns a simulated coalescence tree for k leaves from a populations n
    with fixed maximum time t.  The simulation is conditioned on returning
    a tree that completely coaleces before time T.

    This works, but is very inefficient.  Use sample_coal_tree_bounded
    instead.
    """

    # sample times with rejection sampling
    while True:
        times = [0]
        for j in xrange(k, 1, -1):
            times.append(times[-1] + sample_coal(j, n))
        if times[-1] < t:
            break

    return make_tree_from_times(times, t=T, capped=capped)[0]


def sample_coal_tree_fixed(k, n, t, capped=False):
    """
    Returns a simulated coalescence tree for 'k' leaves from a population size
    'n' with a fixed maximum time 't'.

    The return value is the tuple (tree, lineages) where lineages is a set
    of lineages that have not yet coalesced.

    capped -- if True, remaining lineages are added as children to a artificial
              tree root.
    """

    times = [0]
    for j in xrange(k, 1, -1):
        times.append(times[-1] + sample_coal(j, n))
        if times[-1] > t:
            times.pop()
            break
    
    return make_tree_from_times(times, k, t, capped=capped)


def init_popsizes(stree, n):
    """
    Uses 'n' to initialize a population size dict for species tree 'stree'
    """

    if isinstance(n, (int, float)):
        return dict.fromkeys(stree.nodes.keys(), n)
    elif isinstance(n, dict):
        return n
    else:
        raise Exception("n must be a int or dict.")



def sample_multicoal_tree(stree, n, leaf_counts=None,
                          namefunc=None):
    """
    Returns a gene tree from a multi-species coalescence process

    stree       -- species tree
    n           -- population size (int or dict)
                   If n is a dict it must map from species name to
                   population size.
    leaf_counts -- dict of species names to a starting gene count.
                   Default is 1 gene per extant species.
    namefunc    -- a function that generates new gene names given a species
                   name.
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
    popsizes = init_popsizes(stree, n)

    # init gene counts
    counts = dict((n.name, 0) for n in stree)
    counts.update(leaf_counts)

    # init reconciliation
    recon = {}

    # subtrees
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




def make_tree_from_times(times, k=None, t=None, leaves=None, capped=False):
    """
    Returns a Tree from a list of divergence times.

    The topology is choosen by randomly choosing pairs of leaves.
    """

    # initialize k
    if k is None:
        if leaves is not None:
            k = len(leaves)
        else:
            k = len(times)
            
    tree = treelib.Tree()

    # initialize k children
    if leaves is None:
        children = set(treelib.TreeNode(tree.new_name()) for i in xrange(k))
    else:
        children = set(treelib.TreeNode(name) for name in leaves)
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
            if t is not None:
                node.dist = t - node.data["time"]
            else:
                node.dist = 0.0
        else:
            node.dist = node.parent.data["time"] - node.data["time"]

    # for convenience cap the tree for easy drawing/manipulation
    if capped:
        tree.make_root()
        for node in children:
            tree.add_child(tree.root, node)
    else:
        # set root
        if len(children) == 1:
            tree.root = list(children)[0]
    
    # return tree and remaining lineages
    return tree, children
    



#=============================================================================
# allele frequency

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
    return ((i * (i+1)) / 2.0 * hypergeo(i+2, 1 - i, 2, (1 - r) / 2.0))


def gegenbauer2(i, r):
    return ((i * (i+1)) / float((2*i+1)*(1-r*r)) *
            (legendre(i-1, r) - legendre(i+1, r)))

def gegenbauer3(n, a, z):

    tot = 0
    for k in xrange(int(n/2)+1):
        tot += ((-1)**k * stats.gamma(n - k + a) / (
                stats.gamma(a) * stats.factorial(k) * stats.factorial(n - 2*k))
                * ((2*z) ** (n - 2*k)))
    return tot




def prob_fix(p, n, t, k=8, esp=0.001):
    """Probability of fixation"""
    r = 1 - 2*p
    
    prob = p
    for i in xrange(1, k+1):
        term = (.5 * (-1)**i * (legendre(i-1, r) - legendre(i+1, r)) *
                 exp(-t * i * (i+1) / (4 * n)))
        if term != 0.0 and abs(term) < esp:
            return prob + term
        prob += term

    return prob


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
    return sum(s*exp(i) for s, i in zip(signs, terms))


def loghypergeo(a, b, c, z, k=100):
    """
    Hypergeometric function

    Performs computation in log-space
    """
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




#=============================================================================
    
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


    if 0:
        p0 = .5
        k=30

        p = plotfunc(lambda x: freq_pdf(x, p0, 1000, 30, k=k),
                     .01, .99, .01, style="lines")
        p.enableOutput(True)
        p.replot()




#=============================================================================
# old versions


def hypergeo_old(a, b, c, z, k=100):
    """Hypergeometric function"""
    terms = [1.0]
    for i in xrange(1, k+1):
        terms.append(float((i+a-1)*(i+b-1)*z)/(i+c-1)/i * terms[i-1])
    return sum(terms)


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


def cdf_mrca2(t, k, n):
    """
    Cumulative probability density of the age 't' of the most recent common
    ancestor (MRCA) of 'k' lineages in a population size 'n'
    """

    if k == 1:
        return 1.0
    
    s = 0.0
    for i in xrange(1, k):
        lam = (i+1) * i / 2.0 / n
        s += (1 - exp(- lam * t)) * mrca_const(i, 1, k-1)
    return s



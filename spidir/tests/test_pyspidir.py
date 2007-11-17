#!/usr/bin/env python

import gc
import sys


import pyspidir


#gc.set_debug(gc.DEBUG_LEAK)

#gc.collect()

#print sys.getrefcount(a)

ptree = [3, 3, 4, 4, -1]
seqs = ["AACC", "ATCC", "AAGC"]
bgfreq = [.25, .25, .25, .25]
ratio = .5
maxiter = 10

for i in xrange(100):
    dists, logl = pyspidir.mlhkydist(ptree, seqs, bgfreq, ratio, maxiter)
    print dists, logl
    print "ptree", sys.getrefcount(ptree)
    print "dists", sys.getrefcount(dists)
    print "logl", sys.getrefcount(logl)

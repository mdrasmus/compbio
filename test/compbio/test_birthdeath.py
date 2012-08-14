
import unittest
from math import *

from rasmus.common import *
from rasmus import stats
from rasmus.testing import *

from compbio import birthdeath



#=============================================================================
# test coalescence times (normal, censored, bounded)

def rannala1996_prob_birth_death1(n, t, birth, death):

    ert = exp(-(birth - death)*t)
    p0 = (death - death * ert) / (birth - death * ert)
    p1 = (birth - death)**2 * ert / (birth - death * ert)**2

    if n == 0:
        return p0
    elif n == 1:
        return p1
    else:
        return (birth/death)**n * p1 * p0**(n-1)


def rannala1996_prob_birth_death1_fix(n, t, birth, death):

    ert = exp(-(birth - death)*t)
    p0 = (death - death * ert) / (birth - death * ert)
    p1 = (birth - death)**2 * ert / (birth - death * ert)**2

    if n == 0:
        return p0
    elif n == 1:
        return p1
    else:
        return p1 * (birth/death * p0)**(n-1)

    
class BD (unittest.TestCase):

    def test_prob_birth_death1(self):

        t = 1.0
        birth = 0.5
        death = 0.2

        counts = [birthdeath.sample_birth_death_count(1, t, birth, death)
                  for i in xrange(10000)]

        p = plotdistrib(counts, width=1, low=-.5)
        p.plot([birthdeath.prob_birth_death1(i, t, birth, death)
                for i in xrange(0, max(counts))])
        pause()


    def test_rannala1996_prob_birth_death1(self):

        t = 1.0
        birth = 0.5
        death = 0.2

        counts = [birthdeath.sample_birth_death_count(1, t, birth, death)
                  for i in xrange(10000)]

        p = plotdistrib(counts, width=1, low=-.5)
        p.plot([rannala1996_prob_birth_death1(i, t, birth, death)
                for i in xrange(0, max(counts))])
        pause()


        p = plotdistrib(counts, width=1, low=-.5)
        p.plot([rannala1996_prob_birth_death1_fix(i, t, birth, death)
                for i in xrange(0, max(counts))])
        pause()

        #print 


    def test_prob_birth_death1_eq(self):

        t = 1.0
        n = 5
        birth = 2.0
        death = 2.0

        for i in frange(.6, .99, .05):
            print birthdeath.prob_birth_death1(n, t, birth, death*i)
        print birthdeath.prob_birth_death1(n, t, birth, death*.999)
        v = birthdeath.prob_birth_death1(n, t, birth, death)
        print "==", v

    def test_prob_birth_death1_eq_null(self):

        t = 1.0
        n = 0
        birth = 2.0
        death = 2.0

        for i in frange(.6, .99, .05):
            print birthdeath.prob_birth_death1(n, t, birth, death*i)
        print birthdeath.prob_birth_death1(n, t, birth, death*.999)
        v = birthdeath.prob_birth_death1(n, t, birth, death)
        print "==", v

    def test_birth_wait_time_eq(self):

        t = 0.5
        n = 2
        T = 1.0
        birth = 2.0
        death = 2.0

        for i in frange(.6, .99, .05):
            print birthdeath.birth_wait_time(t, n, T, birth, death*i)
        print birthdeath.birth_wait_time(t, n, T, birth, death*.999)
        v = birthdeath.birth_wait_time(t, n, T, birth, death)
        print "==", v


    def test_prob_no_birth_eq(self):

        n = 2
        T = 1.2
        birth = 2.0
        death = 2.0

        for i in frange(.6, .99, .05):
            print birthdeath.prob_no_birth(n, T, birth, death*i)
        print birthdeath.prob_no_birth(n, T, birth, death*.999)
        v = birthdeath.prob_no_birth(n, T, birth, death)
        print "==", v


    def test_sample_birth_wait_time_eq(self):

        n = 2
        T = 1.2
        birth = 2.0
        death = 2.0

        print "samples"
        for i in xrange(10):
            print birthdeath.sample_birth_wait_time(n, T, birth, death)


#=============================================================================


if __name__ == "__main__":
    test_main()

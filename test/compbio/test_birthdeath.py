
from math import exp
import unittest

from compbio import birthdeath

from rasmus.testing import eq_sample_pmf

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
        """Sampling and PDF for birth-death from single lineage."""
        t = 1.0
        birth = 0.5
        death = 0.2

        counts = [birthdeath.sample_birth_death_count(1, t, birth, death)
                  for i in xrange(10000)]
        eq_sample_pmf(
            counts,
            lambda i: birthdeath.prob_birth_death1(i, t, birth, death),
            pval=0.01)

    def test_rannala1996_prob_birth_death1(self):

        t = 1.0
        birth = 0.5
        death = 0.2

        counts = [birthdeath.sample_birth_death_count(1, t, birth, death)
                  for i in xrange(10000)]

        # original equation should fail
        try:
            eq_sample_pmf(
                counts,
                lambda i: rannala1996_prob_birth_death1(i, t, birth, death))
        except:
            pass
        else:
            raise AssertionError

        eq_sample_pmf(
            counts,
            lambda i: rannala1996_prob_birth_death1_fix(i, t, birth, death))

    def test_prob_birth_death1_eq(self):
        """
        Sampling and PDF for birth-death from single lineage birth=death rate.
        """
        t = 1.0
        birth = 0.5
        death = 0.5

        counts = [birthdeath.sample_birth_death_count(1, t, birth, death)
                  for i in xrange(10000)]
        eq_sample_pmf(
            counts,
            lambda i: birthdeath.prob_birth_death1(i, t, birth, death))

    def test_birth_wait_time_eq(self):

        t = 0.5
        n = 2
        T = 1.0
        birth = 2.0
        death = 2.0
        self.assertAlmostEqual(
            birthdeath.birth_wait_time(t, n, T, birth, death*.9999),
            birthdeath.birth_wait_time(t, n, T, birth, death),
            places=4)

    def test_prob_no_birth_eq(self):

        n = 2
        T = 1.2
        birth = 2.0
        death = 2.0
        self.assertAlmostEqual(
            birthdeath.prob_no_birth(n, T, birth, death*.9999),
            birthdeath.prob_no_birth(n, T, birth, death),
            places=4)

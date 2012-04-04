
from itertools import izip
from . import util
from . import stats


def fequal(f1, f2, rel=.0001, eabs=1e-12):
    """assert whether two floats are approximately equal"""
    
    if f1 == f2:
        return
    
    if f2 == 0:
        err = f1
    elif f1 == 0:
        err = f2
    else:
        err = abs(f1 - f2) / abs(f2)
    x = (err < rel)

    if abs(f1 - f2) < eabs:
        return

    assert x, "%e != %e  [rel=%f, abs=%f]" % (f1, f2, err, abs(f1 - f2))


def fequals(f1, f2, rel=.0001, eabs=1e-12):
    for i, j in izip(f1, f2):
        fequal(i, j, rel=rel, eabs=eabs)


def integrate(func, a, b, step):
    return sum(func(i) * step for i in util.frange(a, b, step))


def eq_sample_pdf(samples, pdf,
                  ndivs=20, start=-util.INF, end=util.INF, pval=.05,
                  step=None):
    """Returns true if a sample matches a distribution"""

    if step is None:
        step = (max(samples) - min(samples)) / float(ndivs)

    cdf = lambda x, params: integrate(pdf, x, x+step, step/10.0)

    chi2, p = stats.chi_square_fit(cdf, [], samples,
                                   ndivs=ndivs, start=start, end=end)
    
    assert p >= pval, p


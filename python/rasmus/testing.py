
import sys, os, shutil, unittest
import optparse
from itertools import izip
from . import util
from . import stats

#=============================================================================
# common utility functions for testing


def clean_dir(path):
    if os.path.exists(path):
        shutil.rmtree(path)

def makedirs(path):
    if not os.path.exists(path):
        os.makedirs(path)


def make_clean_dir(path):
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)




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



#=============================================================================
# common unittest functions


class OptionParser (optparse.OptionParser):
    def __init__(self, *args):
        optparse.OptionParser.__init__(self, *args)

    def exit(self, code, text):
        pass


def list_tests(stack=0):
    
    # get environment
    var = __import__("__main__").__dict__

    for name, obj in var.iteritems():
        if isinstance(obj, type) and issubclass(obj, unittest.TestCase):
            for attr in dir(obj):
                if attr.startswith("test"):
                    print "%s.%s" % (name, attr),
                    doc = getattr(obj, attr).__doc__
                    if doc:
                        print "--", doc.split("\n")[0]
                    else:
                        print


def test_main():

    o = OptionParser()
    o.add_option("-l", "--list_tests", action="store_true")

    conf, args = o.parse_args(sys.argv)

    if conf.list_tests:
        list_tests(1)
        return

    unittest.main()


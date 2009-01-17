
from rasmus.common import *

from spidir import *

def birthDeath(ngenes, time, birthRate, deathRate):
    l = birthRate
    u = deathRate
    r = l - u
    a = u / l

    ut = (1.0 - exp(-r*time)) / (1.0 - a * exp(-r*time))
    p0 = a*ut
    
    if ngenes == 0:
        return p0
    
    # (1.0 - p0)*(1.0 - ut) * ut^{ngenes-1}
    return (1.0 - p0)*(1.0 - ut) * (ut**ngenes-1)


if 0:
    s = 0
    t = 1.0
    l = 3
    u = .5

    plotfunc(lambda s: birthDeathCount(s, t, l, u), 0, 10, 1)


if 1:
    for i in range(1, 9):
        n = numHistories(i)
        n2 = int(factorial(i) * factorial(i-1) / 2**(i-1))
        print i, n, n2
        assert n == n2
        

if 1:
    ctree = spidir.tree2ctree(tree)

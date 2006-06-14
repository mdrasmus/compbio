import util, stats
import random

def metropolis(proposal, density, initState, steps):
    states = [initState]
    p = density(initState)
    for i in xrange(steps):
        while True:
            state = proposal(states[-1])
            q = density(state)
            if q >= p or p == 0:
                break
            elif random.random() < q / p:
                break
        p = q
        states.append(state)
    return states
    
    
def randomWalk(state, sigma):
    return state + random.gauss(0, sigma)

def walk(transition, initState, steps):
    states = [initState]
    for i in xrange(steps):
        states.append(transition(states[-1]))
    return states





# testing
if __name__ == "__main__":
    # gauss_sample = map(lambda x: random.gauss(0, 1), range(10000))
    # d, p = util.plotdistrib(gauss_sample)
    
    func = lambda x: (stats.normal(x, 5, 1) + \
                      stats.normal(x, -2, 1) + \
                      stats.normal(x, 2, 1)) * .33333
    #func = lambda x: stats.poisson(x, .1)
    #func = lambda x: stats.gamma(x, 20, 1)
    #func = lambda x: stats.bionomial(10, x, .6)
    
    if True:
        p = util.plotfunc(func, -10, 10, .1, plab="distrib")
        states = metropolis(lambda x: randomWalk(x, 1), func, 0, 10000)
        d = util.distrib(states, 100)
        p.plot(d[0], d[1], style="boxes", plab="sampling")

    x = walk(lambda x: randomWalk(x, 1), 0, 10000)
    y = walk(lambda x: randomWalk(x, 1), 0, 10000)
    p2 = util.plot(x,y,style="lines")

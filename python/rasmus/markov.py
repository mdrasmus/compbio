from rasmus import util, stats

import random

def metropolis(proposal, density, initState, steps):
    state = initState
    p = density(initState)
    
    for i in xrange(steps):
        nextState = proposal(state)
        p2 = density(nextState)
        if p2 >= p or random.random() < p2 / p:
            # accept new state
            p = p2
            state = nextState
        yield state

        
def metropolisHastings(proposal, probProposal, density, initState, steps):
    state = initState
    p = density(initState)
    
    for i in xrange(steps):
        nextState = proposal(state)
        p2 = density(nextState)
        q2 = probProposal(nextState, state)
        q = probProposal(state, nextState)

        if p * q2 == 0.0:
            accept = 1.0
        else:
            accept = min(1.0, (p2 * q) / (p * q2))
        
        if accept == 1.0 or random.random() <= accept:
            # accept new state
            p = p2
            state = nextState
        yield state
 
    
def randomWalkNormal(sigma):
    def func(state):
        return state + random.gauss(0, sigma)
    return func

def randomWalkGamma(alpha, beta):
    def func(state):
        return state + random.gammavariate(alpha, 1.0/beta) - alpha / beta
    return func

def condProbGamma(alpha, beta):
    def func(nextState, state):
        return gammaPdf(nextState - state + alpha/beta, [alpha, beta])
    return func

def walk(transition, initState, steps):
    states = [initState]
    for i in xrange(steps):
        states.append(transition(states[-1]))
    return states




# testing
if __name__ == "__main__":
    # gauss_sample = map(lambda x: random.gauss(0, 1), range(10000))
    # d, p = util.plotdistrib(gauss_sample)
    
    func = lambda x: (stats.normalPdf(x, [5, 1]) + \
                      stats.normalPdf(x, [-2, 1]) + \
                      stats.normalPdf(x, [2, 1])) * .33333
    #func = lambda x: stats.poisson(x, .1)
    #func = lambda x: stats.gamma(x, 20, 1)
    #func = lambda x: stats.binomial(10, x, .6)
    
    if True:
        p = util.plotfunc(func, -10, 10, .1, plab="distrib")
        states = metropolis(randomWalkNormal(1.0), func, 0, 10000)
        d = util.distrib(states, 100)
        p.plot(d[0], d[1], style="boxes", plab="sampling")

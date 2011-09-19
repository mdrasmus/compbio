"""

  Light weight generic HMM algorithms


  Methods need the following user-defined functions:

    get_num_states(pos)
    prob_prior(state)
    prob_emission(pos, state)
    prob_transition(pos1, state1, pos2, state2)


    probs = [[j for j in get_num_states(i)] 
             for i in xrange(npositions)]


"""

import random

from rasmus import util, stats
from stats import logadd


class HMM (object):
    def __init__(self):
        pass

    def set_callbacks(self, get_num_states=None,
                      prob_prior=None,
                      prob_emission=None,
                      prob_transition=None,
                      emit=None):
        if get_num_states:
            self.get_num_states = get_num_states
        if prob_prior:
            self.prob_prior = prob_prior
        if prob_emission:
            self.prob_emission = prob_emission
        if prob_transition:
            self.prob_transition = prob_transition
        if emit:
            self.emit = emit


    def get_num_states(self, pos):
        return 0

    def prob_prior(self, state):
        return 0.0

    def prob_emission(self, pos, state):
        return 0.0

    def prob_transition(self, pos1, pos2, state2):
        return 0.0

    def emit(self, state):
        return None



def sample_hmm_first_state(model):
    state = 0
    nstates = model.get_num_states(0)
    p = model.prob_prior(state)
    pick = log(random.random())
    while pick > p and state < nstates:
        state += 1
        p = logadd(p, model.prob_prior(state))
    return state


def sample_hmm_next_state(model, pos, state):
    nstates = model.get_num_states(pos)
    state2 = 0
    p = model.prob_transition(pos-1, state, pos, state2)
    pick = log(random.random())
    while pick > p and state2 < nstates:
        state2 += 1
        p = logadd(p, model.prob_transition(pos-1, state, pos, state2))
    return state2
        


def sample_hmm_states(model):
    
    # sample first state
    pos = 0
    state = sample_hmm_first_state(model)
    yield state

    # sample next states
    pos = 1
    while True:
        state = sample_hmm_next_state(model, pos, state)
        yield state
        pos += 1


def sample_hmm_data(model, states=None):
    
    if states is None:
        states = sample_hmm_states(model)

    for state in states:
        yield model.emit(state)



def viterbi(model, n, verbose=False):

    probs = []
    ptrs = []

    # calc first position
    nstates = model.get_num_states(0)
    probs.append([model.prob_prior(j) + model.prob_emission(0, j)
                  for j in xrange(nstates)])
    ptrs.append([-1] * nstates)
    
    if n > 20:
        step = (n // 20)
    else:
        step = 1
    
    # loop through positions
    for i in xrange(1, n):
        if verbose and i % step == 0:
            print " viterbi iter=%d/%d, lnl=%f" % (i+1, n, max(probs[-1]))

        nstates1 = model.get_num_states(i-1)
        nstates2 = model.get_num_states(i)
        col1 = probs[i-1]

        # find max transition and emission
        col2 = []
        col2_ptr = []
        for k in xrange(nstates2):
            top = -util.INF
            ptr = -1
            emit = model.prob_emission(i, k)
            for j in xrange(nstates1):
                p = col1[j] + model.prob_transition(i-1, j, i, k) + emit
                if p > top:
                    top = p
                    ptr = j
            col2.append(top)
            col2_ptr.append(ptr)
                
        probs.append(col2)
        ptrs.append(col2_ptr)

    # find max traceback
    j = util.argmax(probs[-1])
    traceback = [0] * n
    traceback[n-1] = j
    for i in xrange(n-1, 0, -1):
        j = ptrs[i][j]
        traceback[i-1] = j

    return traceback



def forward_algorithm(model, n, verbose=False):

    probs = []

    # calc first position
    nstates = model.get_num_states(0)
    probs.append([model.prob_prior(j) + model.prob_emission(0, j)
                  for j in xrange(nstates)])
    
    if n > 20:
        step = (n // 20)
    else:
        step = 1
    
    # loop through positions
    for i in xrange(1, n):
        if verbose and i % step == 0:
            print " forward iter=%d/%d, lnl=%f" % (i+1, n, max(probs[i-1]))

        nstates1 = model.get_num_states(i-1)
        nstates2 = model.get_num_states(i)
        col1 = probs[i-1]

        # find total transition and emission
        col2 = []
        for k in xrange(nstates2):
            tot = -util.INF
            emit = model.prob_emission(i, k)
            for j in xrange(nstates1):
                p = col1[j] + model.prob_transition(i-1, j, i, k) + emit
                tot = logadd(tot, p)
            col2.append(tot)
                
        probs.append(col2)

    return probs



def backward_algorithm(model, n, verbose=False):

    probs = []

    # calc last position
    nstates = model.get_num_states(n-1)
    for i in xrange(n):
        probs.append(None)
    probs[n-1] = [model.prob_prior(j) + model.prob_emission(n-1, j)
                  for j in xrange(nstates)]
    
    if n > 20:
        step = (n // 20)
    else:
        step = 1
    
    # loop through positions
    for i in xrange(n-2, -1, -1):
        if verbose and i % step == 0:
            print " backward iter=%d/%d, lnl=%f" % (i+1, n, max(probs[i+1]))

        nstates1 = model.get_num_states(i)
        nstates2 = model.get_num_states(i+1)
        col2 = probs[i+1]

        # find total transition and emission
        col1 = []
        for j in xrange(nstates1):
            tot = -util.INF
            for k in xrange(nstates2):
                p = (col2[k] + model.prob_emission(i+1, k) + 
                     model.prob_transition(i, j, i+1, k))
                tot = logadd(tot, p)
            col1.append(tot)
                
        probs[i] = col1

    return probs


def get_posterior_probs(model, n, verbose=False):

    probs_forward = forward_algorithm(model, n, verbose=verbose)
    probs_backward = backward_algorithm(model, n, verbose=verbose)

    total_prob = -util.INF
    for j in xrange(model.get_num_states(n-1)):
        total_prob = logadd(total_prob, probs_forward[n-1][j] +
                            model.prob_prior(j))

    probs_post = [
        [probs_forward[i][j] + probs_backward[i][j] - total_prob
         for j in xrange(model.get_num_states(i))]
        for i in xrange(n)]

    return probs_post


        
if __name__ == "__main__":
    from rasmus.common import *

    def trans(pos1, state1, pos2, state2):
        if state1 == state2:
            return log(.9)
        else:
            return log(.1)

    def emit(state):
        if state == 0:
            if random.random() < .9:
                return "H"
            else:
                return "T"
        elif state == 1:
            if random.random() < .1:
                return "H"
            else:
                return "T"
    
    def prob_emission(state, data):
        if state == 0:
            if data == "H":
                return log(.9)
            else:
                return log(.1)
        elif state == 1:
            if data == "H":
                return log(.1)
            else:
                return log(.9)
            

    model = HMM()
    model.set_callbacks(get_num_states=lambda pos: 2,
                        prob_prior=lambda state: log(.5),
                        prob_transition=trans,
                        emit=emit)

    ndata = 100

    s = sample_hmm_states(model)
    states = [s.next() for i in xrange(ndata)]
    p = plot(states, style="lines")

    s = sample_hmm_data(model, states)
    data = [s.next() for i in xrange(len(states))]
    print data

    model.prob_emission = lambda pos, state: prob_emission(state, data[pos])
    states2 = viterbi(model, len(data))

    probs = get_posterior_probs(model, len(data))
    states3 = [exp(probs[i][1]) for i in xrange(len(data))]

    p.plot(util.vadds(states2, 1.5), style="lines", miny=-1, maxy=4)
    p.plot(util.vadds(states3, 2.5), style="lines", miny=-1, maxy=4)

    


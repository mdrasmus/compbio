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


class HMM (object):
    def __init__(self, get_num_states=None,
                 prob_prior=None,
                 prob_emission=None,
                 prob_transition=None,
                 emit=None):
        self.get_num_states = get_num_states
        self.prob_prior = prob_prior
        self.prob_emission = prob_emission
        self.prob_transition = prob_transition
        self.emit = emit


def sample_hmm_first_state(model):
    state = 0
    nstates = model.get_num_states(0)
    p = model.prob_prior(state)
    pick = log(random.random())
    while pick > p and state < nstates:
        state += 1
        p = stats.logadd(p, model.prob_prior(state))
    return state


def sample_hmm_next_state(model, pos, state):
    nstates = model.get_num_states(pos)
    state2 = 0
    p = model.prob_transition(pos-1, state, pos, state2)
    pick = log(random.random())
    while pick > p and state2 < nstates:
        state2 += 1
        p = stats.logadd(p, model.prob_transition(pos-1, state, pos, state2))
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



def viterbi(model, n):

    probs = []
    ptrs = []

    # calc first position
    nstates = model.get_num_states(0)
    probs.append([model.prob_prior(j) for j in xrange(nstates)])
    ptrs.append([-1] * nstates)
    
    
    # loop through positions
    for i in xrange(1, n):
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
    for i in xrange(n-2, 0, -1):
        j = ptrs[i][j]
        traceback[i-1] = j

    return traceback
        
        
        
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
            

    model = HMM(get_num_states=lambda pos: 2,
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

    p.plot(util.vadds(states2, 1.5), style="lines", miny=-1, maxy=4)

    


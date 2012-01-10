
import unittest

from rasmus.common import *
from rasmus.testing import *

from rasmus import hmm


#=============================================================================

def make_coin_model(t=.1, e=.9):

    
    def trans(pos1, state1, pos2, state2):
        if state1 == state2:
            return log(1.0 - t)
        else:
            return log(t)

    def emit(pos, state):
        if state == 0:
            if random.random() < e:
                return "H"
            else:
                return "T"
        elif state == 1:
            if random.random() < 1 - e:
                return "H"
            else:
                return "T"
    
    def prob_emission_data(state, data):
        if state == 0:
            if data == "H":
                return log(e)
            else:
                return log(1-e)
        elif state == 1:
            if data == "H":
                return log(1-e)
            else:
                return log(e)


    model = hmm.HMM()
    model.set_callbacks(get_num_states=lambda pos: 2,
                        prob_prior=lambda state: log(.5),
                        prob_transition=trans,
                        emit=emit)
    model.prob_emission_data = prob_emission_data

    return model


    
class Test (unittest.TestCase):

    def test_coin(self):

        model = make_coin_model()

        # sample states
        ndata = 100
        states = list(islice(hmm.sample_hmm_states(model), ndata))
        p = plot(states, style="lines")

        # sample data
        data = list(hmm.sample_hmm_data(model, states))
        print data

        # viterbi
        model.prob_emission = (lambda pos, state:
            model.prob_emission_data(state, data[pos]))
        states2 = hmm.viterbi(model, len(data))

        # posterior
        probs = hmm.get_posterior_probs(model, len(data))
        states3 = [exp(probs[i][1]) for i in xrange(len(data))]

        # plot inference
        p.plot(util.vadds(states2, 1.5), style="lines", miny=-1, maxy=4)
        p.plot(util.vadds(states3, 2.5), style="lines", miny=-1, maxy=4)

        pause()

    def test_coin_sample_post(self):
        
        print "test sample posterior"
        model = make_coin_model()

        # sample states and data
        ndata = 100
        states = list(islice(hmm.sample_hmm_states(model), ndata))
        data = list(hmm.sample_hmm_data(model, states))
        model.prob_emission = (lambda pos, state:
            model.prob_emission_data(state, data[pos]))

        
        p = plot(states, style="lines")
        p.enableOutput(False)

        probs = hmm.get_posterior_probs(model, len(data))
        states2 = [exp(probs[i][1]) for i in xrange(len(data))]
        p.plot(util.vadds(states2, 1.5), style="lines", miny=-1, maxy=12)

        for i in range(2, 10):
            print i
            states2 = hmm.sample_posterior(model, ndata)
            p.plot(util.vadds(states2, 1.5*i), style="lines", miny=-1, maxy=12)
        p.enableOutput(True)
        p.replot()

        pause()



#=============================================================================
if __name__ == "__main__":   
    test_main()




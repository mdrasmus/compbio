
from itertools import islice
from math import exp
from math import log
import random
import unittest

from rasmus import stats
from rasmus import util
from rasmus.gnuplot import Gnuplot
from rasmus.testing import make_clean_dir

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
                        prob_prior=lambda pos, state: log(.5),
                        prob_transition=trans,
                        emit=emit)
    model.prob_emission_data = prob_emission_data

    return model


class Test (unittest.TestCase):

    def test_coin(self):
        """Test that viterbi and posterior coding work well."""

        outdir = 'test/tmp/test_hmm/test_coin/'
        make_clean_dir(outdir)

        model = make_coin_model()

        # sample states
        ndata = 100
        states = list(islice(hmm.sample_hmm_states(model), ndata))
        p = Gnuplot()
        p.enableOutput(False)
        p.plot(states, style="lines")

        # sample data
        data = list(hmm.sample_hmm_data(model, states))

        # viterbi
        model.prob_emission = (lambda pos, state:
                               model.prob_emission_data(state, data[pos]))
        states2 = hmm.viterbi(model, len(data))

        # posterior
        probs = hmm.get_posterior_probs(model, len(data))
        states3 = [exp(probs[i][1]) for i in xrange(len(data))]

        # assert that inferences correlates with true state
        self.assertTrue(stats.corr(states, states2) > .5)
        self.assertTrue(stats.corr(states, states3) > .5)

        # plot inference
        p.plot(util.vadds(states2, 1.5), style="lines", miny=-1, maxy=4)
        p.plot(util.vadds(states3, 2.5), style="lines", miny=-1, maxy=4)
        p.enableOutput(True)
        p.save(outdir + 'plot.png')

    def test_coin_sample_post(self):
        """Test sampling from posterior distribution"""

        outdir = 'test/tmp/test_hmm/test_coin_sample_post/'
        make_clean_dir(outdir)
        model = make_coin_model()

        # sample states and data
        ndata = 100
        states = list(islice(hmm.sample_hmm_states(model), ndata))
        data = list(hmm.sample_hmm_data(model, states))
        model.prob_emission = (lambda pos, state:
                               model.prob_emission_data(state, data[pos]))

        p = Gnuplot()
        p.enableOutput(False)
        p.plot(states, style="lines")

        probs = hmm.get_posterior_probs(model, len(data))
        states2 = [exp(probs[i][1]) for i in xrange(len(data))]
        p.plot(util.vadds(states2, 1.5), style="lines", miny=-1, maxy=12)

        for i in range(2, 10):
            states2 = hmm.sample_posterior(model, ndata)
            self.assertTrue(stats.corr(states, states2) > .5)

            p.plot(util.vadds(states2, 1.5*i), style="lines", miny=-1, maxy=12)
        p.enableOutput(True)
        p.save(outdir + 'plot.png')

    def test_coin_post(self):
        """Test that posterior decoding."""

        model = make_coin_model()

        # sample states and data
        ndata = 100
        states = list(islice(hmm.sample_hmm_states(model), ndata))
        data = list(hmm.sample_hmm_data(model, states))
        model.prob_emission = (lambda pos, state:
                               model.prob_emission_data(state, data[pos]))

        probs = hmm.get_posterior_probs(model, len(data))
        for col in probs:
            p = sum(map(exp, col))
            self.assertAlmostEqual(p, 1.0)

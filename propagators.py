from scipy.stats import norm, poisson
from scipy.misc import comb
from constants import *
from numpy import exp, log


class propagator:
    def __init__(self):
        ## initialization
        # p0 = round(signal[1]/alpha)            # initial concentration of the protein
        self.p0 = 15
        self.k_ON = 1. / 10.       # rate of switching to the ON state
        self.k_OFF = 1. / 10.       # rate of switching to the OFF state
        self.K = self.k_ON + self.k_OFF #
        self.k_s = 1.432        # synthesis (production) rate
        self.k_m = 0.0154       # degradation rate
        self.exp_km_q = exp(- self.k_m*delta_t*p_range)
        self.exp_km = np.power((1 - exp(- self.k_m*delta_t)), p_range)
        self.combs = np.zeros(numb_proteins**2).reshape(numb_proteins, numb_proteins)
        for n, p in enumerate(p_range):
            self.combs[n, ...] = map(lambda x: comb(n, x), p_range)
        self.pois = poisson.pmf(p_range, (self.k_s/self.k_m)*self.exp_km[1])
        self.vec_p_m_for_p = np.vectorize(self.p_m_for_p)
        self.P0 = np.zeros(numb_proteins**2, dtype=np.float).reshape(numb_proteins, numb_proteins)
        self.P1 = np.zeros(numb_proteins**2, dtype=np.float).reshape(numb_proteins, numb_proteins)
        for index, p_old in enumerate(p_range):
            self.P0[index, ...] = self.calculate_transitions(p_old, 0)
            self.P1[index, ...] = self.calculate_transitions(p_old, 1)
        self.G = np.zeros(4, dtype=np.float).reshape(2,2)
        for new_state in xrange(2):
            self.G[new_state, ...] = log(map(lambda g: self.p_g(new_state, g, delta_t), np.arange(0, 2)))


    def __init__(self, p0, k_ON, k_OFF, k_s, k_m):
        ## initialization
        # p0 = round(signal[1]/alpha)            # initial concentration of the protein
        self.p0 = p0
        self.k_ON = k_ON       # rate of switching to the ON state
        self.k_OFF = k_OFF       # rate of switching to the OFF state
        self.K = self.k_ON + self.k_OFF #
        self.k_s = k_s        # synthesis (production) rate
        self.k_m = k_m       # degradation rate
        self.exp_km_q = exp(- self.k_m*delta_t*p_range)
        self.exp_km = np.power((1 - exp(- self.k_m*delta_t)), p_range)
        self.combs = np.zeros(numb_proteins**2).reshape(numb_proteins, numb_proteins)
        for n, p in enumerate(p_range):
            self.combs[n, ...] = map(lambda x: comb(n, x), p_range)
        self.pois = poisson.pmf(p_range, (self.k_s/self.k_m)*self.exp_km[1])
        self.vec_p_m_for_p = np.vectorize(self.p_m_for_p)
        self.P0 = np.zeros(numb_proteins**2, dtype=np.float).reshape(numb_proteins, numb_proteins)
        self.P1 = np.zeros(numb_proteins**2, dtype=np.float).reshape(numb_proteins, numb_proteins)
        for index, p_old in enumerate(p_range):
            self.P0[index, ...] = self.calculate_transitions(p_old, 0)
            self.P1[index, ...] = self.calculate_transitions(p_old, 1)
        self.G = np.zeros(4, dtype=np.float).reshape(2,2)
        for new_state in xrange(2):
            self.G[new_state, ...] = log(map(lambda g: self.p_g(new_state, g, delta_t), np.arange(0, 2)))


    def p_g(self, g_new, g_old, t):
        return exp(- self.K * t) * (1 if g_new == g_old else 0) + \
            (1 - exp(- self.K * t)) * ( (self.k_ON if g_new == 1 else self.k_OFF) / self.K)


    def p_g_markov(self, t):
      n = np.power(1 - self.k_ON - self.k_OFF, t)
      x = self.k_ON + self.k_OFF
      return np.matrix( [[(self.k_OFF + self.k_ON*n)/x, (self.k_ON - self.k_ON*n)/x],
                         [(self.k_OFF - self.k_OFF*n)/x, (self.k_ON + self.k_OFF*n)/x]] )


    def p_m_for_p(self, q, n_new, n_old):
        return self.pois[n_new - q] * self.exp_km[n_old - q] * self.exp_km_q[q] * self.combs[n_old, q]


    def p_m(self, n_new, n_old, g, t):
        if g == 1:    # if the promoter is ON
            return log(np.sum(self.vec_p_m_for_p(np.arange(0, n_old+1), n_new, n_old)))
        else:         # when the promoter is OFF
            if n_new <= n_old:
                return log(self.combs[n_old, n_new]*self.exp_km_q[n_new]*self.exp_km[n_old - n_new])
            else:
                return -np.inf


    def calculate_transitions(self, p_old, g):
        return [self.p_m(p_new, p_old, g, delta_t) for p_new in p_range]

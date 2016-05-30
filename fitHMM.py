import numpy as np
from numpy import log, exp
from scipy.stats import norm, poisson
from scipy.misc import comb
import csv
from propagators import *
from constants import *
from Loglikelihood import loglikelihood
import scipy.optimize as optimize

T = 100
numb = np.zeros(T, np.int)
signal = np.zeros(T, dtype=np.float64)
time = np.zeros(T, dtype=np.int)
state = np.zeros(T, dtype=np.int)
with open("data.simulated") as inf:
    for i, rec in enumerate(csv.DictReader(inf, delimiter='\t')):
        numb[i] = int(rec["numb"])
        signal[i] = float(rec["signal"])
        time[i] = int(rec["time"])
        state[i] = int(rec["states"])

# propg = propagator(round(signal[1]/alpha), 1./5, 1./5, 2.5, 0.0154)

def p_s(s, p):
    return norm.pdf(s, p*alpha, np.sqrt(sigma_b2 + beta*p))

Obs = np.zeros(T*numb_proteins, dtype=np.float).reshape(T, numb_proteins)  # Observation matrix

for t in xrange(T):
    Obs[t, ...] = log(map(lambda p: p_s(signal[t], p), p_range))

## initial state distribution and transition matrix
pi  = np.zeros(numb_states, dtype=np.float128) - log(numb_states)     # uniform probability for initial states

def objective(x):
    print x
    if(x[0] > 1 or x[0] < 0):
        return np.Inf
    if(x[1] > 1 or x[1] < 0):
        return np.Inf
    if(x[2] > 5 or x[2] < 0):
        return np.Inf
    propg = propagator(round(signal[1]/alpha), x[0], x[1], x[2], 0.0154)   # half-life is fixed to 0.01
    logL, F = loglikelihood(T, pi, Obs, propg)
    print logL
    return -logL

res = optimize.fmin(objective, [1./5, 1./5, 2.5])
print res

exit()

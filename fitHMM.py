import numpy as np
from scipy.stats import norm, poisson
from scipy.misc import comb
import csv

p_min = 0
p_max = 150
numb_proteins = (p_max - p_min + 1)
numb_states = numb_proteins*2
T = 100
delta_t = 5
numb = np.zeros(T, np.int)
signal = np.zeros(T, dtype=np.float128)
time = np.zeros(T, dtype=np.int)
state = np.zeros(T, dtype=np.int)
with open("data.simulated") as inf:
    for i, rec in enumerate(csv.DictReader(inf, delimiter='\t')):
        numb[i] = int(rec["numb"])
        signal[i] = float(rec["signal"])
        time[i] = int(rec["time"])
        state[i] = int(rec["states"])


## protein to lumi signal parameters
alpha = 27.63
beta = 34
sigma_b2 = 50^2


## initialization
p0 = round(signal[1]/alpha)            # initial concentration of the protein
k_ON = .5       # rate of switching to the ON state
k_OFF = .5      # rate of switching to the OFF state
K = k_ON + k_OFF #
k_s = 2.5        # synthesis (production) rate
k_m = 0.01       # degradation rate


def p_s(s, p):
    return norm.pdf(s, p*alpha, np.sqrt(sigma_b2 + beta*p))


def p_g(g_new, g_old, t ):
    return np.exp(- K * t) * (1 if g_new == g_old else 0) + \
        (1 - np.exp(- K * t)) * ( (k_ON if g_new == 1 else k_OFF) / K )

def p_m_for_p(q, n_new, n_old, mean_n, k_m):
    return poisson.pmf(n_new - q, mean_n) * np.exp(- k_m*t*q) * \
           np.power(1 - np.exp(- k_m * t), n_old - q) * \
           comb(n_old, q)
vec_p_m_for_p = np.vectorize(p_m_for_p)

def p_m(n_new, n_old, g, t):
    k = k_s if g == 1 else 0
    mean_n = (k / k_m) * (1 - np.exp(- k_m * t ))
    return np.sum( vec_p_m_for_p(np.arange(0, n_old+1), n_new, n_old, mean_n, k_m) )
vec_p_m = np.vectorize(p_m)

## initial state distribution and transition matrix
pi  = np.ones(numb_states, dtype=np.float128) / numb_states     # assuming 2000 states, 1000 protein number when gene is OFF, and 1000 for when gene is ON
F = np.zeros(numb_states*T, dtype=np.float128).reshape(numb_states, T)            # Forward equations for the HMM
B = np.zeros(numb_states*T, dtype=np.float128).reshape(numb_states, T)            # Backward equations for the HMM
tmp = map(lambda p: p_s(signal[0], p), np.arange(p_min, p_max+1))
F[..., 1] = np.log(pi * np.tile(tmp, 2))                       # Forward equation for the first time-point
B[..., T-1] = np.zeros(numb_states)

G = np.zeros(4, dtype=np.float).reshape(2,2)
P0 = np.zeros(numb_proteins**2, dtype=np.float).reshape(numb_proteins,numb_proteins)
P1 = np.zeros(numb_proteins**2, dtype=np.float).reshape(numb_proteins,numb_proteins)
S = np.zeros(T*numb_proteins, dtype=np.float).reshape(T,numb_proteins)

for t in xrange(T):
    S[t, ...] = map(lambda p: p_s(signal[t], p), np.arange(p_min, p_max+1, dtype=np.float128))

for new_state in xrange(2):
    G[new_state, ...] = map(lambda g: p_g(new_state, g, delta_t), np.arange(0, 2))

for p_old in np.arange(p_min, p_max + 1):
    # P0[p_old, ...] = vec_p_m( np.arange(p_min, p_max+1), p_old, 0, delta_t )
    # P1[p_old, ...] = vec_p_m( np.arange(p_min, p_max+1), p_old, 1, delta_t )
    P0[p_old, ...] = [p_m(p_new, p_old, 0, delta_t) for p_new in np.arange(p_min, p_max+1)]
    P1[p_old, ...] = [p_m(p_new, p_old, 1, delta_t) for p_new in np.arange(p_min, p_max+1)]
    # P0[p_old, ...] = map(lambda p_new: p_m(p_new, p_old, 0, delta_t), np.arange(p_min, p_max+1))
    # P1[p_old, ...] = map(lambda p_new: p_m(p_new, p_old, 1, delta_t), np.arange(p_min, p_max+1))
    print p_old
print P0[10, ...]

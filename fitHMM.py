import numpy as np
from numpy import log, exp
from scipy.stats import norm, poisson
from scipy.misc import comb
import csv
import threading


p_min = 0
p_max = 100
numb_proteins = (p_max - p_min + 1)
numb_states = numb_proteins*2
p_range =np.arange(p_min, p_max+1)
T = 200
delta_t = 5
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


## protein to lumi signal parameters
alpha = 27.63
beta = 34
sigma_b2 = 50**2

## initialization
p0 = round(signal[1]/alpha)            # initial concentration of the protein
k_ON = .5       # rate of switching to the ON state
k_OFF = .5      # rate of switching to the OFF state
K = k_ON + k_OFF #
k_s = 2.4032        # synthesis (production) rate
k_m = 0.0154       # degradation rate


def p_s(s, p):
    return norm.pdf(s, p*alpha, np.sqrt(sigma_b2 + beta*p))


def p_g(g_new, g_old, t):
    return exp(- K * t) * (1 if g_new == g_old else 0) + \
        (1 - exp(- K * t)) * ( (k_ON if g_new == 1 else k_OFF) / K)

exp_km_q = exp(- k_m*delta_t*p_range)
exp_km = np.power((1 - exp(- k_m*delta_t)), p_range)
combs = np.zeros(numb_proteins**2).reshape(numb_proteins, numb_proteins)
for n, p in enumerate(p_range):
    combs[n, ...] = map(lambda x: comb(n, x), p_range)

pois = poisson.pmf(p_range, (k_s/k_m)*exp_km[1]) * exp_km


def p_m_for_p(q, n_new, n_old):
    return pois[n_new - q] * exp_km_q[q] * combs[n_old, q]
vec_p_m_for_p = np.vectorize(p_m_for_p)


def p_m(n_new, n_old, g, t):
    if g == 1:    # if the promoter is ON
        return log(np.sum(vec_p_m_for_p(np.arange(0, n_old+1), n_new, n_old)))
    else:         # when the promoter is OFF
        if n_new <= n_old:
            return log(combs[n_old, n_new]*exp_km_q[n_old]*exp_km[n_old - n_new])
        else:
            return -np.inf



G = np.zeros(4, dtype=np.float).reshape(2,2)
P0 = np.zeros(numb_proteins**2, dtype=np.float).reshape(numb_proteins,numb_proteins)
P1 = np.zeros(numb_proteins**2, dtype=np.float).reshape(numb_proteins,numb_proteins)
S = np.zeros(T*numb_proteins, dtype=np.float).reshape(T,numb_proteins)

for t in xrange(T):
    S[t, ...] = log(map(lambda p: p_s(signal[t], p), p_range))

for new_state in xrange(2):
    G[new_state, ...] = log(map(lambda g: p_g(new_state, g, delta_t), np.arange(0, 2)))

def calculate_transitions(p_old, g):
    return [p_m(p_new, p_old, g, delta_t) for p_new in p_range]

## initial state distribution and transition matrix
pi  = np.zeros(numb_states, dtype=np.float128) - log(numb_states)     # uniform probability for initial states
F = np.zeros(numb_states*T, dtype=np.float128).reshape(numb_states, T)            # Forward equations for the HMM
B = np.zeros(numb_states*T, dtype=np.float128).reshape(numb_states, T)            # Backward equations for the HMM
F[..., 0] = pi + np.tile(S[0, ...], 2)                       # Forward equation for the first time-point
B[..., T-1] = np.zeros(numb_states)


class myThread (threading.Thread):
    def __init__(self, p_old, P, g):
        threading.Thread.__init__(self)
        self.p_old = p_old
        self.P = P
        self.g = g
    def run(self):
        tmp = calculate_transitions(self.p_old, self.g)
        # Get lock to synchronize threads
        threadLock.acquire()
        # print "thread %d when g=%d" % (self.p_old, self.g)
        self.P[self.p_old, ...] = tmp
        # Free lock to release next thread
        threadLock.release()

# threadLock = threading.Lock()
# threads_P0 = []
# threads_P1 = []

for index, p_old in enumerate(p_range):
    # P0[p_old, ...] = [p_m(p_new, p_old, 0, delta_t) for p_new in np.arange(p_min, p_max+1)]
    # P1[p_old, ...] = [p_m(p_new, p_old, 1, delta_t) for p_new in np.arange(p_min, p_max+1)]
    P0[index, ...] = calculate_transitions(p_old, 0)
    P1[index, ...] = calculate_transitions(p_old, 1)
    # threads_P0.append(myThread(p_old, P0, 0))
    # threads_P0[index].start()
    # threads_P1.append(myThread(p_old, P1, 1))
    # threads_P1[index].start()
    # print index

# Wait for all threads to complete
# for t1, t2 in zip(threads_P0, threads_P1):
#     t1.join()
#     t2.join()


for t in xrange(1, T):
    for p in p_range:
        # to avoid considering states that are highly improbably, we only
        # consider a boundary that is set by p_std
        p_std = 50  # should be changed!
        p_lower_bound = 0 if (p-p_std) < 0 else (p-p_std)
        p_upper_bound = p_max if (p+p_std) > p_max else (p+p_std)
        # calculating Forward matrix
        F[p, t] = np.sum(exp(np.array(G[0, 0] + P0[p_lower_bound:(p+1), p] + F[p_lower_bound:(p+1), (t-1)],
                                      dtype=np.float128)))
        F[p, t] = log(np.sum(exp(np.array(G[1, 0] + P1[p_lower_bound:p_upper_bound, p] + \
                                          F[(p_lower_bound+p_max+1):(p_upper_bound+p_max+1), (t-1)],
                                          dtype=np.float128))) + F[p, t]) + S[t, p]
        F[p + p_max + 1, t] = np.sum(exp(np.array(G[0, 1] + P0[p_lower_bound:(p+1), p] + F[p_lower_bound:(p+1), (t-1)],
                                                  dtype=np.float128)))
        F[p + p_max + 1, t] = log(np.sum(exp(np.array(G[1, 1] + P1[p_lower_bound:p_upper_bound, p] +
                                                      F[(p_lower_bound+p_max+1):(p_upper_bound+p_max+1), (t-1)],
                                                      dtype=np.float128))) + F[p + p_max + 1, t]) + S[t, p]
        # calculating Backward matrix
        B[p, (T - t - 1)] = np.sum(exp(np.array(G[0, 0] + P0[p, p_lower_bound:(p+1)] + B[p_lower_bound:(p+1), T-t] +
                                                S[T-t, p_lower_bound:(p+1)]), dtype=np.float128))
        B[p, (T - t - 1)] = log(np.sum(exp(np.array(G[0, 1] + P0[p, p_lower_bound:(p+1)] +
                                                    B[(p_lower_bound+p_max+1):(p+p_max+2), T-t] +
                                                    S[T-t, p_lower_bound:(p+1)]), dtype=np.float128)) + B[p, (T - t - 1)])
        B[p + p_max + 1, (T - t - 1)] = np.sum(exp(np.array(G[1, 0] + P1[p, p_lower_bound:p_upper_bound] +
                                                            B[p_lower_bound:p_upper_bound, T-t] +
                                                            S[T-t, p_lower_bound:p_upper_bound]), dtype=np.float128))
        B[p + p_max + 1, (T - t - 1)] = log(np.sum(exp(np.array(G[1, 1] + P1[p, p_lower_bound:p_upper_bound] +
                                                                B[(p_lower_bound+p_max+1):(p_upper_bound+p_max+1), T-t] +
                                                                S[T-t, p_lower_bound:p_upper_bound]), dtype=np.float128))
                                            + B[p+p_max, (T - t - 1)])


print np.sum(np.exp(F[:, T-1])) , np.log(np.sum(np.exp(F[:, T-1])))
# print p_range[-1] + p_max
exit()

# for (t in 2:T) {
#   # print(t)
#   for (p in p_min:p_max) {
#     # to reduce the complexity, I only go one standard deviation around p value
#     # the underlying assumtion is that a very big change from one time point to the next is highly improbable
#     p_std = ceiling(sqrt(sigma_b2 + beta*p))
#     p_std = 30  # should be changed!
#     p_lower_bound = ifelse((p-p_std)<0, 0,  p-p_std)
#     p_upper_bound = ifelse((p+p_std)>p_max, p_max,  p+p_std)
#     # states where the gene is OFF
#     F[(p+1) ,t] = log( sum( exp(P0[(p+1), (p_lower_bound + 1):(p+1)] + G[1,1] +
#                                F[ (p_lower_bound + 1):(p+1) ,(t-1)] + S[t, (p+1)] )) +
#                        sum( exp(P0[(p+1), (p_lower_bound + 1):(p+1)] + G[1,2] +
#                                  F[ (p_lower_bound + p_max + 1):(p+p_max+1) ,(t-1)] + S[t, (p+1)] )))
#
#     B[(p+1), T-t+1] = log(sum( exp(P0[(p_lower_bound + 1):(p+1), (p+1)] + G[1,1] +
#                                B[ (p_lower_bound + 1):(p+1) ,(T-t+2) ] +S[(T-t+2), (p+1)] ) ) +
#                           sum( exp(P0[(p_lower_bound + 1):(p+1), (p+1)] + G[2,1] +
#                                  B[ (p_lower_bound + p_max + 1):(p+p_max+1) ,(T-t+2) ] + S[(T-t+2), (p+1)])))
#     # states where the gene is ON
#     F[(p + 2 + p_max) ,t] = log(sum( exp(P1[(p+1), (p_lower_bound + 1):(p_upper_bound + 1)] + G[2,1] +
#                                      F[ (p_lower_bound +1):(p_upper_bound + 1) ,(t-1)] + S[t, (p+1)] )) +
#                             sum( exp(P1[(p+1), (p_lower_bound + 1):(p_upper_bound + 1)] + G[2,2] +
#                                      F[ (p_lower_bound + p_max + 1):(p_upper_bound + p_max + 1) ,(t-1)] + S[t, (p+1)]) ))
#     B[(p+2+p_max), T-t+1] = log(sum( exp(P1[(p_lower_bound + 1):(p_upper_bound + 1), (p+1)] + G[1,2] +
#                                    B[ (p_lower_bound + 1):(p_upper_bound + 1) ,(T-t+2) ] + S[(T-t+2), (p+1)] ) ) +
#                             sum( exp(P1[(p_lower_bound + 1):(p_upper_bound + 1), (p+1)] + G[2,2] +
#                                    B[ (p_lower_bound + p_max + 1):(p_upper_bound + p_max + 1) ,(T-t+2) ] + S[(T-t+2), (p+1)] ) ) )
#
#   }
# }
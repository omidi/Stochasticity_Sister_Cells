import numpy as np
from scipy.linalg import expm
from scipy.stats import norm
from constants import *
import scipy.optimize as optimize

def state_rate_matrix(k_ON, k_OFF, k_s, k_d, p_max):
    I = np.identity(p_max)
    Y = - k_d * np.diag(np.arange(0, p_max)) + k_d * np.diag(np.arange(1, p_max), 1)
    K = - k_s * np.diag(np.concatenate([np.ones(p_max-1), [0]])) + \
        k_s * np.diag(np.ones(p_max-1), -1)


    M1 = np.hstack([-k_OFF*I + K + Y, k_ON*I])
    M2 = np.hstack([k_OFF*I, -k_ON*I + Y])
    return np.vstack([M1, M2])


def steady_state(M):
    M[-1, :] = np.ones(M.shape[1])
    b = np.zeros(M.shape[1])
    b[-1] = 1
    return np.linalg.solve(M , b)


def observation_matrix(signal, alpha, beta, sigma_b2, p_max):
    def p_s(s, p):
        return norm.pdf(s, p*alpha, np.sqrt(sigma_b2 + beta*p))

    T = len(signal)
    p_range = np.arange(0, p_max)
    Po = np.zeros(T*(p_max), dtype=np.float).reshape(p_max, T)  # Observation matrix
    for t in xrange(T):
        Po[:, t] = map(lambda p: p_s(signal[t], p), p_range)
    return np.vstack([Po, Po])


def log_likelihood_for_model(D, Po, Pt, P0):
    Nd = len(D)
    Ln = np.zeros(Nd, dtype=np.float128)
    Fi = np.multiply(Po[:, 0], P0)
    Ln[0] = np.sum(Fi)
    Fi = Fi / Ln[0]
    Ln[0] = np.log(Ln[0])

    for i in xrange(1, Nd):
        Fj = np.multiply(Po[:, i], np.dot(Pt, Fi))
        Ln[i] = np.sum(Fj)
        Fj = Fj / Ln[i];
        Ln[i] = np.log(Ln[i])
        Fi = Fj
    return np.sum(Ln)


def log_likelihood(D, x):
    p_max = 250
    delta_time = 5
    deg_rate = .00456
    Po = observation_matrix(signal, alpha, beta, sigma_b2, p_max)
    M = state_rate_matrix(x[0], x[1], x[2], deg_rate, p_max)
    # M = state_rate_matrix(.5, .5, 2, 0.0145, p_max)
    Pt = expm(M*delta_time)
    P0 = steady_state(M)
    logL = log_likelihood_for_model(signal, Po, Pt, P0)
    return logL



import csv
# T = 100
# numb = np.zeros(T, np.int)
# signal = np.zeros(T, dtype=np.float64)
# time = np.zeros(T, dtype=np.int)
# state = np.zeros(T, dtype=np.int)
# with open("data.simulated") as inf:
#     for i, rec in enumerate(csv.DictReader(inf, delimiter='\t')):
#         if i>=T:
#             break
#         numb[i] = int(rec["numb"])
#         signal[i] = float(rec["signal"])
#         time[i] = int(rec["time"])
#         state[i] = int(rec["states"])

input = "4.11 with mitosis.csv"
column = 90
column -= 90
signal = []
time = []
delta_t = 5
t = 0
with open(input, 'rU') as inf:
    for rec in csv.reader(inf, dialect=csv.excel_tab, delimiter=','):
        if rec[column] == "":
            break
        signal.append(float(rec[column]))
        time.append(t)
        t += delta_t
signal = np.array(signal)
time = np.array(time)
T = len(time)


def objective(x):
    print x
    if(x[0] > 5 or x[0] < 0):
        return np.Inf
    if(x[1] > 5 or x[1] < 0):
        return np.Inf
    if(x[2] > 10 or x[2] < 0):
        return np.Inf
    logL = log_likelihood(signal, x)
    print logL
    # res_file.write('\t'.join([
    #     str(x[0]),
    #     str(x[1]),
    #     str(x[2]),
    #     str(deg_rate),
    #     str(logL) + '\n',
    # ]))
    return -logL

x = [1./10, 1./5, 3.563]
print objective(x)
res = optimize.fmin(objective, x)
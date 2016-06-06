import numpy as np
from numpy import log, exp
from scipy.linalg import expm
from scipy.stats import norm, multivariate_normal
from constants import *
import scipy.optimize as optimize
import csv
import os
# from simanneal import Annealer
import random
import math


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
    Po = observation_matrix(signal, alpha, beta, sigma_b2, p_max)
    M = state_rate_matrix(x[0], x[1], x[2], deg_rate, p_max)
    # M = state_rate_matrix(.5, .5, 2, 0.0145, p_max)
    Pt = expm(M*delta_time)
    P0 = steady_state(M)
    logL = log_likelihood_for_model(signal, Po, Pt, P0)
    return logL

def arguments():
    """adding arguments for the script"""
    import argparse
    parser = argparse.ArgumentParser(description="""
    Take an input file that contains all the data, the column for the cell that is
    interested, then it performs a MCMC sampling over the parameter space. By this,
    we can find the posterior of each parameter.
    """)
    parser.add_argument('-i', '--input', dest='input', action='store',
                        type=str, required=True,
                        help='The input file, should be in CSV format')
    parser.add_argument('-c', '--column', dest='column', action='store',
                        type=int, required=True,
                        help="""The column for the cell. The first column is 1.""")
    parser.add_argument('-d', '--deg_rate', dest='deg_rate', action='store',
                        type=float, required=True,
                        help="""The mRNA degradation rate for the clone.""")
    parser.add_argument('-o', '--output', dest='dest', action='store',
                        type=str, required=False, default="",
                        help="""The destination directory""")
    parser.add_argument('-n', '--niter', dest='niter', action='store',
                        type=int, required=False, default=2000,
                        help="""Number of iteration in the MCMC sampling. Default is 2000.""")
    parser.add_argument('-a', '--accep_ratio', dest='accep_ratio', action='store',
                        type=float, required=False, default=.25,
                        help="""The acceptance ration, must be within (0, 1]. Default 0.2.""")
    args = parser.parse_args()

    if args.accep_ratio > 1 or args.accep_ratio<=0:
        print "Invalid value for the acceptance ratio!"
        exit()

    return args


def MCMC_sampling(x0, niter, accept_ratio, signal):
    nvar = 3
    lamb = 1. / 10000.
    Loglamb = log(lamb)
    def log_prior(x):
        return nvar*Loglamb - lamb*np.sum(x)

    X = np.zeros(nvar*niter, dtype=np.float).reshape(niter, nvar)
    L = np.zeros(niter, dtype=np.float)
    A = np.zeros(niter, dtype=np.float)

    Xprev = x0
    Lprev = log_likelihood(signal, x0)
    Lprop = Lprev

    l = 1
    S = np.matrix(np.identity(nvar) * [.1, .1, .1])
    M = np.matrix(log(x0))

    X[0,:] = x0
    L[0] = Lprev
    A[0] = 1

    print "iter %d:  %f %f %d" % (1, Lprop, Lprev, 1)

    for i in xrange(1, niter):
        rv = multivariate_normal(log(Xprev), np.dot(l, S))
        Xprop = exp(rv.rvs())
        print Xprop
        Lprop = log_likelihood(signal, Xprop)

        if math.isnan(Lprop):
            continue

        Qprop = -np.sum(log(Xprop))
        Qprev = -np.sum(log(Xprev))

        aprob = min(1, exp(Lprop+ log_prior(Xprop)+ Qprev - Lprev - log_prior(Xprev) - Qprop))

        r = random.random()

        if r < aprob:
            print "iter %d:  %f %f %d" % (i + 1, Lprop, Lprev, 1)
            Xprev = Xprop
            Lprev = Lprop

            X[i, :] = Xprop
            L[i] = Lprop
            A[i] = 1
        else:
            print "iter %d:  %f %f %d" % (i + 1, Lprop, Lprev, 0)
            X[i, :] = Xprev
            L[i] = Lprev
            A[i] = 0
        y = 1./(3.*np.sqrt(i + accept_ratio))
        l = exp(log(l) + y*(aprob - accept_ratio))
        S = S + y*( np.dot(np.transpose(log(Xprev) - M), (log(Xprev) - M)) - S)
        M = M + y * (log(Xprev) - M)
        print np.dot(l , S)
    return L, X, A, np.dot(l , S)


def generate_random_state(state):
    while True:
        r = random.uniform(-.01, .01)
        if (state[0] + r) > 0:
            state[0] = state[0] + r
            break

    while True:
        r = random.uniform(-.01, .01)
        if (state[1] + r) > 0:
            state[1] = state[1] + r
            break

    while True:
        r = random.uniform(-.05, .05)
        if (state[2] + r) > 0:
            state[2] = state[2] + r
            break
    return state


if __name__ == '__main__':
    args = arguments()
    signal = []
    time = []
    args.column -= 1
    delta_time = 5
    t = 0
    with open(args.input, 'rU') as inf:
        for rec in csv.reader(inf, dialect=csv.excel_tab, delimiter=','):
            if rec[args.column] == "":
                break
            signal.append(float(rec[args.column]))
            time.append(t)
            t += delta_time
    signal = np.array(signal)
    time = np.array(time)
    deg_rate = args.deg_rate
    T = len(time)
    p_max = int(max(signal)/alpha) + 30
    x0 = [0.05, 0.05, 0.15]

    while True:
        if math.isnan(log_likelihood(signal, x0)):
            x0 = generate_random_state(x0)
        else:
            break

    print x0, log_likelihood(signal, x0)
    L, X, A, C = MCMC_sampling(x0, args.niter, .25, signal)
    res_file_name = os.path.join(args.dest, "cell_%d_%s_mcmc_params" % (args.column+1, os.path.basename(args.input)))
    res_file = open(res_file_name, 'w')
    res_file.write('\t'.join(['k_ON', 'k_OFF', 'k_s', 'deg', 'logL', 'acc\n']))
    for i, state in enumerate(X):
        res_file.write('\t'.join([
            str(state[0]),
            str(state[1]),
            str(state[2]),
            str(args.deg_rate),
            str(L[i]),
            str(int(A[i])) + '\n',
        ]))

    res_file.close()

    res_file_name = os.path.join(args.dest, "cell_%d_%s_mcmc_cov_matrix" % (args.column+1, os.path.basename(args.input)))
    res_file = open(res_file_name, 'w')
    for i in xrange(C.shape[0]):
        res_file.write('\t'.join([str(C[i, j]) for j in xrange(C.shape[1])]))
        res_file.write('\n')
    res_file.close()
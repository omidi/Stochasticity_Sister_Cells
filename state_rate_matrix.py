import numpy as np
from scipy.linalg import expm
from scipy.stats import norm
from constants import *
import scipy.optimize as optimize
import csv
import os

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
    interested, and it tries to fit the parameters for the birth-death process.
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
    args = parser.parse_args()
    return args

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
    res_file_name = os.path.join(args.dest, "cell_%d_%s" % (args.column+1, os.path.basename(args.input)))
    res_file = open(res_file_name, 'w')
    res_file.write('\t'.join(['k_ON', 'k_OFF', 'k_s', 'deg', 'logL\n']))
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
        res_file.write('\t'.join([
            str(x[0]),
            str(x[1]),
            str(x[2]),
            str(deg_rate),
            str(logL) + '\n',
        ]))
        return -logL

    x = [1./2, 1./2, 3.]
    print objective(x)
    res = optimize.fmin(objective, x)
    res_file.close()
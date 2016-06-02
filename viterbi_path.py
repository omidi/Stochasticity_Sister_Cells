import numpy as np
from numpy import log, exp
from scipy.stats import norm, poisson
from scipy.misc import comb
import csv
from propagators import *
from constants import *
from Loglikelihood import loglikelihood
import scipy.optimize as optimize
import os
from state_rate_matrix import *


def arguments():
    """adding arguments for the script"""
    import argparse
    parser = argparse.ArgumentParser(description="""
    Takes as input a data file, a column for which we have modelled the expression,
    and the output file of the fitHMM.py.
    It then calculates the most likely path of states using the viterbi algorithm.
    """)
    parser.add_argument('-i', '--input', dest='input', action='store',
                        type=str, required=True,
                        help='The input file, should be in CSV format')
    parser.add_argument('-c', '--column', dest='column', action='store',
                        type=int, required=True,
                        help="""The column for the cell. The first column is 1.""")
    parser.add_argument('-f', '--hmm_fit', dest='hmm_fit', action='store',
                        type=str, required=True,
                        help="""The output file of the fitHMM.py script.""")
    parser.add_argument('-o', '--output', dest='dest', action='store',
                        type=str, required=False, default="",
                        help="""The destination directory""")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = arguments()
    with open(args.hmm_fit) as inf:
        last_line = inf.readlines()[-1].split()

    k_ON = float(last_line[0])
    k_OFF = float(last_line[1])
    k_s = float(last_line[2])
    k_d = float(last_line[3])

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
    T = len(time)
    p_max = int(max(signal)/alpha) + 30

    Po = observation_matrix(signal, alpha, beta, sigma_b2, p_max)
    M = state_rate_matrix(k_ON, k_OFF, k_s, k_d, p_max)
    Pt = expm(M*delta_time)
    P0 = steady_state(M)
    total_states = Pt.shape[0]
    Pi = np.ones(total_states, dtype=np.float128) / float(total_states)
    V = np.zeros(total_states*T, dtype=np.float128).reshape(total_states, T)
    P = np.zeros(total_states*T).reshape(total_states, T)
    V[:, 0] = np.multiply(Po[:, 0], P0)

    for t in np.arange(1, T):
        for s in np.arange(total_states):
            tmp = Po[:, t]*Pt[:, s]*V[:, t-1]
            x = np.argmax(tmp)
            P[s, t] = x
            V[s, t] = tmp[x]
    X = []
    x_t = np.argmax(V[:, -1])
    X.append(int(x_t))
    for t in np.arange(1, T)[::-1]:
        x_t = P[x_t, t]
        X.append(int(x_t))

    res_filename = os.path.join(args.dest, "viter_path_cell_%d_%s" % (args.column+1, os.path.basename(args.input)))
    with open(res_filename, 'w') as outf:
        outf.write('\t'.join([
            "time",
            "signal",
            "p_n",
            "state\n",
        ]))

        for t, x_t in enumerate(X[::-1]):
            outf.write('\t'.join([
                str(time[t]),
                str(signal[t]),
                str(x_t - (x_t / p_max)*p_max),
                str(x_t / p_max) + '\n',
            ]))

    os.system('Rscript plot_viterbi_res.R "%s" %d' % (res_filename, args.column + 1))

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


def p_s(s, p):
    return norm.pdf(s, p*alpha, np.sqrt(sigma_b2 + beta*p))


def main():
    args = arguments()
    signal = []
    time = []
    args.column -= 1
    t = 0
    with open(args.input, 'rU') as inf:
        for rec in csv.reader(inf, dialect=csv.excel_tab, delimiter=','):
            if rec[args.column] == "":
                break
            signal.append(float(rec[args.column]))
            time.append(t)
            t += delta_t
    signal = np.array(signal)
    time = np.array(time)
    T = len(time)

    # T = 100
    # numb = np.zeros(T, np.int)
    # signal = np.zeros(T, dtype=np.float64)
    # time = np.zeros(T, dtype=np.int)
    # state = np.zeros(T, dtype=np.int)
    # with open("data.simulated") as inf:
    #     for i, rec in enumerate(csv.DictReader(inf, delimiter='\t')):
    #         numb[i] = int(rec["numb"])
    #         signal[i] = float(rec["signal"])
    #         time[i] = int(rec["time"])
    #         state[i] = int(rec["states"])

    # propg = propagator(round(signal[1]/alpha), 1./5, 1./5, 2.5, 0.0154)

    Obs = np.zeros(T*numb_proteins, dtype=np.float).reshape(T, numb_proteins)  # Observation matrix

    for t in xrange(T):
        Obs[t, ...] = log(map(lambda p: p_s(signal[t], p), p_range))

    ## initial state distribution and transition matrix
    res_file_name = os.path.join(args.dest, "cell_%d_%s" % (args.column+1, os.path.basename(args.input)))
    deg_rate = args.deg_rate
    pi  = np.zeros(numb_states, dtype=np.float128) - log(numb_states)     # uniform probability for initial states
    deg_rate = 0.0154
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
        propg = propagator(round(signal[1]/alpha), x[0], x[1], x[2], deg_rate)   # half-life is fixed to 0.01
        logL, F = loglikelihood(T, pi, Obs, propg)
        print logL
        res_file.write('\t'.join([
            str(x[0]),
            str(x[1]),
            str(x[2]),
            str(deg_rate),
            str(logL) + '\n',
        ]))
        return -logL


    res = optimize.fmin(objective, [1./2, 1./2, 2.])
    res_file.close()

    print res


if __name__ == '__main__':
    main()



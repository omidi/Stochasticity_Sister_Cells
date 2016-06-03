import numpy as np
from scipy.linalg import expm
from scipy.stats import norm
from constants import *
import scipy.optimize as optimize
import csv
import os
# from simanneal import Annealer
import random


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
    def objective(x):
        # print x
        if(x[0] > 5 or x[0] < 0):
            return np.Inf
        if(x[1] > 5 or x[1] < 0):
            return np.Inf
        if(x[2] > 10 or x[2] < 0):
            return np.Inf
        logL = log_likelihood(signal, x)
        # print logL
        # res_file.write('\t'.join([
        #     str(x[0]),
        #     str(x[1]),
        #     str(x[2]),
        #     str(deg_rate),
        #     str(logL) + '\n',
        # ]))
        return -logL

    # np.random.seed(555)
    # print objective(x0)
    # res = optimize.fmin(objective, x0)
    # res = optimize.basinhopping(objective, x0, niter=100,
    #                             disp=True, niter_success=10)
    # rranges = (slice(0.0001, 4, 0.001), slice(0.001, 4, 0.001), slice(0.0001, 7, 0.05), )
    # res = optimize.brute(objective, rranges, full_output=True,
    #                      finish=optimize.fmin)

    class TelegramModel:
        def __init__(self, x0):
            self.temp = 500.0
            self.n_iter = 150.
            self.min_temp = .5
            # self.step = (self.temp - self.min_temp) / float(self.n_iter)
            self.step = self.n_iter
            self.state = x0
            self.energy = self.calculate_energy(x0)
            while self.energy == np.nan:
                self.state = self.move(self.state)
                self.energy = self.calculate_energy(self.state)
            self.best_state = self.state
            self.min_energy = self.energy

        def set_temp(self, new_temp):
            self.temp = new_temp

        def set_step(self, new_step):
            self.step = new_step

        def move(self):
            new_state = [0, 0, 0]
            while True:
                r = random.uniform(-.01, .01)
                if (self.state[0] + r) > 0:
                    new_state[0] = self.state[0] + r
                    break

            while True:
                r = random.uniform(-.01, .01)
                if (self.state[1] + r) > 0:
                    new_state[1] = self.state[1] + r
                    break

            while True:
                r = random.uniform(-.05, .05)
                if (self.state[2] + r) > 0:
                    new_state[2] = self.state[2] + r
                    break
            return new_state


        def calculate_energy(self, state):
            e = objective(state)
            return e


        def update(self, new_state):
            new_energy = self.calculate_energy(new_state)
            is_accepted = 0

            if new_energy == np.nan:
                return is_accepted

            if new_energy < self.energy:
                self.energy = new_energy
                self.state = new_state
                is_accepted = 1
                if new_energy < self.min_energy:
                    self.best_state = new_state
                    self.min_energy = new_energy
            else:
                r = random.random()
                if np.exp(- (new_energy - self.energy)/self.temp ) >= r:
                    self.energy = new_energy
                    self.state = new_state
                    is_accepted = 1

            print new_state, '%0.3f' % self.temp, \
                "Accepted" if is_accepted==1 else "Not_accepted", \
                '%0.3f' % self.energy
            self.temp = self.temp - self.temp/self.step
            return is_accepted


        def optimize(self):
            while self.temp >= self.min_temp:
                new_state = self.move()
                self.update(new_state)
            return self.state

        def give_best_fit(self):
            return self.best_state, self.min_energy

    x0 = [.05, .05, .15]

    model = TelegramModel(x0)
    state = model.optimize()
    best_state, min_energy = model.give_best_fit()

    res_file_name = os.path.join(args.dest, "cell_%d_%s" % (args.column+1, os.path.basename(args.input)))
    res_file = open(res_file_name, 'w')
    res_file.write('\t'.join(['k_ON', 'k_OFF', 'k_s', 'deg', 'logL\n']))
    res_file.write('\t'.join([
        str(best_state[0]),
        str(best_state[1]),
        str(best_state[2]),
        str(args.deg_rate),
        str(min_energy) + '\n',
    ]))
    print best_state, min_energy

    res = optimize.fmin(objective, best_state)
    res_file.write('\t'.join([
        str(res[0]),
        str(res[1]),
        str(res[2]),
        str(args.deg_rate),
        str(objective(res)) + '\n',
    ]))
    res_file.close()
    print best_state, min_energy
    # model.Tmax = 40000.0
    # model.steps = 2500
    # model.Tmin = 2.5
    # params, logL = model.anneal()
    #
    # print params
    # print logL
    # res = optimize.fmin(objective, params)


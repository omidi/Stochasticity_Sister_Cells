from constants import *
import numpy as np
from numpy import exp, log
from forward import forward

def loglikelihood(T, pi, S, prop):
    F = forward(T, pi, S, prop.P0, prop.P1, prop.G)
    return np.log(np.sum(np.exp(F[:, T-1]))), F


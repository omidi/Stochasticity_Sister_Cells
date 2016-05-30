from constants import *
import numpy as np
from numpy import exp, log


def forward(T, pi, S, P0, P1, G):
    F = np.zeros(numb_states*T, dtype=np.float128).reshape(numb_states, T)            # Forward equations for the HMM
    F[..., 0] = pi + np.tile(S[0, ...], 2)    # Forward equation for the first time-point
    for t in xrange(1, T):
        for p in p_range:
            # to avoid considering states that are highly improbably, we only
            # consider a boundary that is set by p_std
            p_std = 60  # should be changed!
            p_lower_bound = 0 if (p-p_std) < 0 else (p-p_std)
            p_upper_bound = p_max if (p+p_std) > p_max else (p+p_std)
            # calculating Forward matrix
            F[p, t] = np.sum(exp(np.array(G[0, 0] + P0[p_lower_bound:p_upper_bound, p] +
                                          F[p_lower_bound:p_upper_bound, (t-1)],
                                          dtype=np.float128)))
            F[p, t] = log(np.sum(exp(np.array(G[0, 1] + P1[p_lower_bound:p_upper_bound, p] + \
                                              F[(p_lower_bound+p_max+1):(p_upper_bound+p_max+1), (t-1)],
                                              dtype=np.float128))) + F[p, t]) + S[t, p]
            F[p + p_max + 1, t] = np.sum(exp(np.array(G[1, 0] + P0[p_lower_bound:p_upper_bound, p] +
                                                      F[p_lower_bound:p_upper_bound, (t-1)],
                                                      dtype=np.float128)))
            F[p + p_max + 1, t] = log(np.sum(exp(np.array(G[1, 1] + P1[p_lower_bound:p_upper_bound, p] +
                                                          F[(p_lower_bound+p_max+1):(p_upper_bound+p_max+1), (t-1)],
                                                          dtype=np.float128))) + F[p + p_max + 1, t]) + S[t, p]
            # calculating Backward matrix
            # B[p, (T - t - 1)] = np.sum(exp(np.array(G[0, 0] + P0[p, p_lower_bound:(p+1)] + B[p_lower_bound:(p+1), T-t] +
            #                                         S[T-t, p_lower_bound:(p+1)]), dtype=np.float128))
            # B[p, (T - t - 1)] = log(np.sum(exp(np.array(G[0, 1] + P0[p, p_lower_bound:(p+1)] +
            #                                             B[(p_lower_bound+p_max+1):(p+p_max+2), T-t] +
            #                                             S[T-t, p_lower_bound:(p+1)]), dtype=np.float128)) + B[p, (T - t - 1)])
            # B[p + p_max + 1, (T - t - 1)] = np.sum(exp(np.array(G[1, 0] + P1[p, p_lower_bound:p_upper_bound] +
            #                                                     B[p_lower_bound:p_upper_bound, T-t] +
            #                                                     S[T-t, p_lower_bound:p_upper_bound]), dtype=np.float128))
            # B[p + p_max + 1, (T - t - 1)] = log(np.sum(exp(np.array(G[1, 1] + P1[p, p_lower_bound:p_upper_bound] +
            #                                                         B[(p_lower_bound+p_max+1):(p_upper_bound+p_max+1), T-t] +
            #                                                         S[T-t, p_lower_bound:p_upper_bound]), dtype=np.float128))
            #                                     + B[p+p_max, (T - t - 1)])
    return F
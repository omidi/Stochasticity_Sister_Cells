## set working directory
setwd("~/PycharmProjects/Stochasticity_Sister_Cells/")

## protein to lumi signal parameters
alpha = 27.63
beta = 34
sigma_b2 = 50^2
#
# ## protein to lumi signal parameters
# alpha = 7.63
# beta = 54
# sigma_b2 = 100^2

p_s = function(s, p) {
  dnorm(s, p*alpha, sqrt(sigma_b2 + beta*p) )
}

p_g = function(g_new, g_old, t) {
  exp(- K * t) * ifelse(g_new == g_old, 1, 0) +
    (1 - exp(- K * t)) * ( ifelse(g_new == 1, k_ON, k_OFF) / K)
}

p_m = function(n_new, n_old, g, t) {
  k = ifelse(g==1, k_s, 0)
  mean_n = (k / k_m) * (1 - exp(- k_m * t))
  f = function(q) {
    dpois(n_new - q, mean_n) * exp(- k_m * t)^q *
      (1 - exp(- k_m * t))^(n_old - q) *
      choose(n_old, q)
  }
  sum(sapply(seq(0, n_old), f))
}

# loading the raw data
# df = read.csv('~/Dropbox/Communication/Gene.Expression.Inheritance/Data 29.04.16/4.11 with mitosis.csv',
#               header=FALSE)
# signal = df[, 1]
# signal = s[! is.na(signal) ]
# times = seq(1 , length(signal)*10, by = 10)
df = read.table("data.simulated", header=TRUE)
times = df$time
signal = df$signal

## initialization
p0 = round(signal[1]/alpha)            # initial concentration of the protein
k_ON = .5       # rate of switching to the ON state
k_OFF = .5      # rate of switching to the OFF state
K = k_ON + k_OFF #
k_s = 2.5        # synthesis (production) rate
k_m = 0.01       # degradation rate

# ## initialization
# p0 = round(signal[1]/alpha)         # initial concentration of the protein
# g0 = 0                              # initial gene/promoter state
# k_ON = 2                            # rate of switching to the ON state
# k_OFF = 2                           # rate of switching to the OFF state
# K = k_ON + k_OFF
# k_s = 4.43                          # synthesis (production) rate
# k_m = .254                          # degradation rate


p_max = 100
p_min = 0
T = length(times)
delta_t = times[2] - times[1]
T = 105
total_states = (p_max - p_min + 1)*2
## initial state distribution and transition matrix
pi  = rep(1/total_states , total_states)                # assuming 2000 states, 1000 protein number when gene is OFF, and 1000 for when gene is ON
F = matrix(0, nr=total_states, nc=T)        # Forward equations for the HMM
B = matrix(0, nr=total_states, nc=T)        # Backward equations for the HMM
F[ ,1] = log(pi * sapply(1:total_states, p_s, s=signal[1]))   # Forward equation for the first time-point
B[ ,T] = rep(0, total_states)

G = matrix(0 , nr=2 , nc=2)
total_protein = (p_max - p_min + 1)
P0 = matrix(0, nr=total_protein, nc=total_protein)
P1 = matrix(0, nr=total_protein, nc=total_protein)
S = matrix(0, nr=T, nc=total_protein)
for(t in 1:T) {
  S[t, ] = log(sapply(p_min:p_max, p_s, s=signal[t]))
}

for(new_state in 1:2)
  G[new_state, ] = log(sapply(0:1, p_g, g_new=new_state, t=delta_t))
for(p in p_min:p_max) {
  P0[(p+1), ] = log(sapply(p_min:p_max, p_m, n_old=p, t=delta_t, g=0))
  P1[(p+1), ] = log(sapply(p_min:p_max, p_m, n_old=p, t=delta_t, g=1))
}

for (t in 2:T) {
  print(t)
  for (p in p_min:p_max) {
    # to reduce the complexity, I only go one standard deviation around p value
    # the underlying assumtion is that a very big change from one time point to the next is highly improbable
    p_std = ceiling(sqrt(sigma_b2 + beta*p))
    p_std = 50
    p_lower_bound = ifelse((p-p_std)<0, 0,  p-p_std)
    p_upper_bound = ifelse((p+p_std)>p_max, p_max,  p+p_std)
    # states where the gene is OFF
    F[(p+1) ,t] = log( sum( exp(P0[(p+1), (p_lower_bound + 1):(p+1)] + G[1,1] +
                               F[ (p_lower_bound + 1):(p+1) ,(t-1)] + S[t, (p+1)] )) +
                      sum( exp(P0[(p+1), (p_lower_bound + 1):(p+1)] + G[1,2] +
                                 F[ (p_lower_bound + p_max + 1):(p+p_max+1) ,(t-1)] + S[t, (p+1)] )))

    B[(p+1), T-t+1] = log(sum( exp(P0[(p_lower_bound + 1):(p+1), (p+1)] + G[1,1] +
                               B[ (p_lower_bound + 1):(p+1) ,(T-t+2) ] +S[(T-t+2), (p+1)] ) ) +
                          sum( exp(P0[(p_lower_bound + 1):(p+1), (p+1)] + G[2,1] +
                                 B[ (p_lower_bound + p_max + 1):(p+p_max+1) ,(T-t+2) ] + S[(T-t+2), (p+1)])))
    # states where the gene is ON
    F[(p + 1 + p_max) ,t] = log(sum( exp(P1[(p+1), (p_lower_bound + 1):(p_upper_bound + 1)] + G[2,1] +
                                     F[ (p_lower_bound +1):(p_upper_bound + 1) ,(t-1)] + S[t, (p+1)] )) +
                            sum( exp(P1[(p+1), (p_lower_bound + 1):(p_upper_bound + 1)] + G[2,2] +
                                     F[ (p_lower_bound + p_max + 1):(p_upper_bound + p_max + 1) ,(t-1)] + S[t, (p+1)]) ))
    B[(p+1+p_max), T-t+1] = log(sum( exp(P1[(p_lower_bound + 1):(p_upper_bound + 1), (p+1)] + G[1,2] +
                                   B[ (p_lower_bound + 1):(p_upper_bound + 1) ,(T-t+2) ] + S[(T-t+2), (p+1)] ) ) +
                            sum( exp(P1[(p_lower_bound + 1):(p_upper_bound + 1), (p+1)] + G[2,2] +
                                   B[ (p_lower_bound + p_max + 1):(p_upper_bound + p_max + 1) ,(T-t+2) ] + S[(T-t+2), (p+1)] ) ) )

  }
}







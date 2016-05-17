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
k_m = 0.25       # degradation rate

# ## initialization 
# p0 = round(signal[1]/alpha)         # initial concentration of the protein
# g0 = 0                              # initial gene/promoter state
# k_ON = 2                            # rate of switching to the ON state
# k_OFF = 2                           # rate of switching to the OFF state
# K = k_ON + k_OFF  
# k_s = 4.43                          # synthesis (production) rate
# k_m = .254                          # degradation rate


p_max = 25
p_min = 0 
total_states = (p_max - p_min)*2 
## initial state distribution and transition matrix
pi  = rep(1/total_states , total_states)                # assuming 2000 states, 1000 protein number when gene is OFF, and 1000 for when gene is ON
F = matrix(0, nr=total_states, nc=length(times))        # Forward equations for the HMM
F[,1] = pi * sapply(1:total_states, p_s, s=signal[1])   # Forward equation for the first time-point

for (p in p_min:p_max) {
  p_std = ceiling(sqrt(sigma_b2 + beta*p))
  p_min = ifelse((p-p_std)<0, 0,  p-p_std)
  p_max = ifelse((p+p_std)>p_max, p_max,  p+p_std)
  x = sapply( p_min:p_max, p_m, n_new=p, g=1, t=2)
  F[p+1 ,2] = sum( (sapply( p_min:p_max, p_m, n_new=p, g=1, t=2) * p_g(1, 0, 2)) * F[(p_min+1):(p_max+1), 2-1] * sapply(p_min:p_max, p_s, s = signal[2] ))
}





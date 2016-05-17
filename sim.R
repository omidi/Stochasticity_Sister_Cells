## set working directory 
setwd("~/PycharmProjects/Stochasticity_Sister_Cells/")

## protein to lumi signal parameters 
alpha = 27.63
beta = 34 
sigma_b2 = 50^2

## initialization 
p0 = 14            # initial concentration of the protein
k_ON = .35       # rate of switching to the ON state
k_OFF = .05      # rate of switching to the OFF state
K = k_ON + k_OFF # 
k_s = 3.432      # synthesis (production) rate
k_m = 0.254       # degradation rate

p_g = function(g1, g2, t) {
  exp(- K * t) * ifelse(g1 == g2, 1, 0) + 
    (1 - exp(- K * t)) * ( ifelse(g2 == 1, k_ON, k_OFF) / K)
}

p_m = function(n1, n2, g, t) {
  k = g * k_s
  mean_n = (k / k_m) * (1 - exp(- k_m * t))
  f = function(q) {
    dpois(n1 - q, mean_n) * exp(- k_m * t)^q * 
      (1 - exp(- k_m * t))^(n2 - q) * 
      choose(n2, q)
  }
  sum(sapply(seq(0, n2), f))
}

curr_state = 1
curr_num = p0
states = c()
numbs = c()
t = 0
times = c()
while (t < 100*10) {
  states = c(states, curr_state)
  numbs = c(numbs, curr_num)
  times = c(times, t)
  r = runif(1)
  if(r < p_g(1, curr_state, t)) {
    curr_state = 1
  } else {
    curr_state = 0
  }
  n = seq(0, 30)
  s = sapply(n, p_m, n2 = curr_num, g = curr_state, t = 1)
  r = runif(1)
  curr_num = n[cumsum(s) >= r][1]
  t = t + 5
}

# states[which(states == 0)] = -20
# states[which(states == 1)] = 0
plot(states, ylim = c(0, 30), pch=15, cex= 3.1)
lines(numbs, type = 'l')
signal = sapply(numbs, function(x) rnorm(1, alpha*x, sqrt(sigma_b2 + beta*x)) )

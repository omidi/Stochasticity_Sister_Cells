## set working directory 
setwd("~/PycharmProjects/Stochasticity_Sister_Cells/")

## initialization 
p0 = 6            # initial concentration of the protein
k_ON = .3       # rate of switching to the ON state
k_OFF = .12      # rate of switching to the OFF state
K = k_ON + k_OFF # 
k_s = 4.43        # synthesis (production) rate
k_m = .254      # degradation rate

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
while (t < 100) {
  states = c(states, curr_state)
  numbs = c(numbs, curr_num)
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
  t = t + 1
}

states[which(states == 0)] = -20
states[which(states == 1)] = 0
plot(states, ylim = c(0, 30), pch=15, cex= 3.1)
lines(numbs, type = 'l')

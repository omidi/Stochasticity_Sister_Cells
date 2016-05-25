## set working directory
setwd("~/PycharmProjects/Stochasticity_Sister_Cells/")

## protein to lumi signal parameters
alpha = 27.63
beta = 34
sigma_b2 = 50^2

## initialization
p0 = 14            # initial concentration of the protein
k_ON = .6       # rate of switching to the ON state
k_OFF = .8      # rate of switching to the OFF state
K = k_ON + k_OFF
k_s = 2.432       # synthesis (production) rate
k_m = 0.0154       # degradation rate

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

curr_state = 1
curr_num = p0
states = c()
numbs = c()
t = 0
times = c()
while (t < 200*5) {
  states = c(states, curr_state)
  numbs = c(numbs, curr_num)
  times = c(times, t)
  r = runif(1)
  if(r < p_g(1, curr_state, 5)) {
    curr_state = 1
  } else {
    curr_state = 0
  }
  n = seq(0, 200)
  s = sapply(n, p_m, n_old = curr_num, g = curr_state, t = 5)
  r = runif(1)
  curr_num = n[cumsum(s) >= r][1]
  t = t + 5
}

# states[which(states == 0)] = -50
# states[which(states == 1)] = 0
plot(states, ylim = c(0, max(numbs)), pch=15, cex=.5)
lines(numbs, type = 'l')
signal = sapply(numbs, function(x) rnorm(1, alpha*x, sqrt(sigma_b2 + beta*x)))
df = data.frame(time=times, signal=signal, numb=numbs, states=states)

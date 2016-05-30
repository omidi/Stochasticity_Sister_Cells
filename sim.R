## set working directory
setwd("~/PycharmProjects/Stochasticity_Sister_Cells/")
require(scales)

## protein to lumi signal parameters
alpha = 27.63
beta = 34
sigma_b2 = 50^2

## initialization
p0 = 14            # initial concentration of the protein
k_ON = 1/10       # rate of switching to the ON state
k_ON = 0.17571506  
k_OFF = 1/5      # rate of switching to the OFF state
k_OFF = 0.36363595  
K = k_ON + k_OFF
k_s = 2.432       # synthesis (production) rate
k_s = 3.563   
k_m = 0.0154       # degradation rate

p_s = function(s, p) {
  dnorm(s, p*alpha, sqrt(sigma_b2 + beta*p) )
}

# p_g = function(g_new, g_old, t) {
#   exp(- K * t) * ifelse(g_new == g_old, 1, 0) +
#     (1 - exp(- K * t)) * ( ifelse(g_new == 1, k_ON, k_OFF) / K)
# }

# P_g = function(n) {
#   t = (1 - k_ON - k_OFF)^n
#   x = k_ON + k_OFF
#   matrix(c((k_ON + k_OFF*t)/x, (k_ON - k_ON*t)/x,
#            (k_OFF - k_OFF*t)/x, (k_ON + k_OFF*t)/x), nr=2, nc=2, byrow = TRUE)
# }

p_n = function(g_new, g_old, t) {
  exp(- K * t) * ifelse(g_new == g_old, 1, 0) +
    (1 - exp(- K * t)) * ( ifelse(g_new == 1, k_ON, k_OFF) / K)
}

P_g = function(n) {
  matrix(c(c(p_n(0, 0, n), p_n(1, 0, n), 
             c(p_n(0, 1, n) , p_n(1, 1, n) ))) , nr =2, nc=2, byrow=TRUE)
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

T= 100
curr_state = 0
curr_num = p0
states = rep(0, T)
numbs = rep(0, T)
signal = rep(0, T)
t = 0
times = c()
G = P_g(5)
print(G)
while (t < T*5) {
  states[t/5 + 1] = curr_state
  numbs[t/5 + 1] = curr_num
  signal[t/5 + 1] =  rnorm(1, curr_num*alpha, sqrt(sigma_b2 + beta*curr_num)) 
  times[t/5 + 1] = t
  r = runif(1)
  if(r < G[curr_state + 1, 2]) {
    curr_state = 1
  } else {
    curr_state = 0
  }
  n = seq(0, 300)
  s = sapply(n, p_m, n_old = curr_num, g = curr_state, t = 5)
  r = runif(1)
  curr_num = n[cumsum(s) >= r][1]
  t = t + 5
  print(t)
}

# states[which(states == 0)] = -50
# states[which(states == 1)] = 0
par(mar=c(5,5,5,2))
plot(times, states, ylim = c(0, max(signal)), pch=15, cex=.5, col=0, 
     xlab= "time (min)", ylab="signal", 
     cex.lab=1.5, cex.axis=1.4)
# abline(v=times[which(states==1)])
lines(smooth.spline(times, signal, spar = .4), type = 'l', lwd=3)
points(times, signal)

i = 1
while (TRUE) {
  if(states[i] == 0) {
    i = i + 1
  } else {
    start = times[i] 
    while(states[i] == 1) {
      i = i+1
      if(i >= length(states))   break
    } 
    end = times[i] 
    if(is.na(end)) {
      end = times[i - 1]
    }
    polygon(c(start, start, end, end), c(-20, 10000, 10000, -5), 
            col=alpha("grey", .3), border = alpha("grey", .3))
  }
  if(i >= length(times))  break
}

# signal = sapply(numbs, function(x) rnorm(1, alpha*x, sqrt(sigma_b2 + beta*x)))
df = data.frame(time=times, signal=signal, numb=numbs, states=states)
# write.table(df, "data.simulated", row.names = FALSE, col.names = TRUE, quote = FALSE, sep='\t')

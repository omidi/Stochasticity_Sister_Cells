a = .5  # Off to On
b = .1  # On to Off
P = matrix(c(1-a, a, b, 1-b), nr=2, nc=2, byrow = TRUE)

P_n = function(n) {
  t = (1 - a - b)^n
  x = a + b
  matrix(c((b + a*t)/x, (a - a*t)/x,
           (b - b*t)/x, (a + b*t)/x), nr=2, nc=2, byrow = TRUE)
}

G_t = function(t) {
  k = a + b
  e = exp(-k * t)
  matrix(c(e + (1 - e)*(1-a)/k,
           (1 - e)*a/k,
           (1 - e)*b/k,
           e + (1 - e)*(1-b)/k), nr=2, nc=2, byrow = TRUE)
}


# curr_state = 1
# states = c(curr_state)
# for(t in seq(1, 200)) {
#   r = runif(1)
#   P = P_n(3)
#   if(r < P[curr_state, 1]) {
#     curr_state = 1
#   } else {
#     curr_state = 2
#   }
#   states = c(states, curr_state)
# }
#
# plot(states, col=0)
# points(which(states==2),  rep(1, length(which(states==2))), pch=15, cex=1.)
# points(which(states==1),  rep(1, length(which(states==1))), pch=15, cex=1., col=0)
#

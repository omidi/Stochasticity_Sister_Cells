require(GenSA)

log_likelihood_data = function(X, vi, ve, corr){
  S = matrix(c(vi, 0, 0, vi), nc=2, nr=2)
  C = matrix(c(ve, corr, corr, ve), nr=2, nc=2)
  C.inv = solve(C)
  S.inv = solve(S)
  M = solve(S + S %*% C.inv %*% S)
  L = -.5*apply(X, 1, function(x) t(x) %*% S.inv %*% x) - .5*log(det(S)) - .5*log(det(C)) + 
    .5*log(det(solve(C.inv + S.inv))) + .5* apply(X, 1, function(x) (t(x) %*% M %*% x))  
  sum(L)
}


maximum.likelihood.simple.model <- function(data) {
  n.it=5000   # number of iterations in the Simulated Annealing
  epsilon = 1e-3   # criterion for stopping the fitting algorithm
  temperature = 3000  # the initial temperature
  acc.rate = .25  # acceptance rate
  fn <- function(params, names) {
    -log_likelihood_data(data, params[1], params[2], params[3])
  }
  C = cov(data)
  lower = c(vi=1, ve=min(C[1,1], C[1,2]) - 1000000, corr=C[1,2]-1000000)
  upper = c(vi=var(abs(X[,1] - X[,2])) + 1000000, ve=max(C[1,1], C[2,2])+1000000, corr=(C[1,2]) + 1000000)
  params = c(vi=var(abs(X[,1] - X[,2])), ve=max(C[1,1], C[2,2]), corr=C[1,2])
  print(params)
  out = GenSA(par=params, fn=fn, lower=lower, upper=upper, names=names(params), 
              control=list(maxit=n.it, verbose=TRUE, nb.stop.improvement=epsilon, 
                           acceptance.param=acc.rate, temperature=temperature))
  params = as.list(out$par)
  names(params) = c("vi", "ve", "corr")
  print(params)
  list(sigma=matrix(c(params[[2]],params[[3]],params[[3]],params[[2]]),nrow=2), intr=params[[1]])
}


plot_contour_for_one_time_point = function(df, noise.df, time, name='Jarid2') {
  sigma = matrix(c( noise.df[time, 2], noise.df[time, 3], 
                    noise.df[time, 3], noise.df[time, 2]), nc=2, nr=2) 
  print(sigma)
  m = matrix(colMeans(df[time:(time+1), ]), nc=2, byrow = TRUE)  
  m = m - mean(df)
  # sigma = cov(m)
  x.points <- seq(-7000,7000,length.out=100)
  y.points <- x.points
  z = matrix(0,nrow=100,ncol=100)
  for (i in 1:100) {
    for (j in 1:100) {
      z[i,j] <- dmvnorm(c(x.points[i],y.points[j]),
                        mean=c(0, 0),sigma=sigma)
    }
  }
  contour(x.points,y.points,z, col=1, lwd=2, xlab="Sister cell 1", ylab="Sister cell 2", 
          cex.axis=1.6, cex.lab=1.7, main=sprintf('%s time: %i', name, time), cex.main=2, 
          nlevels=10)
  points(m, pch=16, cex=1.3, col=2)
  abline(h=0, v=0)
  abline(0, 1, col='darkgrey')
  return(sigma)  
}

plot_contour_all_timepoints = function(df, noise.df) {
  df.n = normalize_time(df)
  
  for (index in seq(1, 99)) {
    if(index>=10) {
      fname = 'images_contours/frame_%i.png'
    } else {
      fname = 'images_contours/frame_0%i.png'
    }
    print(index)
    png(sprintf(fname, index), width = 700, height = 700) 
    par(mar=c(5,5,5,2))
    plot_contour_for_one_time_point(df.n, noise.df, index)
    dev.off()
  }
}






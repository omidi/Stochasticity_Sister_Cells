m = matrix(colMeans(df, na.rm = TRUE), nc=2, byrow = TRUE)
df.norm = sweep(df.n, 2, apply(df.n, 2, max), FUN='/')


plot_mean_expression_correlation = function(df, name=NA){
  m = matrix(colMeans(df, na.rm = TRUE), nc=2, byrow = TRUE) 
  par(mar=c(7,7,4,2))
  plot(m, pch=19, ylim = range(m), xlim = range(m), 
       cex.axis=1.7, cex.lab=1.8, xlab="average signal sister 1", ylab="average signal sister 2", 
       bty='n', cex=1.4, main=name, cex.main=2)
  
  # mtext(side = 1, text = "Mean signal sister 2", line = 5, cex = 1.6)
  # mtext(side = 2, text = "Mean signal sister 2", line = 5, cex = 1.6)
  
  text(max(m)/1.3, min(m) + (max(m) - min(m))/10, sprintf("Correlation coeff.: %0.2f", cor(m)[1, 2]), cex=1.7)
  abline(0, 1, lwd=2, lty=2)
  # segments(2000, 2000, 10000, 10000)
}


plot_mean_expression_correlation_all_times = function(df, interval=10, name=NA){
  cols = brewer.pal(8, "Set1")
  m = matrix(colMeans(df, na.rm = TRUE), nc=2, byrow = TRUE) 
  avg_corr = cor(m)[1, 2]
  cell_cycle_len = sapply(1:dim(df)[2], function(col) length(df[!is.na( df[, col]), col]))
  cc_class = ifelse(cell_cycle_len < 168, 1, ifelse(cell_cycle_len<216, 2, 3))
  df = normalize_time(df)
  require(RColorBrewer)
  print(cell_cycle_len)
  index = 1
  corr.vals = c()
  xvals = c()
  for (i in seq(1, 100 - interval/2 - 1, interval/2)) {
    # pdf(sprintf('images/frame_%i.pdf', index), height = 6, width = 6)
    if(index>=10) {
      fname = 'images/frame_%i.png'
    } else {
      fname = 'images/frame_0%i.png'
    }
    png(sprintf(fname, index), width = 1000, height = 500)
    par(mfrow=c(1, 2), mar=c(5,5,5,2))
    df2 = df[i:(i+interval-1), ]
    m = matrix(colMeans(df2, na.rm = TRUE), nc=2, byrow = TRUE)
    corr.vals = c(corr.vals, cor(m)[1, 2])
    xvals = c(xvals, i+interval/2)
    plot( xvals, corr.vals, pch=16, cex=1.6, xlab="% Cell cycle" , ylab = "Correlation coeff.", 
          cex.axis=1.6, cex.lab= 1.7, main=name, cex.main=2, bty='n', xlim=c(0,100), ylim=c(.2, 1), col=0)
    points(xvals[1], corr.vals[1], pch=16, cex=1.2)
    lines(xvals, corr.vals, lty=1, lwd=4)
    abline(h=avg_corr, lty=2, lwd=2, col="grey")
    plot(m, pch=19, ylim = range(df), xlim = range(df), 
        cex.axis=1.7, cex.lab=1.8, xlab="Mean signal sister 1", ylab="Mean signal sister 2", 
        bty='n', cex=1.4, cex.main=2, main=sprintf("%s  (%0.0f-%0.0f)%% of cell cycle", name, (index-1)*interval/2 + 1, (index-1)*interval/2 + interval), 
        cex.main=2, col=cols[cc_class])
    # mtext(side = 1, text = "Mean signal sister 2", line = 5, cex = 1.6)
    # mtext(side = 2, text = "Mean signal sister 2", line = 5, cex = 1.6)
    # text(max(m)/1.3, min(m) + (max(m) - min(m))/10, sprintf("Correlation coeff.: %0.2f", cor(m)[1, 2]), cex=1.3)
    text(13000, 1000, sprintf("Corr: %0.2f", cor(m)[1, 2]), cex=1.6)
    abline(0, 1, lwd=2, lty=2)
    # segments(2000, 2000, 10000, 10000)
    legend("topleft", c("< 14 hr", "< 18 hr & > 14 hr", "> 18 hr"), col = cols[1:3], pch=16, cex=1.6)
    dev.off()
    index = index + 1
    
  }
}


plot_eigne_decomp_expression_correlation_all_time = function(df, interval=10, name=NA){
  cols = brewer.pal(8, "Set1")
  m = matrix(colMeans(df, na.rm = TRUE), nc=2, byrow = TRUE) 
  avg_corr = cor(m)[1, 2]
  cell_cycle_len = sapply(1:dim(df)[2], function(col) length(df[!is.na( df[, col]), col]))
  cc_class = ifelse(cell_cycle_len < 168, 1, ifelse(cell_cycle_len<216, 2, 3))
  df = normalize_time(df)
  require(RColorBrewer)
  print(cell_cycle_len)
  index = 1
  corr.vals = c()
  xvals = c()
  for (i in seq(1, 100 - interval/2 - 1, interval/2)) {
    # pdf(sprintf('images/frame_%i.pdf', index), height = 6, width = 6)
    if(index>=10) {
      fname = 'images/eigen_frame_%i.png'
    } else {
      fname = 'images/eigen_frame_0%i.png'
    }
    png(sprintf(fname, index), width = 1000, height = 500)
    par(mfrow=c(1, 2), mar=c(5,5,5,2))
    df2 = df[i:(i+interval-1), ]
    m = matrix(colMeans(df2, na.rm = TRUE), nc=2, byrow = TRUE)
    vals = eigen(cov(m))$values
    corr.vals = c(corr.vals, vals[2])
    xvals = c(xvals, vals[1])
    plot( xvals, corr.vals, pch=16, cex=1.6, xlab="% Cell cycle" , ylab = "Correlation coeff.", 
          cex.axis=1.6, cex.lab= 1.7, main=name, cex.main=2, bty='n', col=0, 
          xlim=c(4.8e+06, 1.2e+07))
    points(xvals[1], corr.vals[1], pch=16, cex=1.2)
    lines(xvals, corr.vals, lty=1, lwd=4)
    abline(h=avg_corr, lty=2, lwd=2, col="grey")
    plot(m, pch=19, ylim = range(df), xlim = range(df), 
         cex.axis=1.7, cex.lab=1.8, xlab="Mean signal sister 1", ylab="Mean signal sister 2", 
         bty='n', cex=1.4, cex.main=2, main=sprintf("%s  (%0.0f-%0.0f)%% of cell cycle", name, (index-1)*interval/2 + 1, (index-1)*interval/2 + interval), 
         cex.main=2, col=cols[cc_class])
    # mtext(side = 1, text = "Mean signal sister 2", line = 5, cex = 1.6)
    # mtext(side = 2, text = "Mean signal sister 2", line = 5, cex = 1.6)
    # text(max(m)/1.3, min(m) + (max(m) - min(m))/10, sprintf("Correlation coeff.: %0.2f", cor(m)[1, 2]), cex=1.3)
    text(13000, 1000, sprintf("Corr: %0.2f", cor(m)[1, 2]), cex=1.6)
    abline(0, 1, lwd=2, lty=2)
    # segments(2000, 2000, 10000, 10000)
    legend("topleft", c("< 14 hr", "< 18 hr & > 14 hr", "> 18 hr"), col = cols[1:3], pch=16, cex=1.6)
    dev.off()
    index = index + 1
    
  }
}

COV<- function(x,y) {
  if(length(x)!=length(y)) {stop('x must have the same length as y ')}
  x.bar <- mean(x)
  y.bar <- mean(y)
  N <- length(x)
  
  Cov <- (sum((x-x.bar)*(y-y.bar))) / (N-1)
  return(Cov)
}


mean_exp_correlation_all_clones = function(){
  fnames = c("4.12 with mitosis.csv", "4.2 with mitosis.csv", 
             "4.11 with mitosis.csv", "4.10 with mitosis.csv", 
             "PGK 19 all.csv", "8.2 with mitosis.csv")
  clones = c("Dstn", "Nono", "Jarid2", "RbpJ", 'Pgk', 'Sproty4')
  par(mfrow=c(3,2), mar=c(5,5,5,2)) 
  for(i in 1:length(fnames)) {
    df = read.csv(fnames[i], header=FALSE)
    plot_mean_expression_correlation(df, clones[i])
  }
}

diff_related_cells = function(df, name=NULL) {
  require(RColorBrewer)
  reds = brewer.pal(9, "Reds")
  df.n = t(normalize_time(df))
  res = matrix(nr=dim(df.n)[1]/2, nc=dim(df.n)[2])
  for(i in seq(1, dim(df.n)[1]/2)) {
    res[i, ] = abs((df.n[(i-1)*2 + 1, ]) - (df.n[(i-1)*2 + 2, ]))
  }
  v = colMeans(res, na.rm = TRUE)
  s = apply(res, 2, function(x) ( sd(x, na.rm = TRUE)) )
  plot(v, ylim=c(0, .4), type='l', lwd=3, ylab="Absolute FC sister-1 / sister-2", xlab="Cell cycle percentage",
       main=paste("", name),
       cex.axis=1.5, cex.lab=1.6)
  polygon(c(1:100, 100:1), c(v-s/2, rev(v+s/2)), col=alpha("grey80", .3), border="grey80")
  lines(v, lwd=4)
  rand.res = matrix(0, nr=500, nc=100)
  for(index in 1:dim(rand.res)[1]) {
    rand.order = sample(seq(1, dim(df.n)[1]))
    res.rand = matrix(nr=dim(df.n)[1]/2, nc=dim(df.n)[2])
    for(i in seq(1, dim(df.n)[1]/2)) {
      res.rand[i, ] = abs(df.n[rand.order[(i-1)*2 + 1], ] - df.n[rand.order[(i-1)*2 + 2], ])
    }
    rand.res[index, ] = colMeans(res.rand, na.rm = TRUE)
  }
  v = colMeans(rand.res, na.rm = TRUE)
  s = apply(rand.res, 2, function(x) ( sd(x, na.rm = TRUE)) )
  polygon(c(1:100, 100:1), c(v-s/2, rev(v+s/2)), col=alpha(reds[4], .3), border=reds[3])
  lines(v, lwd=4, col=2)
}


diff_over_time = function(df, interval=10) {
  df.n = normalize_time(df)
  mean.diff = c()
  sd.diff = c()
  for (i in seq(1, 100 - interval/2 - 1, by=interval/2)) {
    m = matrix(colMeans(df.n[i:(i+interval-1), ], na.rm = TRUE), nc=2, byrow = TRUE)
    mean.diff = c(mean.diff, mean(abs(m[, 1] - m[, 2])))
    sd.diff = c(sd.diff, sd(abs(m[, 1] - m[, 2])))
  }
  plot(seq(1, 100 - interval/2 - 1, by=interval/2), mean.diff)
  lines(seq(1, 100 - interval/2 - 1, by=interval/2), mean.diff, lty=2)
}


fc_over_time = function(df, interval=10) {
  df.n = normalize_time(df)
  mean.diff = c()
  sd.diff = c()
  for (i in seq(1, 100 - interval/2 - 1, by=interval/2)) {
    m = matrix(colMeans(df.n[i:(i+interval-1), ], na.rm = TRUE), nc=2, byrow = TRUE)
    mean.diff = c(mean.diff, mean(abs(log2(m[, 1]) - log2(m[, 2]))))
    sd.diff = c(sd.diff, sd(abs(log2(m[, 1]) - log2(m[, 2]))))
  }
  plot(seq(1, 100 - interval/2 - 1, by=interval/2), mean.diff, ylim = c(0, max(mean.diff) + max(sd.diff)/2), 
       pch=16, cex=1.6, 
       cex.axis=2.4, cex.lab=1.7, xlab="% Cell cycle", ylab="absolute Log2 change")
  lines(seq(1, 100 - interval/2 - 1, by=interval/2), mean.diff, lty=2)
  x = seq(1, 100 - interval/2 - 1, by=interval/2)
  for (i in 1:length(sd.diff)) {
    segments(x[i], mean.diff[i] - sd.diff[i]/2, x[i], mean.diff[i] + sd.diff[i]/2, lwd=2)
    segments(x[i] - 1, mean.diff[i] - sd.diff[i]/2, x[i] + 1, mean.diff[i] - sd.diff[i]/2, lwd=2)
    segments(x[i] - 1, mean.diff[i] + sd.diff[i]/2, x[i] + 1, mean.diff[i] + sd.diff[i]/2, lwd=2)
  }
}

expression_corr_over_time = function(df, interval=10, name=NA) {
  df.n = normalize_time(df)
  cor.coefs = c()
  for (i in seq(1, 100 - interval/2 - 1, by=interval/2)) {
    m = matrix(colMeans(df.n[i:(i+interval-1), ], na.rm = TRUE), nc=2, byrow = TRUE)
    cor.coefs = c(cor.coefs, cor(m)[1, 2])
  }
  par(mar=c(7,7,3,2))
  plot(seq(1, 100 - interval/2 - 1, by=interval/2) + interval/2, cor.coefs, ylim = c(.2, 1), xlim=c(0,100), 
       col=0, xlab = "% Cell cycle", ylab="Correlation Coeff.", 
       cex.axis=2.4, cex.lab=1.7, main=name, cex.main=2, 
       bty='n')
  # mtext(side = 1, text = "% Cell cycle", line = 5, cex = 1.4)
  # mtext(side = 2, text = "Correlation coefficient", line = 5, cex = 1.4)
  
  abline(h=mean(cor.coefs), lty=2, col="grey", lwd=2)
  points(seq(1, 100 - interval/2 - 1, by=interval/2) + interval/2, cor.coefs, pch=16, cex=2)
  lines(seq(1, 100 - interval/2 - 1, by=interval/2) + interval/2, cor.coefs, lty=2, lwd=3)
}

expression_corr_over_time_all_clones2 = function(interval=10) {
  fnames = c("4.12 with mitosis.csv", "4.2 with mitosis.csv", 
             "4.11 with mitosis.csv", "4.10 with mitosis.csv", 
             "PGK 19 all.csv", "8.2 with mitosis.csv")
  clones = c("Dstn", "Nono", "Jarid2", "RbpJ", 'Pgk', 'Sproty4')
  par(mfrow=c(3,2), mar=c(5,5,5,2)) 
  
  for (i in 1:length(fnames)) {
    df = read.csv(fnames[i], header=FALSE)
    expression_corr_over_time(df, 10, clones[i])
  }
}


expression_corr_over_time_all_clones = function(interval=10) {
  fnames = c("4.12 with mitosis.csv", "4.2 with mitosis.csv", 
             "4.11 with mitosis.csv", "4.10 with mitosis.csv", 
             "PGK 19 all.csv", "8.2 with mitosis.csv")
  clones = c("Dstn", "Nono", "Jarid2", "RbpJ", 'Pgk', 'Sproty4')
  
  par(mar=c(7,7,3,2))
  plot(seq(1, 100 - interval/2 - 1, by=interval/2) + interval/2, seq(1, 100 - interval/2 - 1, by=interval/2), 
       ylim = c(0., 1), xlim=c(0,100), 
       col=0, 
       cex.axis=2.4, cex.lab=1.7, 
       bty='n', ann=FALSE, yaxt = 'n'
  )
  mtext(side = 1, text = "% Cell cycle", line = 4, cex = 2)
  mtext(side = 2, text = "Correlation coefficient", line = 5, cex = 2)
  axis(2, at=seq(0, 1, by=.2), labels = seq(0, 1, by=.2), cex.axis=2.4, las=2)
  require(RColorBrewer)
  cols = brewer.pal(7, "Set1") 
  require(scales)
  
  for(j in 1:length(fnames)) {
    df = read.csv(fnames[j], header=FALSE)
    df.n = normalize_time(df)
    cor.coefs = c()
    for (i in seq(1, 100 - interval/2 - 1, by=interval/2)) {
      m = matrix(colMeans(df.n[i:(i+interval-1), ], na.rm = TRUE), nc=2, byrow = TRUE)
      cor.coefs = c(cor.coefs, cor(m)[1, 2])
    }
    
    corr = cor(matrix(colMeans(df, na.rm = TRUE), nc=2, byrow = TRUE) )[1, 2]
    
    abline(h=corr, lty=2, col=alpha(cols[j], .65), lwd=2)
    points(seq(1, 100 - interval/2 - 1, by=interval/2) + interval/2, cor.coefs, pch=j+14, cex=2, 
           col=cols[j])
    lines(seq(1, 100 - interval/2 - 1, by=interval/2) + interval/2, cor.coefs, lty=1, lwd=3, col=cols[j])
  }
  
  legend("bottomleft", clones, col=cols[1:length(clones)], pch=1:length(clones)+14, lty=1, lwd=3, bty = 'n', 
         cex=1.2)
}


expression_fc_over_time = function(df, interval=10) {
  df.n = normalize_time(df)
  df.norm = sweep(df.n, 2, apply(df.n, 2, max), FUN='/')
  fc.mean = c()
  ncell = dim(df)[2]
  fc.sd = c()
  for (i in seq(1, 100 - interval/2 - 1, by=interval/2)) {
    m = log2(df.norm[i:(i+interval-1), ])
    fc = colMeans(abs(m[, seq(1, ncell-1, by=2)] - m[, seq(2, ncell, by=2)]))
    fc.mean = c(fc.mean, mean(fc))
    fc.sd = c(fc.sd, sd(fc))
  }
  par(mar=c(7,7,3,2))
  xaxis = seq(1, 100 - interval/2 - 1, by=interval/2) + interval/2
  plot(xaxis, fc.mean, ylim = c(0, max(fc.mean)+max(fc.sd)), xlim=c(0,100), 
       col=0, 
       cex.axis=2.4, cex.lab=1.7, 
       bty='n', ann=FALSE
  )
  
  mtext(side = 1, text = "% Cell cycle", line = 5, cex = 2)
  mtext(side = 2, text = "Correlation coefficient", line = 5, cex = 2)
  # abline(h=fc.mean, lty=2, col="grey", lwd=2)
  points(xaxis, fc.mean, pch=16, cex=2)
  lines(xaxis, fc.mean, lty=1, lwd=3)
  
  for(i in 1:length(fc.sd)) {
    segments( xaxis[i], fc.mean[i] - fc.sd[i] , xaxis[i],  fc.mean[i] + fc.sd[i], lwd=2)
    segments( xaxis[i]-.8, fc.mean[i] - fc.sd[i] , xaxis[i]+.8,  fc.mean[i] - fc.sd[i], lwd=2)
    segments( xaxis[i]-.8, fc.mean[i] + fc.sd[i] , xaxis[i]+.8,  fc.mean[i] + fc.sd[i], lwd=2)
  }
}


# pdf('exp_diff_expression_levels.pdf', height = 6, width = 9)
# t = seq(0, 99)
# par(mar=c(5,5,3,2))
# plot(t, df.n[, 1], ylim = range(df.n), type='l', lwd=2, col=cols[1], xlab="% Cell cycle", ylab="Luminescence signal",
#      cex.axis=1.7, cex.lab=1.8, bty='n')
# abline(h=mean(df.n[,2]), col=cols[1], lty=1) 
# lines(t, df.n[, 2], ylim = range(df.n), type='l', lwd=3, col=cols[1], lty=2)
# abline(h=mean(df.n[,2]), col=cols[1], lty=2)
# lines(t, df.n[, 78], ylim = range(df.n), type='l', lwd=3, col=cols[2], lty=1)
# abline(h=mean(df.n[,78]), col=cols[2], lty=1)
# lines(t, df.n[, 77], ylim = range(df.n), type='l', lwd=3, col=cols[2], lty=2)
# abline(h=mean(df.n[,77]), col=cols[2], lty=2)
# lines(t, df.n[, 41], ylim = range(df.n), type='l', lwd=3, col=cols[3], lty=1)
# abline(h=mean(df.n[,41]), col=cols[3], lty=1)
# lines(t, df.n[, 42], ylim = range(df.n), type='l', lwd=3, col=cols[3], lty=2)
# abline(h=mean(df.n[,42]), col=cols[3], lty=2)
# legend("topleft", c("sister 1", 'sister 2'), lty=c(1, 2), lwd=3, bty='n', cex=1.4)
# dev.off()

# par(mar=c(5,5,5,2))
# plot(t, df.n[, 1], ylim = range(df.n), type='l', lwd=2, col=cols[1], xlab="% Cell cycle", ylab="Luminescence signal",
#      cex.axis=1.7, cex.lab=1.8, bty='n')
# interval =10
# for (i in seq(1, 100 - interval/2 - 1, by=interval/2)) {
#   print(i)
#   m = mean(df.n[i:(i+interval-1), 1])
#   segments(i, m, i+interval, m, col=cols[1], lwd=1)
#   m2 = mean(df.n[(i + interval/2):(i+interval-1), 1])
# }
# 
# lines(t, df.n[, 78], ylim = range(df.n), type='l', lwd=3, col=cols[2], lty=2)
# for (i in seq(1, 100 - interval/2 - 1, by=interval/2)) {
#   print(i)
#   m = mean(df.n[i:(i+interval-2), 78])
#   segments(i, m, i+interval, m, col=cols[2], lwd=1)
#   m2 = mean(df.n[(i + interval/2):(i+interval-1), 1])
# }
# 

# df = read.csv('Mother-daughter pgk21.csv', header=FALSE)
# df.n = normalize_time(df)
# df.m = df.n[, seq(1, dim(df.n)[2], 4)]
# df.d1 = df.n[, seq(2, dim(df.n)[2], 4)]
# df.d2 = df.n[, seq(4, dim(df.n)[2], 4)]
# 
# m1 = matrix(c(colMeans(df.m), colMeans(df.d1)), byrow = FALSE, nc=2)
# m2 = matrix(c(colMeans(df.m), colMeans(df.d2)), byrow = FALSE, nc=2)

# pdf("mother_daughter_corr_coeff_mean_exp.pdf", height = 6, width = 6)
# par(mar=c(5,5,5,2))
# plot(m1, ylim=range(colMeans(df.n)), xlim=range(colMeans(df.n)), 
#      cex.axis=1.6, cex.lab=1.5, xlab="Mother cell", ylab = "Daughter cell", pch=16, cex=1.3, bty='n')
# abline(lm(m1[,2]~m1[,1]), lty=2, col=1, lwd=2)
# points(m2, col=cols[1], pch=17, cex=1.3)
# abline(lm(m2[,2]~m2[,1]), lty=2, col=cols[1], lwd=2)
# text(7000, 4000, paste("Correlation coeff.:", round(cor(m1)[1, 2], 2) ) , cex=1.4)
# text(7000, 3500, paste("Correlation coeff.:", round(cor(m2)[1, 2], 2) ) , cex=1.4, col=cols[1])
# legend('topleft', c("Daughter #1", "Daughter #2"), pch=c(16, 17), cex=1.1, col=c(1, cols[1]))
# dev.off()

# v1 = c()
# v2 = c()
# for(i in seq(1, dim(df.n)[2], 4)) {
#     v1 = c(v1, i)
#     v1 = c(v1, i + 1)
#     v2 = c(v2, i)
#     v2 = c(v2, i+3)
# }
# 
# df1 = df.n[, v1]
# df2 = df.n[, v2]
# 
# pdf("mother_daughter_corr_coeff_mean_exp_in_time.pdf", height = 6, width = 16)
# par(mfrow=c(1, 2))
# expression_corr_over_time(df1)
# expression_corr_over_time(df2)
# dev.off()
# 
# pdf("mother_daughter_fc.pdf", height = 6, width = 16)
# par(mfrow=c(1, 2), mar=c(5,5,2,2))
# fold_change_related_cells(df1)
# fold_change_related_cells(df2)
# dev.off()
# 
# pdf("mother_daughter_fc_mean_norm.pdf", height = 6, width = 16)
# par(mfrow=c(1, 2), mar=c(5,5,2,2))
# fold_change_normalized_for_mean(df1)
# fold_change_normalized_for_mean(df2)
# dev.off()

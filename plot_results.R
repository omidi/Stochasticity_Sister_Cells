folder = "results_4.11_simann"
# files = list.files(folder)
# ncells = length(files)
# cell_ids = sapply(files, function(x) as.integer(strsplit(x, "_")[[1]][2]) )
# 
# M = matrix(0, nr=ncells/2, nc=6)
# for(p in 1:(ncells/2)){
#    c1 = files[which(cell_ids==(p*2 - 1))]
#    c2 = files[which(cell_ids==p*2)]
#    df1 = read.table(paste(folder, c1, sep="/"), header=TRUE)
#    df2 = read.table(paste(folder, c2, sep="/"), header=TRUE)
#    M[p, 1:3] = unlist(df1[dim(df1)[1], 1:3])
#    M[p, 4:6] = unlist(df2[dim(df2)[1], 1:3])
# }
# 
# 
# colnames(M) = c("k_on1", "k_off1", "k_s1",
#                 "k_on2", "k_off2", "k_s2")
# 
# 
# M2 = M[which(M[, 3]>0.5), ]

cols =  c("k_ON", "k_OFF", "k_s")
# M = M[-4, ]
par(mfrow = c(2,2))

ind = seq(1, dim(M)[1])
rnd_ind = sample(ind)
for (i in 1:3) {
  x.range = range(M[, c(i, i+3)])
  plot(M[, i], M[, i+3], pch=19, xlim = x.range, ylim = x.range, 
       main = paste( cols[i], "corr:", round(cor(M[,i], M[, i+3]), 3)), 
       xlab="sister 1", ylab="sister 2")
  abline(lm(M[, i+3] ~ M[, i]))
  abline(0, 1, lty=2)
  # points(M[, i], M[rnd_ind, i+3], col=2, pch=2)
}

e1 = M[,3]*(M[, 1] / (M[, 1] + M[, 2])) / 0.004620981
e2 = M[,6]*(M[, 4] / (M[, 4] + M[, 5])) / 0.004620981
y.range = range(cbind(e1, e2))
plot(e1, e2, ylim=y.range, xlim=y.range, 
     main = paste("corr:", round(cor(e1,e2), 3)), pch=19, 
     xlab= "sister 1", ylab="sister 2")
abline(0, 1, lty=2)
abline(lm(e2 ~ e1))
# points(e1, e2[sample(ind)], pch=2, col=2)

nrand = 10000
rand_corr = matrix(0, nr=nrand, nc=4)
colnames(rand_corr) = c("k_ON", "k_OFF", "k_s", "E")
for (k in 1:3) {
  for(i in 1:nrand) {
    rand_corr[i,k] = cor(M[, k], M[sample(ind), k+3]) 
  }  
}


for(i in 1:nrand) {
  rand_corr[i,4] = cor(e1, e2[sample(ind)]) 
}  

# boxplot(rand_corr, ylim=c(-.6, .6))


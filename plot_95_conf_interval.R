# options(echo=TRUE) # if you want see commands in output file
# args <- commandArgs(trailingOnly = TRUE)
# dirname = args[1]

setwd("~/PycharmProjects/Stochasticity_Sister_Cells/")
dirname = "results_4.12_mcmc"
files = list.files(dirname, pattern = "mcmc_params$")

require(ellipse)
require(scales)

plot(1,1, ylim=c(0, 300), xlim=c(0, 1), col=0)
for(fname in files){
  d = read.table(paste(dirname, fname, sep='/'), header = TRUE)
  dd = subset(d[100:dim(d)[1], ], acc==1)
  x = cbind(dd[,1]/(dd[, 2] + dd[,1]), dd[, 3] / dd[, 1])
  lines(ellipse(cov(x), centre=colMeans(x), level=.95), lwd=1, col=alpha("grey", .4))
}

for(fname in files[19:20]){
  d = read.table(paste(dirname, fname, sep='/'), header = TRUE)
  dd = subset(d[100:dim(d)[1], ], acc==1)
  x = cbind(dd[,1]/(dd[, 2] + dd[,1]), dd[, 3] / dd[, 1])
  lines(ellipse(cov(x), centre=colMeans(x), level=.95), lwd=2, col=1)
}

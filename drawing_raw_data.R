setwd('~/PycharmProjects/Stochasticity_Sister_Cells/')
# df = read.csv('~/Dropbox/Communication/Gene.Expression.Inheritance/Data 29.04.16/4.11 with mitosis.csv', 
#               header=FALSE)
require(RColorBrewer)
cols = brewer.pal(7, "Set1")
dest = "4.11_raw_data_plots" 

# for (i in seq(1, dim(df)[2], by=2)) {
# 
#   if(dest!="") {
#     fname = sprintf('%s/%d_cell_pairs.pdf', dest, i)
#   } else {  
#     fname = sprintf('%d_cell_pairs.pdf', i)
#   }
#   
#   pdf(fname, height = 6, width = 8)
#   df2 = df[, c(i, i+1)]
#   colnames(df2) = c("V1", "V2")
#   y.range = range(df2, na.rm = TRUE)
#   df2 = subset(df2, !(is.na(V1) & is.na(V2)) )
#   times = seq(0, dim(df2)[1]*5 - 5, by=5 )
#   par(mar=c(5,5,5,2))
#   plot(times, df2$V1, col=cols[1], 
#        xlab = "time (hr)", ylab="Luminescence signal", 
#        cex.lab=1.5, cex.axis=1.4, xaxt='n', cex.main=2, 
#        main = i, ylim = y.range) 
#   axis(1, at=seq(0, max(times), by=120), labels=seq(0, max(times), by=120)/60 , cex.axis=1.5)
#   lines(smooth.spline(times[!is.na(df2$V1)], df2$V1[!is.na(df2$V1)], spar = .3), col=cols[1], lwd=3)
#   points(times, df2$V2, pch=2, col=cols[2])
#   lines(smooth.spline(times[!is.na(df2$V2)], df2$V2[!is.na(df2$V2)], spar = .3), col=cols[2], lwd=3)
#   dev.off()
# }
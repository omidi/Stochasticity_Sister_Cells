options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
fname = args[1]
column = args[2]

alpha = 27.63

require(scales)

df = read.table(fname, header=TRUE)
states = df[, "state"]
signal = df[, "signal"]
times = df[, "time"]
numb = df[, "p_n"]

res_file = paste(fname, ".pdf", sep="")
pdf(res_file, height = 6, width = 8)
par(mar=c(5,5,5,5))
plot(times, states, ylim = c(0, max(signal)), pch=15, cex=.5, col=0,
     xlab= "time (hr)", ylab="Luminescence signal",
     cex.lab=1.5, cex.axis=1.4, xaxt='n', cex.main=2,
     main = column)

axis(1, at=seq(0, max(times), by=120), labels=seq(0, max(times), by=120)/60 , cex.axis=1.5)
lines(smooth.spline(times, signal, spar = .4), type = 'l', lwd=3)
points(times, signal)

lines(times, numb*(alpha/2), col=2, lwd=3)
axis(4, at=seq(0, max(numb)+50, by=50)*(alpha/2), labels = seq(0, max(numb)+50, by=50),
     col.ticks=2, col.axis=2, cex.axis = 1.4, cex.lab=1.5, las=2)
mtext("Protein number", line=3.7, side=4, cex=1.5, las=0, col=2)

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
    polygon(c(start, start, end, end), c(-200, 100000, 100000, -200),
            col=alpha("grey", .3), border = alpha("grey", .3))
  }
  if(i >= length(times))  break
}

lines(times, numb*(alpha/2), col=2, lwd=3)

dev.off()



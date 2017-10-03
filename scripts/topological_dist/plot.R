#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
data_fn=args[1]
plot_fn=args[2]

data = read.table(data_fn)
tot = nrow(data)
sep1 = tot/3+1
sep2 = tot/3*2+1

if(length(grep('png', plot_fn))>0) {
	png(plot_fn, width=1000, height=500)
} else {
	pdf(plot_fn, width=12, height=6)
}
plot(data[,2], lwd=2, ylab="% Peaks", type='l')
abline(v=sep1, col='red', lwd=2, lty=2)
abline(v=sep2, col='red', lwd=2, lty=2)
dev.off()

# plot the Topological distribution of m6A peaks
# Zijun Zhang
# 9.18.2017
# revised 1.8.2018: better visualization by ggplot2

#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=T)
data_fn=args[1]
plot_fn=args[2]

data = read.table(data_fn)
tot = nrow(data)
sep1 = tot/3+1
sep2 = tot/3*2+1

# if(length(grep('png', plot_fn))>0) {
	# png(plot_fn, width=1000, height=500)
# } else {
	# pdf(plot_fn, width=12, height=6)
# }
#plot(data[,2], lwd=2, ylab="% Peaks", type='l')
#abline(v=sep1, col='red', lwd=2, lty=2)
#abline(v=sep2, col='red', lwd=2, lty=2)
#dev.off()

library(ggplot2)
colnames(data) = c('Peak_Tx', 'Peak_Tx_WithPeak', 'NumPeak')
data$coord = seq(1, nrow(data))

p = ggplot(data, aes(x=coord, y=Peak_Tx_WithPeak, fill=NumPeak, colour=NumPeak)) + geom_point(colour='darkblue', size=2) +
	#scale_colour_gradient(low = "orange",  high = "darkblue", breaks=c(0,10,50,100)) +
	geom_line() +
	guides(fill=F, colour=F) +
	theme_bw() +
	geom_vline(xintercept=c(sep1, sep2), linetype='dashed') +
	annotate('text', label=c("3'UTR", "CDS", "5'UTR"), 
		x=c( sep1/2, (sep1+sep2)/2, sep2+sep1/2 ), y=max(data[,2])*0.9)

#ggsave(file=plot_fn, plot=p, width=10, height=5, units='in')  # NOT WORKING in the Anaconda R. ZZJ 8.11.2018
png(plot_fn, width=500, height=250)
p
dev.off()
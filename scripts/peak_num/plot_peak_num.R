## Compute the overlapping status of peak files
## and plot a heatmap
## Zijun Zhang
## 10.23.2017
## revised 1.8.2018: better visualization
## revised 3.3.2018: generalize to all peak in a project

library(ggplot2)

argv = commandArgs(trailingOnly=T)
project = argv[1]
outfn = argv[2]
peak_idx = argv[3:length(argv)]
peak_idx = basename(peak_idx)
num_cmpr = length(peak_idx)

r1 = as.vector(sapply(peak_idx, function(x) rep(x, num_cmpr)))
r2 = rep(peak_idx, num_cmpr)
peak_df = data.frame(
	x=r1,
	y=r2,
	value=NA,
	stringsAsFactors=F
	)

for(i in 1:nrow(peak_df))
{
	#p1 = file.path('projects', project, 'clam', peak_df$x[i], 'narrow_peak.unique.bed')
	p1 = file.path('projects', project, 'macs2', peak_df$x[i], paste0(peak_df$x[i],'--nomodel_peaks.narrowPeak'))
	#p2 = file.path('projects', project, 'clam', peak_df$y[i], 'narrow_peak.unique.bed')
	p2 = file.path('projects', project, 'macs2', peak_df$y[i], paste0(peak_df$y[i],'--nomodel_peaks.narrowPeak'))
	cmd = paste0('bedtools intersect -a ', p1, ' -b ', p2, ' -u | wc -l')
	resp = system(cmd, intern=T)
	peak_df$value[i] = as.numeric(resp)
}


p = ggplot(data=peak_df, aes(x=x, y=y)) + 
	geom_tile(aes(fill=value)) + 
	scale_fill_gradient(low = "white",  high = "steelblue", guide=FALSE) +
	geom_text(aes(label = value), size=3) +
	xlab('Reference') + ylab('Target') +
	theme_bw() +
	theme(axis.text.x=element_text(size=9,angle = 90, hjust=0.5), axis.title=element_text(size=10), axis.text.y=element_text(size=9)) 
	
ggsave(outfn, plot=p, width=6, height=6)
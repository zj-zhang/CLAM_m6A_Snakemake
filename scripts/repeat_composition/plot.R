## Plot a pie chart for top repetitive elements
## Modified from "Dropbox/CLAM_m6A_meta_analysis/evaluator/evaluator_CLAM.R"
## Zijun Zhang
## 10.23.2017

library("ggplot2")

myTheme= theme(legend.text = element_text(size = 13),
						plot.title = element_text(size=13, face="bold"), axis.title.y = element_text(size=13), axis.title.x = element_text(size=13),
						axis.text.y = element_text(size=13, angle = 90, hjust = 0.5, vjust=0.5), 
						axis.text.x = element_text(size=13, angle=0, hjust=0.5, vjust=0.5),
						panel.border = element_rect(colour = "black", fill=NA, size=1), legend.background = element_rect(fill = "transparent", colour = "transparent"),
						text = element_text(size=13))

no_border = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.background = element_blank(), axis.line = element_line(colour = "black")) +
	theme(panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))						


argv = commandArgs(trailingOnly=T)

df = read.table(argv[1], header=F, stringsAsFactors=F)

colnames(df) = c('Type', 'Count')
top_types = df[order(-df[,2])[1:5],1]
top_types_idx = which(df[,1]%in%top_types)
dense_df = df[top_types_idx,]
dense_df = rbind.data.frame(dense_df, list(Type='others', Count=sum(df[-top_types_idx,2])))

percent = function(x, digits=1)
{
	paste(round(x*100, digits), sep='')
}


p = ggplot(data=dense_df, aes(x="", y=Count, fill=Type)) + geom_bar(stat='identity') +
		coord_polar("y", start=0) + theme(axis.text.x=element_blank()) +
		geom_text(aes(x=1.65, y = Count, label = Count), size=2, position=position_stack(vjust=0.5)) + 
		geom_text(aes(x=1.4, y = Count, label = percent(Count/sum(Count), 0)), size=2, position=position_stack(vjust=0.5)) + 
		scale_fill_brewer(palette="Dark2", name=NULL) +
		xlab('') + ylab('') +
		theme(legend.position='bottom', legend.text = element_text(size = 3.5), legend.key.size = unit(0.35,"cm"))

ggsave(argv[2], p, width=2.5, height=2.5)

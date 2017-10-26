# use R package to get all SRA filenames and FTP paths, given a GSE id.
# Zijun Zhang, 10.25.2016

library(GEOquery)

argv=commandArgs(trailingOnly=T)
max_try = 100
flag = 0
for(i in 1:max_try) {
	print(paste0("try ",i))
	gse1 <- try(getGEO(argv[1]))
	if(! inherits(gse1, "try-error") ) {
		flag=1
		break
	}
	Sys.sleep(30)
	
}
if(! flag) stop(paste0('i tried ', max_try, 'times and cannot get GEO info for ',argv[1]))
B = length(gse1)
write_out1 = NULL
for(b in 1:B) {
	df1 = as.data.frame(gse1[[b]])
	fn1 = unlist(lapply(df1$title, function(x) paste0(strsplit(as.character(x)," ")[[1]], collapse='_')))
	fn1 = gsub('/','_',fn1)
	fn1 = gsub('-','',fn1)
	fn1 = gsub('\\(','',fn1)
	fn1 = gsub('\\)','',fn1)
	all_path1 = apply(as.data.frame(df1[,colnames(df1)[grepl("supplementary",colnames(df1))]]), 2, as.character)
	ftp_path1 = c()
	for(i in 1:nrow(all_path1))
	{
		j = grepl("SRX", all_path1[i,])
		ftp_path1 = c(ftp_path1, all_path1[i,j])
	}
	write_out1 = rbind(write_out1, cbind(fn1, ftp_path1))
}
write.table(write_out1, file=argv[2], quote=F, row.names=F, col.names=F)

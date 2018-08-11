## read star output for mapping stats
## Zijun Zhang
## 8.11.2017
## revised 9.22.2017: incorporate into Snakemake
## revised 3.3.2018: update exon read criteria

import sys
import os
import subprocess


def read_star_log(fn):
	assert os.path.isfile(fn)
	stats = {}
	stats_list = [
		'Number of input reads',
		#'Uniquely mapped reads number',
		'Uniquely mapped reads %',
		#'Number of reads mapped to multiple loci',
		'% of reads mapped to multiple loci',
		'Number of splices: Total',
		]
	with open(fn, 'r') as f:
		for line in f:
			ele=[x.strip() for x in line.strip().split('|')]
			if ele[0] in stats_list:
				stats[ele[0]] = ele[1]
	return stats

def read_junction_reads(fn):
	assert os.path.isfile(fn)
	cmd = '''samtools view -q 10 %s | awk -F"\t" '$6~"N"&&$6!~"D"&&$6!~"I"&&$0~"NH:i:1"' | wc -l'''%fn
	sj=subprocess.check_output(cmd, shell=True)
	return sj.strip()

def read_exon_reads(fn, read_len=None):
	assert os.path.isfile(fn)
	assert read_len is None or type(read_len)==int
	if read_len is None:
		cmd2 = "samtools view -q 10 %s | head -n 1 | awk '{print length($10)}'"%fn
		print >>sys.stderr, cmd2
		read_len = subprocess.check_output(cmd2, shell=True)
		read_len = read_len.strip()
	#cmd = '''samtools view -q 10 %s | awk '$6=="%sM"&&$0~"NH:i:1"' | wc -l'''%(fn,read_len)
	cmd = '''samtools view -q 10 %s | awk '$0~/NH:i:1\>/&&$6~"M"&&$6!~"N"' | wc -l'''%(fn)
	#print >>sys.stderr, cmd
	sj=subprocess.check_output(cmd, shell=True)
	return sj.strip()

def read_barcode(fn):
	bar_dict = {}
	with open(fn, 'r') as f:
		firstline=True
		for line in f:
			ele = line.strip().split()
			if firstline:
				header = {ele[x]:x for x in range(len(ele))}
				firstline=False
				continue
			bar_dict[ele[header['barcode']]] = ele[1:]
	return bar_dict


def mapping_stat(bam_fn, verbose=True):
	stats_list = [
		'Number of input reads',
		#'Uniquely mapped reads number',
		'Uniquely mapped reads %',
		#'Number of reads mapped to multiple loci',
		'% of reads mapped to multiple loci',
		'Number of splices: Total',
		'Number of splice junction reads',
		'Number of exon reads'
		]
	#if verbose: print '\t'.join(['Dataset']+stats_list)
	source_dir = os.path.dirname(bam_fn)		
	stat_fn = os.path.join(source_dir, 'Log.final.out')
	stats = read_star_log(stat_fn)
	sj = read_junction_reads(bam_fn)
	#ec = str(int(stats['Uniquely mapped reads number'])*2 - int(sj))
	ec = read_exon_reads(bam_fn)
	stats['Number of splice junction reads'] = sj
	stats['Number of exon reads'] = ec
	if verbose:
		for x in stats_list:
			print x + '\t' + str(stats[x])
	return stats
		
		
if __name__ == '__main__':
	if len(sys.argv)>0:
		bam_fn = sys.argv[1]
		mapping_stat(bam_fn)
	else:
		sys.exit(1)


## make the [hg19|mm10]_ensg_genename.txt file
## given a ensGene.txt and a GTF.
## Zijun Zhang
## 10.2.2017

import sys
import re


def read_gtf_tx(fn):
	tx2type = {}
	with open(fn, 'r') as f:
		for line in f:
			if line.startswith('#'):
				continue
			ele=line.strip().split()
			if ele[2]=='transcript':
				tx_id = re.search(r'transcript_id "(.+?)";',line).group(1).split('.')[0]
				tx_type = re.search(r'transcript_type "(.+?)";', line).group(1)
				tx2type[tx_id] = tx_type
	return tx2type



tx2type=read_gtf_tx('/u/nobackup/yxing/NOBACKUP/frankwoe/mm10/gencode.vM13.annotation.gtf')


for line in sys.stdin:
	if line.startswith('#'):
		print >> sys.stdout, "#name\tsource"
		continue
	tx = line.strip().split('\t')[1]
	if tx in tx2type:
		print >> sys.stdout, "{0}\t{1}".format(tx, tx2type[tx])
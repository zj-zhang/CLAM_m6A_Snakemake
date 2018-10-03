'''Use pyBigwig to merge junction counts
from different replicates
Zijun Zhang
June 20, 2018
'''

import sys
import os
import pyBigWig
from collections import defaultdict
from pathlib import Path


outfn = sys.argv[1]
group = sys.argv[2]
genome = sys.argv[3]
in_bb_list = sys.argv[4:]

## DEBUG; COMMET OUT BEFORE USE
#in_bb_list = [os.path.join('Gastro_RCNP_FGC1866/bigwig',x,x+'.junc.bb') for x in os.listdir('Gastro_RCNP_FGC1866/bigwig/')]

# open file handlers in a list
bb_filehandler_list = [pyBigWig.open(x) for x in in_bb_list]

# read in chrom sizes
chrom_size_fn = os.path.join(str(Path(__file__).parents[1]),"UCSC", genome+'.chrom.sizes')
chrom_sizes = [x.strip().split() for x in open(chrom_size_fn, 'r').readlines()]
chrom_sizes = {x[0]:int(x[1]) for x in chrom_sizes}

# merging
MIN_JC_COUNT = 2
with open(outfn+'.tmp', 'w') as fo:
	for chrom in chrom_sizes:
		if '_' in chrom: # skip random chunks
			continue
		# junction dict for this chromosome
		chrom_jc_count = defaultdict(int)
		chrom_jc_info = {}
		for bb in bb_filehandler_list:
				entries = bb.entries(chrom, 1, chrom_sizes[chrom])
				for this in entries:
					id = '\t'.join([chrom, str(this[0]), str(this[1]) ])
					count = int(this[2].split('\t')[0].split('_')[-1])
					info = '\t'.join(this[2].split('\t')[1:])
					chrom_jc_count[id] += count
					chrom_jc_info[id] = info

		## write out this chrom
		for id in chrom_jc_info:
			jc_count = chrom_jc_count[id]
			if jc_count <= MIN_JC_COUNT:
					continue
			jc_info = chrom_jc_info[id]
			fo.write('\t'.join([id, group+'_'+str(jc_count), jc_info ]) + '\n')

# convert to BB
os.system("LC_COLLATE=C sort -k1,1 -k2,2n {} > {}".format(outfn+'.tmp', outfn+'.sorted.tmp') )

convert_bin = os.path.join( str(Path(__file__).parents[1]), 'UCSC', 'bedToBigBed' )
os.system("{bin} {in_bed} {chrom_size} {out_bb}".format(bin=convert_bin, in_bed=outfn+'.sorted.tmp', chrom_size=chrom_size_fn, out_bb=outfn) )

# clean up
os.remove(outfn+'.tmp')
os.remove(outfn+'.sorted.tmp')
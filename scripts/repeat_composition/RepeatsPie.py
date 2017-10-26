#!python2
## Count repetitive elements compotision for a given set of peaks
## Modified from analysis of NAR paper
## Zijun Zhang, 2017.10.23


import sys
import os
import subprocess
from collections import defaultdict

if __name__ == '__main__':
	cur_dir = os.path.dirname(os.path.realpath(__file__))
	pean_fn = sys.argv[1]
	repeat_fn = os.path.join( cur_dir, sys.argv[2] + '_repeatMasker.sorted.bed' )
	
	total_peak = subprocess.check_output('wc -l '+pean_fn, shell=True).split()[0]
	total_peak = int(total_peak)
	
	intersection = subprocess.check_output('bedtools intersect -a ' + repeat_fn +' -b ' + pean_fn + ' -wo', shell=True)
	
	repeat_count = defaultdict(int)
	in_repeat=set()
	for line in intersection.split('\n'):
		try:
			ele=line.split('\t')
			family = ele[4] if ele[5]==ele[11] else 'anti-'+ele[4]
			repeat_count[family] += 1
			#in_repeat.add(ele[-4])
			in_repeat.add(':'.join([ele[6],ele[7],ele[8]]))
		except:
			print >> sys.stderr, line
	
	#repeat_sum = 0
	for repeat in repeat_count:
		print repeat + '\t' + str(repeat_count[repeat])
		#repeat_sum += repeat_count[repeat]
	#print 'Non-repeat\t' + str(total_peak - repeat_sum)
	print 'Non-repeat\t' + str(total_peak - len(in_repeat))

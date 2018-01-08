## convert CLAM output into bigwig
## modified from CLAM evaluator back in 2016
## Zijun Zhang
## 11.6.2017

import sys
import subprocess
import os

# get auxillary scripts path
cur_dir = os.path.dirname(os.path.realpath(__file__))
# bam_script: clam_bam2bed.py aligner_out is_stranded +/-
bam_script = os.path.join(cur_dir, 'clam_bam2bed.py')
# bam_script: clam_bam2bed.py score_column
bed_script = os.path.join(cur_dir, 'clam_bed2bedG.py')

# parse command-line args
unique_bam = sys.argv[1]
aligner_out = sys.argv[2]
out_bg_dir = sys.argv[3]
is_stranded = True if len(sys.argv)>=5 and sys.argv[4].lower().startswith('t') else False
verbose = True

# make unique bedGraph
if is_stranded:
	#cmd = 'bedtools genomecov -bg -split -strand + -ibam ' + unique_bam + ' -g /u/home/f/frankwoe/nobackup/programs/UCSC/hg19.chrom.sizes | LC_COLLATE=C sort -k1,1 -k2,2n > ' + out_bg_dir + '/unique_pos.bdg'
	cmd = 'python2 %s %s %s %s %s'%(bam_script, unique_bam, is_stranded, '+', 'False')+''' | awk '{print $1"\t"$2"\t"$3"\t"$5}' ''' + '| LC_COLLATE=C sort -k1,1 -k2,2n | python2 ' + bed_script + ' 4 > ' + out_bg_dir + '/unique_pos.bdg'
else:
	#cmd = 'bedtools genomecov -bg -split -ibam ' + unique_bam + ' -g /u/home/f/frankwoe/nobackup/programs/UCSC/hg19.chrom.sizes | LC_COLLATE=C sort -k1,1 -k2,2n > ' + out_bg_dir + '/unique_pos.bdg'
	cmd = 'python2 %s %s %s %s %s'%(bam_script, unique_bam, is_stranded, '+', 'False')+''' | awk '{print $1"\t"$2"\t"$3"\t"$5}' ''' + '| LC_COLLATE=C sort -k1,1 -k2,2n | python2 ' + bed_script + ' 4 > ' + out_bg_dir + '/unique_pos.bdg'
if verbose:
	print cmd
subprocess.call(cmd, shell=True)
#cmd = 'bedtools genomecov -bg -split -strand - -ibam ' + unique_bam + ' -g /u/home/f/frankwoe/nobackup/programs/UCSC/hg19.chrom.sizes | LC_COLLATE=C sort -k1,1 -k2,2n > ' + out_bg_dir + '/unique_neg.bdg'
cmd = 'python2 %s %s %s %s %s'%(bam_script, unique_bam, is_stranded, '-', 'False')+''' | awk '{print $1"\t"$2"\t"$3"\t"$5}' ''' + '| LC_COLLATE=C sort -k1,1 -k2,2n | python2 ' + bed_script + ' 4 > ' + out_bg_dir + '/unique_neg.bdg'
if is_stranded:
	if verbose: print cmd
	subprocess.call(cmd, shell=True)

# make multiple bedGraph
if is_stranded:
	cmd = 'python2 %s %s %s %s'%(bam_script, aligner_out, is_stranded, '+')+''' | awk '{print $1"\t"$2"\t"$3"\t"$5}' ''' + '| LC_COLLATE=C sort -k1,1 -k2,2n | python2 ' + bed_script + ' 4 > ' + out_bg_dir + '/multi_pos.bdg'
	if verbose: print cmd
	subprocess.call(cmd, shell=True)

	cmd = 'python2 %s %s %s %s'%(bam_script, aligner_out, is_stranded, '-')+''' | awk '{print $1"\t"$2"\t"$3"\t"$5}' ''' + '| LC_COLLATE=C sort -k1,1 -k2,2n | python2 ' + bed_script + ' 4 > ' + out_bg_dir + '/multi_neg.bdg'
	if verbose: print cmd
	subprocess.call(cmd, shell=True)
else:
	cmd = 'python2 %s %s %s %s'%(bam_script, aligner_out, is_stranded, '+')+''' | awk '{print $1"\t"$2"\t"$3"\t"$5}' ''' + '| LC_COLLATE=C sort -k1,1 -k2,2n | python2 ' + bed_script + ' 4 > ' + out_bg_dir + '/multi_pos.bdg'
	if verbose: print cmd
	subprocess.call(cmd, shell=True)

# take union to get combined bedGraph

cmd = 'bedtools unionbedg -i ' + out_bg_dir + '/multi_pos.bdg' + ' ' + out_bg_dir + '''/unique_pos.bdg | awk '{print $1"\t"$2"\t"$3"\t"$4+$5}' >''' + out_bg_dir + '/combined_pos.bdg'
if verbose: print cmd
subprocess.call(cmd, shell=True)

cmd = 'bedtools unionbedg -i ' + out_bg_dir + '/multi_neg.bdg' + ' ' + out_bg_dir + '''/unique_neg.bdg | awk '{print $1"\t"$2"\t"$3"\t"$4+$5}' >''' + out_bg_dir + '/combined_neg.bdg'
if is_stranded:
	if verbose: print cmd
	subprocess.call(cmd, shell=True)

# convert combined bedGraph to BigWig
cmd = 'bedGraphToBigWig ' + out_bg_dir + '/combined_pos.bdg /u/home/f/frankwoe/nobackup/programs/UCSC/hg19.chrom.sizes ' + out_bg_dir + '/combined_pos.bw'
if verbose: print cmd
subprocess.call(cmd, shell=True)

cmd = 'bedGraphToBigWig ' + out_bg_dir + '/combined_neg.bdg /u/home/f/frankwoe/nobackup/programs/UCSC/hg19.chrom.sizes ' + out_bg_dir + '/combined_neg.bw'
if is_stranded:
	if verbose: print cmd
	subprocess.call(cmd, shell=True)

# convert unique bedGraph to BigWig
cmd = 'bedGraphToBigWig ' + out_bg_dir + '/unique_pos.bdg /u/home/f/frankwoe/nobackup/programs/UCSC/hg19.chrom.sizes ' + out_bg_dir + '/unique_pos.bw'
subprocess.call(cmd, shell=True)

cmd = 'bedGraphToBigWig ' + out_bg_dir + '/unique_neg.bdg /u/home/f/frankwoe/nobackup/programs/UCSC/hg19.chrom.sizes ' + out_bg_dir + '/unique_neg.bw'
if is_stranded:
	subprocess.call(cmd, shell=True)

# clean up
os.remove(os.path.join(out_bg_dir, "unique_pos.bdg"))
os.remove(os.path.join(out_bg_dir, "multi_pos.bdg"))
os.remove(os.path.join(out_bg_dir, "combined_pos.bdg"))
if is_stranded:
	os.remove(os.path.join(out_bg_dir, "unique_neg.bdg"))
	os.remove(os.path.join(out_bg_dir, "multi_neg.bdg"))
	os.remove(os.path.join(out_bg_dir, "combined_neg.bdg"))

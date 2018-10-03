'''
Takes a list of bigWig files, and merge into one
Zijun Zhang
June 19, 2018
'''

import sys
import os

out_fn = sys.argv[1]
genome = sys.argv[2]
bw_list = sys.argv[3:]
script_dir = os.path.dirname(__file__)
bwMerge_bin = os.path.join(os.path.dirname(script_dir), 'UCSC', 'bigWigMerge')
cmd = "{bin} {in_list} {out_bedG}".format(bin=bwMerge_bin, in_list=' '.join(bw_list), out_bedG=out_fn+'.tmp')
print(cmd)
os.system(cmd)

bedGtoBw = os.path.join(os.path.dirname(script_dir), 'UCSC', 'bedGraphToBigWig')
genom_size = os.path.join(os.path.dirname(script_dir), 'UCSC', genome+'.chrom.sizes')
os.system("LC_COLLATE=C sort -k1,1 -k2,2n {} > {}".format(out_fn+'.tmp', out_fn+'.sorted.tmp'))
cmd2 = "{bin} {tmp} {genom_size} {out}".format(bin=bedGtoBw, tmp=out_fn+'.sorted.tmp', genom_size=genom_size, out=out_fn)
print(cmd2)
os.system(cmd2)

os.remove(out_fn+'.tmp')
os.remove(out_fn+'.sorted.tmp')
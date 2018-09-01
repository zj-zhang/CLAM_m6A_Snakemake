'''
Make softlink utils
for current Snakefile
ZZJ
Doc 8.10.2018
Usage:
python2 make_ln_2.py /path/to/source_dir /path/to/target_dir [suffix]
'''

import os
import sys

ref_dir = os.path.realpath(sys.argv[1])
tar_dir = os.path.realpath(sys.argv[2])
suffix = sys.argv[3] if len(sys.argv)>3 else ''

print "ref_dir="+ref_dir
print "tar_dir="+tar_dir

for fn in os.listdir(ref_dir):
    if not ( fn.endswith('fastq.gz') or fn.endswith('.fastq')):
        continue
    if 'IP' in fn.upper():
        type='ip'
    elif 'INP' in fn.upper() or 'CON' in fn.upper():
        type='con'
    else:
        continue
    src = os.path.join(ref_dir, fn)
    new_fn = fn.replace("_1","").replace("-","_")
    #tar = os.path.join(tar_dir, type, new_fn)
    new_fn_prefix = os.path.basename(new_fn).split('.')[0]
    foo_fn = os.path.join(tar_dir, new_fn_prefix, 'foo.txt')
    sub_dir = os.path.join(tar_dir, new_fn_prefix)
    new_fn = new_fn.split('.')[0] + suffix +  '.' +  '.'.join(new_fn.split('.')[1:])
    tar = os.path.join(tar_dir, new_fn_prefix, new_fn)
    #print sub_dir
    if not os.path.isdir(sub_dir):
    	os.makedirs(sub_dir)
    cmd = "ln -sf %s %s"%(src, tar)
    print cmd
    os.system(cmd)
    os.system('touch %s'%foo_fn)

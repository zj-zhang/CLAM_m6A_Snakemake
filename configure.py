#!/usr/bin/env python2

""" Make configuration files for a new project in Snakemake pipeline
Zijun Zhang
10.22.2017
"""

import os
import sys
import re
from collections import defaultdict

#sra_script = 'scripts/geo_downloader/getSRApath.R'
sra_script = 'scripts/geo_downloader/getSRApath.py'


def get_sra_path(gid):
	r = re.compile('GSE\d\d\d\d\d')
	if r.match(gid) is None or len(gid)!=8:
		raise Exception('wrong gid format')
	outfn = os.path.join('projects', gid, 'config', 'sra_dir.txt')
	if not os.path.isdir(os.path.join('projects', gid)):
		os.mkdir(os.path.join('projects', gid))
	if not os.path.isdir(os.path.join('projects', gid, 'config')):
		os.mkdir(os.path.join('projects', gid, 'config'))
	if not os.path.isfile(outfn):
		print "\n------\nwait, connecting to GEO..\n"
		#os.system("Rscript {script} {gid} {outfn} >/dev/null 2>&1".format(script=sra_script, gid=gid, outfn=outfn))
		os.system("python2 {script} {gid} {outfn}".format(script=sra_script, gid=gid, outfn=outfn))
	
	sample_list = []
	with open(outfn, 'r') as f:
		for line in f:
			ele = line.strip().split()
			sample_list.append(ele[0])
	print sample_list
	while True:
		resp = raw_input('\n------\ndo you want to remove samples? y/[n]\n>')
		if resp=='y':
			new_sample_list = []
			for sample in sample_list:
				while True:
					resp2 = raw_input('%s [y]/n:'%sample)
					if resp2=='n':
						break
					else:
						new_sample_list.append(sample)
						break
			sample_list=new_sample_list
			break
		else:
			break
	return sample_list


def parse_sample_type_dict(sample_list, guess=True):
	s = ''
	sample_type_dict = defaultdict(list)
	if guess:
		for sample in sample_list:
			#if 'input' in sample.lower() or 'control' in sample.lower() or 'rna' in sample.lower():
			if 'input' in sample.lower() or 'control' in sample.lower():
				s += '    {0}: {1}\n'.format(sample, 'con')
				sample_type_dict['con'].append(sample)
			else:
				s += '    {0}: {1}\n'.format(sample, 'ip')
				sample_type_dict['ip'].append(sample)
	
	print "\n------\ni guessed the following:"
	print s
	while(True):
		resp = raw_input('is it ok? [y]/n\n>')
		if resp!='n':
			s = 'sample_type_dict:\n' + s
			return s, sample_type_dict
		else:
			break
	sample_type_dict['ip'] = []
	sample_type_dict['con'] = []
	s = 'sample_type_dict:\n'
	for sample in sample_list:
		while True:
			resp = raw_input('%s [ip/con]'%sample)
			if resp=='ip' or resp=='con':
				break
		s +='    {0}: {1}\n'.format(sample, resp)
		sample_type_dict[resp].append(sample)
	
	return s, sample_type_dict


def parse_sample_comparison(sample_dict, sample_type_dict, guess=True):
	s = ''
	num_ip = len(sample_type_dict['ip'])
	num_con = len(sample_type_dict['con'])
	ip_list = sorted(sample_type_dict['ip'])
	con_list = sorted(sample_type_dict['con'])
	if guess:
		min_len = min(len(ip_list), len(con_list))
		for i in range(min_len):
			x1 = i
			x2 = i
			comp = '[ ' + ip_list[x1] + ', ' + con_list[x2] + ' ],'
			s += '        {0}\n'.format(comp)
		s += '    ]\n'
		s = 'sample_comparison:\n    [\n' + s
		print "i guessed the following:\n"
		print s
		resp = raw_input('is it ok? [y]/n\n>')
		if resp!='n':
			return s
	s = ''	
	while True:
		print "ip:\n"+ ', '.join([str(i)+':' + ip_list[i] for i in range(num_ip) ])
		print "con:\n"+ ', '.join([str(i)+':'+ con_list[i] for i in range(num_con) ])
		
		try:
			x1 = int(raw_input('\n------\nselect ip number [0-'+str(num_ip-1)+']: '))
		except ValueError:
			break
		assert 0 <= x1 < num_ip
		try:
			x2 = int(raw_input('select con number [0-'+str(num_con-1)+']: '))
		except ValueError:
			break
		assert 0 <= x2 < num_con
		
		comp = '[ ' + ip_list[x1] + ', ' + con_list[x2] + ' ],'
		print "\n------\nabout to add comparison: \n" + comp
		resp = raw_input('do you confirm [y]/n/e/r, e=yes and end input; r=start over\n>')
		if resp=='y' or resp=='':
			s += '        {0}\n'.format(comp)
		elif resp=='e':
			s += '        {0}\n'.format(comp)
			break
		elif resp=='r':
			s = ''
		print '\n------\ncurrent comparisons:'
		print s
	s += '    ]\n'
	s = 'sample_comparison:\n    [\n' + s
	return s
	

def parse_template():
	s = """
clam:
    max_tags: -1

genome: hg19

hg19:
    star_idx: /u/nobackup/yxing/NOBACKUP/frankwoe/hg19/star_idx_gencode_v19
    gtf: /u/nobackup/yxing/NOBACKUP/frankwoe/hg19/gencode.v19.annotation.gtf
    kallisto_idx: /u/nobackup/yxing/NOBACKUP/frankwoe/hg19/kallisto_pctx_idx_gencodeV19/hg19_pc_tx_kal.idx
mm10:
    star_idx: /u/nobackup/yxing/NOBACKUP/frankwoe/mm10/star_idx_gencodeM13
    gtf: /u/nobackup/yxing/NOBACKUP/frankwoe/mm10/gencode.vM13.annotation.gtf
    kallisto_idx: /u/nobackup/yxing/NOBACKUP/frankwoe/mm10/mm10_kallisto_pctx_idx/mm10_pc_tx_kal_idx
"""
	return s
	
def configure(gid):
	s = ''
	# get sample list by querying GEO
	print "***--- get sample list by querying GEO "+gid+" ---***"
	print "***--- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="+gid+" ---***\n\n"
	sample_list = get_sra_path(gid)
	# split ip and con
	print "***--- split ip and con ---***\n\n"
	s1, sample_type_dict = parse_sample_type_dict(sample_list)
	print sample_type_dict
	s += s1
	# enter comparisons
	print "***--- enter comparisons ---***\n\n"
	s2 = parse_sample_comparison(sample_list, sample_type_dict)
	s += s2
	# add template
	s3 = parse_template()
	s += s3
	# write to file
	print "***--- a final look before writing to disk ---***"
	print "***--- gid: "+gid+' ---***\n\n'
	#print "\n------\na final look before writing to disk.\n------"
	print s
	raw_input('\n------\n press enter to continue')
	config_fp = os.path.join('projects', gid, 'config', 'config.yaml')
	with open(config_fp, 'w') as f:
		f.write(s)
	os.mkdir(os.path.join('projects', gid, 'reads'))
	for x in sample_list:
		os.mkdir(os.path.join('projects', gid, 'reads', x) )

	
if __name__ == '__main__':
	gid = sys.argv[1]
	configure(gid)

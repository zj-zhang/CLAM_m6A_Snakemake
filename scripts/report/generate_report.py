#!/usr/bin/env python2

"""Generate html report for m6A CLAM pipeline
"""

import bs4 as bs
import os
import sys
import pandas as pd



def read_homer_table(fn):
	par_dir = os.path.dirname(os.path.realpath(fn))
	with open(fn, 'r') as f:
		data = ''.join(f.readlines())
		soup = bs.BeautifulSoup(data, 'lxml')
		table = soup.find('table')
		homer_table = str(table).replace('homerResults', os.path.join(par_dir,'homerResults'))
	return homer_table


def count_peak_num(fn):
	peak_counter = 0
	with open(fn, 'r') as f:
		for line in f:
			peak_counter += 1
	return peak_counter


def read_star_mapping_stats(fn):
	stats = {}
	stats_list = [
		'Number of input reads',
		'Uniquely mapped reads number',
		'Uniquely mapped reads %',
		'Number of splices: Total',
		]
	with open(fn, 'r') as f:
		for line in f:
			ele=[x.strip() for x in line.strip().split('|')]
			if ele[0] in stats_list:
				stats[ele[0]] = ele[1]
	return stats


def read_project_mapping_stats(ip_sample_list, con_sample_list, project_name):
	stats_df = pd.DataFrame(columns=['Sample', 'Type', 'Number of input reads',
		'Uniquely mapped reads number',
		'Uniquely mapped reads %',
		'Number of splices: Total',
		'Number of splice junction reads',
		'Number of exon reads'])
	tmp = []
	for ip_sample in ip_sample_list:
		this_dict = {'Sample': ip_sample, 'Type':'IP'}
		with open(os.path.join('star', project_name, 'ip', ip_sample, 'mapping_stats.txt'), 'r') as f:
			for line in f:
				try:
					key, val = line.strip().split('\t')
				except:
					continue
				this_dict[key] = val
		tmp.append(this_dict)
	for con_sample in con_sample_list:
		this_dict = {'Sample': con_sample, 'Type':'Con'}
		with open(os.path.join('star', project_name, 'con', con_sample, 'mapping_stats.txt'), 'r') as f:
			for line in f:
				try:
					key, val = line.strip().split('\t')
				except:
					continue
				this_dict[key] = val
		tmp.append(this_dict)
	stats_df = stats_df.append(tmp, ignore_index=True)
	return stats_df.to_html(border=1)


def generate_report(ip_sample_list, con_sample_list, pardir, outfn, project_name):
	html_str = '<html><head><title>Evaluation of Project "%s"</title></head>\n'%project_name
	html_str += '<body>\n'
	html_str += '<h1>%s</h1>\n'%project_name
	
	#*--- MAPPING STATS FOR EACH BAM FILE ---*
	html_str += '<hr>\n'
	html_str += '<h2>Mapping Stats</h2>\n'
	html_str += read_project_mapping_stats(ip_sample_list, con_sample_list, project_name)
	
	#stats_list = [
	#	'Number of input reads',
	#	'Uniquely mapped reads number',
	#	'Uniquely mapped reads %',
	#	'Number of splices: Total',
	#	]
	#*--- EVALUATION FOR EACH IP-CON COMPARISON ---*	
	for ip_sample in ip_sample_list:
		for con_sample in con_sample_list:
			html_str += '<hr>\n'
			comparison = "{ip_sample}-{con_sample}".format(ip_sample=ip_sample, con_sample=con_sample)
			homer_table = read_homer_table( os.path.join(pardir, 'homer', project_name, comparison, 'homerResults.html') )
			topology_img = os.path.join(pardir, 'topology', project_name, comparison, 'dist.png' )
			#ip_mapping_stat = read_star_mapping_stats( os.path.join(pardir, 'star', project_name, 'ip', ip_sample, 'Log.final.out') )
			#con_mapping_stat = read_star_mapping_stats( os.path.join(pardir, 'star', project_name, 'con', con_sample, 'Log.final.out') )
			tot_peak = count_peak_num(os.path.join(pardir, 'clam', project_name, 'peaks-'+comparison, 'narrow_peak.unique.bed'))
			
			html_str += '<h2>%s</h2>\n'%comparison
			# summary stats
			html_str += '<h3>Summary statistics</h3>\n'
			#html_str += '<h5>IP Mapping</h5>\n'
			#for cat in stats_list:
			#	html_str += '<p>%s = %s</p>\n'%( cat, ip_mapping_stat[cat] )
			#html_str += '<h5>Input Mapping</h5>\n'
			#for cat in stats_list:
			#	html_str += '<p>%s = %s</p>\n'%( cat, con_mapping_stat[cat] )
			html_str += '<h5>Peak calling</h5>\n'
			html_str += '<p>No. peaks = %s</p>\n'%tot_peak
			html_str += '<br>\n'
			# homer motifs
			html_str += '<h3>HOMER motifs</h3>\n'
			html_str += homer_table
			html_str += '<br>\n'
			# topology distribution
			html_str += '<h3>Topology distribution</h3>\n'
			html_str += '<img src="%s" />\n'%topology_img
	html_str += '</body>\n'
	html_str += '</html>\n'
	with open(outfn, 'w') as f:
		f.write(html_str)
	return

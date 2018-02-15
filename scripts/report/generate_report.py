#!/usr/bin/env python2

"""Generate html report for m6A CLAM pipeline
"""

import bs4 as bs
from lxml import html
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
	html_table = html.fragment_fromstring(homer_table)
	top_row = 5
	row_counter = 0
	for row in html_table.iterchildren():
		row_counter += 1
		row.remove(row.getchildren()[-1])
		row.remove(row.getchildren()[-1])
		if row_counter>= top_row:
			row.clear()
	
	html_table = str(html.tostring(html_table, encoding='unicode', with_tail=False))
	return html_table


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
	target_col = ['Sample', 'Type',
		'Number of input reads',
		#'Uniquely mapped reads number',
		'Uniquely mapped reads %',
		#'Number of reads mapped to multiple loci'，
		'% of reads mapped to multiple loci',
		'Number of splices: Total',
		'Number of splice junction reads',
		'Number of exon reads'
		]
	stats_df = pd.DataFrame(columns=target_col)
	tmp = []
	for ip_sample in ip_sample_list:
		this_dict = {'Sample': ip_sample, 'Type':'IP'}
		with open(os.path.join('projects', project_name, 'star', 'ip', ip_sample, 'mapping_stats.txt'), 'r') as f:
			for line in f:
				try:
					key, val = line.strip().split('\t')
				except:
					continue
				this_dict[key] = val
		tmp.append(this_dict)
	for con_sample in con_sample_list:
		this_dict = {'Sample': con_sample, 'Type':'Con'}
		with open(os.path.join('projects', project_name, 'star', 'con', con_sample, 'mapping_stats.txt'), 'r') as f:
			for line in f:
				try:
					key, val = line.strip().split('\t')
				except:
					continue
				this_dict[key] = val
		tmp.append(this_dict)
	stats_df = stats_df.append(tmp, ignore_index=True)
	stats_df = stats_df[target_col]
	stats_df.index = stats_df['Sample']
	#return stats_df.to_html(border=1)
	return stats_df


def generate_report(comparison_list, pardir, outfn, project_name, include_mread_analysis=True):
	ip_sample_list = list(set([x[0] for x in comparison_list]))
	con_sample_list = list(set([x[1] for x in comparison_list]))
	
	html_str = '<html><head><title>Evaluation of Project "%s"</title></head>\n'%project_name
	html_str += '<body>\n'
	html_str += '<h1>%s</h1>\n'%project_name
	
	#*--- MAPPING STATS FOR EACH BAM FILE ---*
	html_str += '<hr>\n'
	html_str += '<h2>Mapping Stats</h2>\n'
	stats_df =  read_project_mapping_stats(ip_sample_list, con_sample_list, project_name)
	html_str += stats_df.to_html(border=1, index=False)
	
	#stats_list = [
	#	'Number of input reads',
	#	'Uniquely mapped reads number',
	#	'Uniquely mapped reads %',
	#	'Number of splices: Total',
	#	]
	#*--- EVALUATION FOR EACH IP-CON COMPARISON ---*	
	for x in comparison_list:
		ip_sample = x[0]
		con_sample = x[1]
		html_str += '<hr>\n'
		comparison = "{ip_sample}-{con_sample}".format(ip_sample=ip_sample, con_sample=con_sample)
		
		homer_table_unique = read_homer_table( os.path.join(pardir, 'projects', project_name,  'homer', comparison, 'clam_unique', 'homerResults.html') )
		topology_img_unique = os.path.join(pardir, 'projects', project_name, 'topology', comparison, 'clam_unique', 'dist.png' )
		repeat_img_unique = os.path.join(pardir, 'projects', project_name, 'repeats', comparison, 'clam_unique', 'dist.png' )

		homer_table_rescue = read_homer_table( os.path.join(pardir, 'projects', project_name, 'homer', comparison, 'clam_rescue', 'homerResults.html') )
		topology_img_rescue = os.path.join(pardir, 'projects', project_name, 'topology', comparison, 'clam_rescue', 'dist.png' )
		repeat_img_rescue = os.path.join(pardir, 'projects', project_name, 'repeats', comparison, 'clam_rescue', 'dist.png' )

		peak_num_img = os.path.join(pardir, 'projects', project_name, 'clam', 'peaks-'+comparison, 'peak_num.png' )
		
		html_str += '<h2>%s</h2>\n'%comparison
		# summary stats
		html_str += '<h3>Summary statistics</h3>\n'
		html_str += stats_df.loc[[ip_sample, con_sample]].to_html(border=1, index=False)
		html_str += '<img src="%s" />\n'%peak_num_img
		html_str += '<br>\n'
		# homer motifs
		html_str += '<h3>HOMER motifs</h3>\n'
		html_str += '<h4>Unique peaks</h4>\n'
		html_str += homer_table_unique
		if include_mread_analysis:
			html_str += '<h4>Rescue peaks</h4>\n'
			html_str += homer_table_rescue
		html_str += '<br>\n'
		# topology distribution
		html_str += '<h3>Topology distribution</h3>\n'
		html_str += '<h4>Unique peaks</h4>\n'
		html_str += '<img src="%s" />\n'%topology_img_unique
		if include_mread_analysis:
			html_str += '<h4>Rescue peaks</h4>\n'
			html_str += '<img src="%s" />\n'%topology_img_rescue
		# repetitive elements
		if include_mread_analysis:
			html_str += '<h3>Repetitive elements</h3>\n'
			html_str += '<h4>Unique peaks</h4>\n'
			html_str += '<img src="%s" />\n'%repeat_img_unique	
			html_str += '<h4>Rescue peaks</h4>\n'
			html_str += '<img src="%s" />\n'%repeat_img_rescue

	
	html_str += '</body>\n'
	html_str += '</html>\n'
	with open(outfn, 'w') as f:
		f.write(html_str)
	return

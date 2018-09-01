'''
read in peaks and compute t-statistic
ZZJ
8.27.2018
'''

import os
import pandas as pd
import numpy as np
import scipy.stats as ss
import pysam
from collections import defaultdict


def read_t2g(fn='/u/nobackup/yxing/NOBACKUP/frankwoe/CLAM_m6A_Snakemake/t2g_mm10.txt'):
	tmp = pd.read_table(fn, header=0)
	t2g = {}
	for i in range(tmp.shape[0]):
		t2g[tmp.iloc[i,0]] = tmp.iloc[i,1]
	return t2g

def read_tx_file(fn, t2g):
	tx = pd.read_table(fn, header=0)
	gene_dict = defaultdict(float)
	unmapped = 0
	for i in range(tx.shape[0]):
		t = tx.loc[i, 'target_id']
		t = t.split('.')[0]
		if t in t2g:
			g = t2g[t]
		else:
			unmapped += 1
			continue
		gene_dict[g] += tx.loc[i, 'tpm']
	if unmapped:
		print('%i tx nont mappable in this t2g index.'%unmapped)
	return gene_dict


def read_project_kallisto(_project_dir, t2g):
	project_dir = os.path.join(_project_dir, 'kallisto')
	tx_fn_list = [os.path.join(project_dir, x, 'abundance.tsv') for x in os.listdir(project_dir)]
	sample_gene_dict = {}
	for tx_fn in tx_fn_list:
		gene_dict = read_tx_file(tx_fn, t2g)
		sample_name = tx_fn.split('/')[-2]
		sample_gene_dict[sample_name] = gene_dict
	return sample_gene_dict


def read_peak_file(fn, peak_to_gene=None):
	peak_dict = {}
	if peak_to_gene is None:
		peak_to_gene = {}
	with open(fn, 'r') as f:
		for line in f:
			ele = line.strip().split()
			peak = ':'.join(ele[0:3])
			score = float(ele[-4])
			gene = ele[3].split('-')[0]
			peak_dict[peak] = score
			peak_to_gene[peak] =  gene
	return peak_dict, peak_to_gene


def read_project_clam(_project_dir):
	project_dir = os.path.join(_project_dir, 'clam')
	project_dict = {}
	peak_fn_list = [os.path.join(x, 'narrow_peak.unique.bed.all') for x in os.listdir(project_dir) if x.startswith('peaks-')]
	for peak_fn in peak_fn_list:
		peak_name = peak_fn.split('/')[-2].split('-')[1].rstrip('_IP')
		project_dict[peak_name] = read_peak_file(peak_fn)[0]

	df = pd.DataFrame.from_dict(project_dict)

	sig_peaks = set()
	peak_to_gene = dict()
	peak_fn_list2 = [os.path.join(x, 'narrow_peak.unique.bed') for x in os.listdir(project_dir) if x.startswith('peaks-')]
	for peak_fn in peak_fn_list2:
		tmp, peak_to_gene = read_peak_file(peak_fn, peak_to_gene)
		_ = map(sig_peaks.add, tmp.keys())

	df_with_sig = df.loc[sig_peaks]
	return df_with_sig, peak_to_gene


def wilcoxon_test(df_with_sig):
	group_comparisons = {
		'KI_gender': [
			['F_KI_06', 'F_KI_07', 'F_KI_08', 'F_KI_09', 'F_KI_10',],
			['M_KI_01', 'M_KI_02', 'M_KI_03', 'M_KI_04', 'M_KI_5',]
		],
		'V_position': [
			['LV_01', 'LV_02', 'LV_03', 'LV_04', 'LV_05', 'LV_06', 'LV_07', 'LV_08','LV_09','LV_10',],
			['RV_01', 'RV_02', 'RV_03', 'RV_04', 'RV_05', 'RV_06', 'RV_07', 'RV_08','RV_09','RV_10',],
		],
	}

	res = pd.DataFrame(1, index=df_with_sig.index, columns=list(group_comparisons.keys()))
	for peak in df_with_sig.index:
		for comparison in group_comparisons:
			a = df_with_sig.loc[peak, group_comparisons[comparison][0]]
			b = df_with_sig.loc[peak, group_comparisons[comparison][1]]
			a = a[~np.isnan(a)]
			b = b[~np.isnan(b)]
			fc = np.exp(abs(np.mean(a) - np.mean(b)))
			if fc>2:
				pv = ss.mannwhitneyu(a,b).pvalue
				res.loc[peak, comparison] = pv


def get_total_reads(bam_filename):
	idxstats  = pysam.idxstats(bam_filename).split('\n')
	tot = 0
	for l in idxstats:
		if not l:
			continue
		ele = l.split('\t')
		tot += int(ele[-2])
	return tot


def get_reads_window(bam, chrom, start, end):
	reads  = [ 1 for x in bam.fetch(chrom ,start, end) ]
	return len(reads)


def make_file_handlers(project_dir):
	bam_dict = {}
	bam_fn_list = [os.path.join(x, 'unique.sorted.bam') for x in os.listdir(project_dir) \
		if not x.startswith('peaks-') and os.path.isdir(x)]
	for bam_fn in bam_fn_list:
		sample_name = bam_fn.split('/')[-2]
		bam_dict[sample_name] = pysam.Samfile(bam_fn, 'rb')
	return bam_dict


def count(df_with_sig, bam_dict):
	'''
	Given a list of peaks and samples, count the RPKM,
	then save to file.
	'''
	input_readcounts = pd.DataFrame(0, index=df_with_sig.index, columns=df_with_sig.columns)
	ip_readcounts = pd.DataFrame(0, index=df_with_sig.index, columns=df_with_sig.columns)
	
	input_total_counts = {
		x+'_Input' : get_total_reads(x+'_Input/unique.sorted.bam')/float(10**6)
		 for x in df_with_sig.columns
		 }
	ip_total_counts = {
		x+'_IP':get_total_reads(x+'_IP/unique.sorted.bam')/float(10**6)
	 	for x in df_with_sig.columns
	 	}

	for peak in df_with_sig.index:
		chrom, start, end = peak.split(':')
		start = int(start)
		end = int(end)
		for sam in df_with_sig.columns:
			input_name = sam+'_Input'
			this_bam = bam_dict[input_name]
			this_input = get_reads_window(this_bam, chrom, start, end)
			input_readcounts.loc[peak, sam] = \
				this_input / float(input_total_counts[input_name]) / ((end-start)/1000.)

			ip_name = sam+'_IP'
			this_bam = bam_dict[ip_name]
			this_ip = get_reads_window(this_bam, chrom, start, end)
			ip_readcounts.loc[peak, sam] = \
				this_ip / float(ip_total_counts[ip_name]) / ((end-start)/1000.)
	return input_readcounts, ip_readcounts


def run_and_save():
	PROJECT_DIR = '.'
	T2G_FILE = '/u/nobackup/yxing/NOBACKUP/frankwoe/CLAM_m6A_Snakemake/t2g_mm10.txt'

	# read in kallisto gene level estimates
	t2g = read_t2g(T2G_FILE)
	gene_dict = read_project_kallisto(PROJECT_DIR, t2g)

	# read in clam peaks
	df_with_sig, peak_to_gene = read_project_clam(PROJECT_DIR)


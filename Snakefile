#!/usr/bin/env python

"""A pipeline for processing m6A data using CLAM
Author: 
	Zijun Zhang <zj.z@ucla.edu>
Date: 
	9.18.2017
Last revised:
	10.1.2017
"""

import os
import re

def get_sample_names(project):
	"""Get the sample names for IP and control in local folder
	"""
	IP_SAMPLE_LIST = []
	FQ_CMD = ''
	for fn in os.listdir('reads/{project}/ip/'.format(project=project)):
		IP_SAMPLE_LIST.append(re.search(r"(.+)\.fastq|\.gz", fn).group(1))
		if fn.endswith('gz'):
			FQ_CMD = '--readFilesCommand zcat'
	CON_SAMPLE_LIST = []
	for fn in os.listdir('reads/{project}/con/'.format(project=project)):
		CON_SAMPLE_LIST.append(re.search(r"(.+)\.fastq|\.gz", fn).group(1))
	if len(IP_SAMPLE_LIST)==0 or len(CON_SAMPLE_LIST)==0:
		raise Exception('either ip or ctrl sample has no raw reads found.')
	return IP_SAMPLE_LIST, CON_SAMPLE_LIST, FQ_CMD


GENOME = 'mm10'
PROJECT_NAME = os.environ.get("PROJECT_NAME")
IP_SAMPLE_LIST, CON_SAMPLE_LIST, FQ_CMD = get_sample_names(PROJECT_NAME)
IP_FQ_PATTERN = "reads/{project}/ip/{ip_sample}.fastq" if FQ_CMD=='' else "reads/{project}/ip/{ip_sample}.fastq.gz"
CON_FQ_PATTERN = "reads/{project}/con/{con_sample}.fastq" if FQ_CMD=='' else "reads/{project}/con/{con_sample}.fastq.gz"


###**--- SNAKEMAKE FILE ---**###

configfile: "config.yaml"

rule all:
	input:
		"archive/{project}.tar.gz".format(project=PROJECT_NAME)
				

rule alignment:
	input:
		expand( "star/{project}/ip/{ip_sample}/Aligned.sortedByCoord.out.bam", project=PROJECT_NAME, ip_sample=IP_SAMPLE_LIST ),
		expand( "star/{project}/con/{con_sample}/Aligned.sortedByCoord.out.bam", project=PROJECT_NAME, con_sample=CON_SAMPLE_LIST )

				
rule star_ip:
	input:
		sample=[IP_FQ_PATTERN]
	output:
		"star/{project}/ip/{ip_sample}/Aligned.sortedByCoord.out.bam"
	log:
		"logs/star/{ip_sample}.log"
	params:
		prefix="star/{project}/ip/{ip_sample}/",
		# path to STAR reference genome index
		index=config[GENOME]['star_idx'],
		fq_cmd = FQ_CMD
	threads: 4
	shell:
		"""
		STAR --genomeDir {params.index} \
		--readFilesIn {input.sample[0]}  --outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix {params.prefix} \
		--runThreadN {threads} \
		{params.fq_cmd}
		"""

rule star_con:
	input:
		sample=[CON_FQ_PATTERN]
	output:
		"star/{project}/con/{con_sample}/Aligned.sortedByCoord.out.bam"
	log:
		"logs/star/{con_sample}.log"
	params:
		prefix="star/{project}/con/{con_sample}/",
		index=config[GENOME]['star_idx'],
		fq_cmd = FQ_CMD
	threads: 4
	shell:
		"""
		STAR --genomeDir {params.index} \
		--readFilesIn {input.sample[0]}  --outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix {params.prefix} \
		--runThreadN {threads} \
		{params.fq_cmd}
		"""


rule mapping_stat_ip:
	input:
		align="star/{project}/ip/{sample}/Aligned.sortedByCoord.out.bam",
	output:
		"star/{project}/ip/{sample}/mapping_stats.txt"
	shell:
		"python2 scripts/report/mapping_stat.py {input.align} > {output}"


rule mapping_stat_con:
	input:
		align="star/{project}/con/{sample}/Aligned.sortedByCoord.out.bam",
	output:
		"star/{project}/con/{sample}/mapping_stats.txt"
	shell:
		"python2 scripts/report/mapping_stat.py {input.align} > {output}"


rule clam_prep_ip:
	input:
		ip_align="star/{project}/ip/{ip_sample}/Aligned.sortedByCoord.out.bam"
	output:
		"clam/{project}/ip/{ip_sample}/unique.sorted.bam"
	log:
		"logs/clam/{ip_sample}.log"
	params:
		outdir="clam/{project}/ip/{ip_sample}",
		tagger_method="median"
	shell:
		"CLAM preprocessor -i {input.ip_align} -o {params.outdir} --read-tagger-method {params.tagger_method}"
		

rule clam_prep_con:		
	input:
		con_align="star/{project}/con/{con_sample}/Aligned.sortedByCoord.out.bam"
	output:
		"clam/{project}/con/{con_sample}/unique.sorted.bam"
	log:
		"logs/clam/{con_sample}.log"
	params:
		outdir="clam/{project}/con/{con_sample}",
		tagger_method="median"
	shell:
		"CLAM preprocessor -i {input.con_align} -o {params.outdir} --read-tagger-method {params.tagger_method}"



rule clam_callpeak:
	input:
		ip_prep="clam/{project}/ip/{ip_sample}/unique.sorted.bam",
		con_prep="clam/{project}/con/{con_sample}/unique.sorted.bam"
	output:
		"clam/{project}/peaks-{ip_sample}-{con_sample}/narrow_peak.unique.bed"
	params:
		outdir="clam/{project}/peaks-{ip_sample}-{con_sample}",
		gtf=config[GENOME]['gtf'],
		binsize=100,
		qval_cutoff=0.005
	shell:
		"""
		CLAM peakcaller -i {input.ip_prep}  -c {input.con_prep} \
			-o {params.outdir} --gtf {params.gtf} --unique-only --unstranded --binsize {params.binsize} \
			--qval-cutoff {params.qval_cutoff}
		"""
		#mv log.CLAM.txt logs/clam/log.{wildcards.ip_sample}.{wildcards.con_sample}.txt


rule homer_motif:
	input:
		peak_fn="clam/{project}/peaks-{ip_sample}-{con_sample}/narrow_peak.unique.bed"
	output:
		"homer/{project}/{ip_sample}-{con_sample}/homerResults.html"
	params:
		outdir="homer/{project}/{ip_sample}-{con_sample}",
		motif_len='8',
		genome=GENOME,
		nthread=4,
		size=100,
		motif_num=5
	log:
		"logs/homer/log.homer.{ip_sample}-{con_sample}.txt"
	shell:
		"findMotifsGenome.pl {input.peak_fn} {params.genome} {params.outdir} "\
		" -rna -len {params.motif_len} "\
		"-p {params.nthread} -size {params.size} -S {params.motif_num}"


rule topology_dist:
	input:
		peak_fn="clam/{project}/peaks-{ip_sample}-{con_sample}/narrow_peak.unique.bed"
	output:
		dist_img="topology/{project}/{ip_sample}-{con_sample}/dist.png"
	log:
		"logs/topology/{ip_sample}-{con_sample}.log"
	params:
		outdir="topology/{project}/{ip_sample}-{con_sample}",
		count_script="scripts/topological_dist/Peak_distribution_on_utr_cds.py",
		plot_script="scripts/topological_dist/plot.R",
		genome=GENOME
	shell:
		"""
		python2 {params.count_script} {input} {params.genome} > {params.outdir}/dist.data
		Rscript {params.plot_script} {params.outdir}/dist.data {output.dist_img}
		"""

rule report:
	input:
		expand("star/{project}/ip/{sample}/mapping_stats.txt", project=PROJECT_NAME, sample=IP_SAMPLE_LIST),
		expand("star/{project}/con/{sample}/mapping_stats.txt", project=PROJECT_NAME, sample=CON_SAMPLE_LIST),
		expand( "homer/{project}/{ip_sample}-{con_sample}/homerResults.html", project=PROJECT_NAME, ip_sample=IP_SAMPLE_LIST, con_sample=CON_SAMPLE_LIST ),
		expand( "topology/{project}/{ip_sample}-{con_sample}/dist.png", project=PROJECT_NAME, ip_sample=IP_SAMPLE_LIST, con_sample=CON_SAMPLE_LIST )
	output:
		"reports/report_{project}.pdf".format(project=PROJECT_NAME)
	params:
		out_html="reports/report_{project}.html".format(project=PROJECT_NAME)
	run:
		from scripts.report import generate_report
		import pdfkit
		pardir = os.getcwd()
		generate_report.generate_report(IP_SAMPLE_LIST, CON_SAMPLE_LIST, pardir, params.out_html, PROJECT_NAME)
		pdfkit.from_file(params.out_html, output[0])

rule archive:
	input:
		peak_fn = expand("clam/{project}/peaks-{ip_sample}-{con_sample}/narrow_peak.unique.bed", 
			project=PROJECT_NAME, ip_sample=IP_SAMPLE_LIST, con_sample=CON_SAMPLE_LIST),
		report="reports/report_{project}.pdf".format(project=PROJECT_NAME),
	output:
		"archive/{project}.tar.gz".format(project=PROJECT_NAME)
	params:
		project=PROJECT_NAME
	shell:
		"""
		tar -czvf {output} clam/{params.project}/peaks-* {input.report}
		"""
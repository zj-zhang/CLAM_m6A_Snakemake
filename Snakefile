#!/usr/bin/env python

"""A pipeline for processing m6A data using CLAM
Author: 
	Zijun Zhang <zj.z@ucla.edu>
Date: 
	9.18.2017
Last revised:
	1.3.2018
Revision history:
	10.1.2017 null
	10.18.2017 add em-related functionalities 
	12.13.2017 changed starting file requirement to avoid re-run and enable starting from read fastq files
	1.3.2018 modified report
	2.14.2018 modified peak-calling filtering
	3.4.2018 added homer motif and topology; not in report yet
"""

import os
import re

###**--- UTILITY FUNCTION ---**###
def get_sample_names(project):
	"""Get the sample names for IP and control in local folder
	"""
	IP_SAMPLE_LIST = []
	FQ_CMD = ''
	for fn in os.listdir('projects/{project}/reads/ip/'.format(project=project)):
		IP_SAMPLE_LIST.append(re.search(r"(.+)\.fastq|\.gz", fn).group(1))
		if fn.endswith('gz'):
			FQ_CMD = '--readFilesCommand zcat'
	CON_SAMPLE_LIST = []
	for fn in os.listdir('projects/{project}/reads/con/'.format(project=project)):
		CON_SAMPLE_LIST.append(re.search(r"(.+)\.fastq|\.gz", fn).group(1))
	if len(IP_SAMPLE_LIST)==0 or len(CON_SAMPLE_LIST)==0:
		raise Exception('either ip or ctrl sample has no raw reads found.')
	return IP_SAMPLE_LIST, CON_SAMPLE_LIST, FQ_CMD


###**--- IN-FILE CONFIG ---**###
	
PROJECT = os.environ.get("PROJECT")
CONFIG_FP = os.path.join("projects", PROJECT, 'config', 'config.yaml')

configfile: 
	CONFIG_FP


PAIRED_END = False if 'paired_end' not in config else config['paired_end']
INCLUDE_MREAD_ANALYSIS = True if 'include_mread_analysis' not in config else config['include_mread_analysis']

FQ_PATTERN = [ 'projects/{project}/reads/{sample_type}/{sample_name}_1.fastq.gz', 'projects/{project}/reads/{sample_type}/{sample_name}_2.fastq.gz'] \
	if PAIRED_END else \
	['projects/{project}/reads/{sample_type}/{sample_name}.fastq.gz']	

GENOME = config['genome']
FQ_CMD = '--readFilesCommand zcat'

if PAIRED_END:
	IP_FQ_PATTERN = ["projects/{project}/reads/ip/{ip_sample}_1.fastq", "projects/{project}/reads/ip/{ip_sample}_2.fastq"] \
	if FQ_CMD=='' else \
	["projects/{project}/reads/ip/{ip_sample}_1.fastq.gz", "projects/{project}/reads/ip/{ip_sample}_2.fastq.gz"]
	CON_FQ_PATTERN = ["projects/{project}/reads/con/{con_sample}_1.fastq", "projects/{project}/reads/con/{con_sample}_2.fastq"] \ 
	if FQ_CMD=='' else \
	["projects/{project}/reads/con/{con_sample}_1.fastq.gz", "projects/{project}/reads/con/{con_sample}_2.fastq.gz"]
else:
	IP_FQ_PATTERN = "projects/{project}/reads/ip/{ip_sample}.fastq" if FQ_CMD=='' else "projects/{project}/reads/ip/{ip_sample}.fastq.gz"
	CON_FQ_PATTERN = "projects/{project}/reads/con/{con_sample}.fastq" if FQ_CMD=='' else "projects/{project}/reads/con/{con_sample}.fastq.gz"

SAMPLE_TYPE_DICT = config['sample_type_dict']
COMPARISON_LIST = config['sample_comparison']
MAX_TAGS = config['clam']['max_tags']


###**--- SNAKEMAKE FILE ---**###

rule all:
	input:
		"projects/{project}/archive/{project}.tar.gz".format(project=PROJECT)

rule all_clam:
	input:
		clam_peak = [
			"projects/{project}/clam/peaks-{ip_sample}-{con_sample}/narrow_peak.unique.bed".format( 
				project=PROJECT, ip_sample=x[0], con_sample=x[1])
				for x in COMPARISON_LIST
				],
		clam_mpeak = [
			"projects/{project}/clam/peaks-{ip_sample}-{con_sample}/narrow_peak.combined.bed".format( 
				project=PROJECT, ip_sample=x[0], con_sample=x[1])
				for x in COMPARISON_LIST
				]

### downloading ###

rule download_sra:
	input:
		#config_fn="projects/{project}/config/config.yaml".format(project=PROJECT),
		sradir_fn="projects/{project}/config/sra_dir.txt".format(project=PROJECT)
	output:
		["projects/{project}/sra/{sample_type}/{sample_name}.sra".format(project=PROJECT, sample_type=SAMPLE_TYPE_DICT[x], sample_name=x)
			for x in SAMPLE_TYPE_DICT.keys()]
	params:
		project=PROJECT,
		downloader_fn="projects/{project}/config/sra_downloader.sh".format(project=PROJECT),
		sradir_fn="projects/{project}/config/sra_dir.txt".format(project=PROJECT),
		outdir="projects/{project}/sra/".format(project=PROJECT)
	shell:
		"""
		#Rscript scripts/geo_downloader/getSRApath.R {params.project} {params.sradir_fn}
		cat {params.sradir_fn} | python2 scripts/geo_downloader/download_SRA.py {params.project} > {params.downloader_fn}
		bash {params.downloader_fn}
		#rm {params.sradir_fn}
		"""


rule dump_sra:
	input:
		lambda wildcards:
			"projects/{project}/sra/{sample_type}/{sample_name}.sra".format(
				project=PROJECT, 
				sample_type=SAMPLE_TYPE_DICT[wildcards.sample_name],
				sample_name=wildcards.sample_name)
	output:
		FQ_PATTERN
	params:
		outdir="projects/{project}/reads/{sample_type}/"
	shell:
		"fastq-dump -O {params.outdir} --gzip --split-3 {input}"
		

### alignment and preprocessing

# NOTE 2017.10.25: the current STAR rule runs on single-end even if the 
# data is paired-end. This is because in CLAM realigner, we cannot handle
# paired-end read yet and two reads with identical ID will cause issues.
# Will be fixed in the next update.		
rule star_ip:
	input:
		sample=[IP_FQ_PATTERN]
	output:
		"projects/{project}/star/ip/{ip_sample}/Aligned.sortedByCoord.out.bam"
	log:
		"projects/{project}/logs/star/{ip_sample}.log"
	params:
		prefix="projects/{project}/star/ip/{ip_sample}/",
		index=config[GENOME]['star_idx'],
		fq_cmd = FQ_CMD,
		max_hits = 100
	threads: 4
	shell:
		"""
		STAR --genomeDir {params.index} \
		--readFilesIn {input.sample[0]}  --outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix {params.prefix} \
		--outFilterMultimapNmax {params.max_hits} \
		--runThreadN {threads} \
		{params.fq_cmd} \
		--outStd Log >{log}		
		"""

rule star_con:
	input:
		sample=[CON_FQ_PATTERN]
	output:
		"projects/{project}/star/con/{con_sample}/Aligned.sortedByCoord.out.bam"
	log:
		"projects/{project}/logs/star/{con_sample}.log"
	params:
		prefix="projects/{project}/star/con/{con_sample}/",
		index=config[GENOME]['star_idx'],
		fq_cmd = FQ_CMD,
		max_hits = 100
	threads: 4
	shell:
		"""
		STAR --genomeDir {params.index} \
		--readFilesIn {input.sample[0]}  --outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix {params.prefix} \
		--outFilterMultimapNmax {params.max_hits} \
		--runThreadN {threads} \
		{params.fq_cmd} \
		--outStd Log >{log}
		"""


rule mapping_stat_ip:
	input:
		align="projects/{project}/star/ip/{sample}/Aligned.sortedByCoord.out.bam",
	output:
		"projects/{project}/star/ip/{sample}/mapping_stats.txt"
	shell:
		"python2 scripts/report/mapping_stat.py {input.align} > {output}"


rule mapping_stat_con:
	input:
		align="projects/{project}/star/con/{sample}/Aligned.sortedByCoord.out.bam",
	output:
		"projects/{project}/star/con/{sample}/mapping_stats.txt"
	shell:
		"python2 scripts/report/mapping_stat.py {input.align} > {output}"


rule clam_prep_ip:
	input:
		ip_align="projects/{project}/star/ip/{ip_sample}/Aligned.sortedByCoord.out.bam"
	output:
		"projects/{project}/clam/ip/{ip_sample}/unique.sorted.collapsed.bam" if MAX_TAGS>0 else \
			"projects/{project}/clam/ip/{ip_sample}/unique.sorted.bam"
	log:
		"projects/{project}/logs/clam/{ip_sample}-prep.log"
	params:
		outdir="projects/{project}/clam/ip/{ip_sample}",
		tagger_method="median",
		max_tags=MAX_TAGS
	shell:
		"CLAM preprocessor -i {input.ip_align} -o {params.outdir} --max-tags {params.max_tags} --read-tagger-method {params.tagger_method} >{log} 2>&1"
		

rule clam_prep_con:		
	input:
		con_align="projects/{project}/star/con/{con_sample}/Aligned.sortedByCoord.out.bam"
	output:
		"projects/{project}/clam/con/{con_sample}/unique.sorted.collapsed.bam" if MAX_TAGS>0 else \
			"projects/{project}/clam/con/{con_sample}/unique.sorted.bam"
	log:
		"projects/{project}/logs/clam/{con_sample}-prep.log"
	params:
		outdir="projects/{project}/clam/con/{con_sample}",
		tagger_method="median",
		max_tags=MAX_TAGS
	shell:
		"CLAM preprocessor -i {input.con_align} -o {params.outdir} --max-tags {params.max_tags} --read-tagger-method {params.tagger_method} >{log} 2>&1"

rule clam_em_ip:
	input:
		"projects/{project}/clam/ip/{ip_sample}/unique.sorted.collapsed.bam" if MAX_TAGS>0 else \
			"projects/{project}/clam/ip/{ip_sample}/unique.sorted.bam"
	output:
		"projects/{project}/clam/ip/{ip_sample}/realigned.sorted.bam"
	log:
		"projects/{project}/logs/clam/{ip_sample}-em.log"
	params:
		outdir="projects/{project}/clam/ip/{ip_sample}",
		max_tags=MAX_TAGS,
		winsize=100
	shell:
		"CLAM realigner -i {input} -o {params.outdir} --winsize {params.winsize} --max-tags {params.max_tags} --unstranded >{log} 2>&1"

rule clam_em_con:
	input:
		"projects/{project}/clam/con/{con_sample}/unique.sorted.collapsed.bam" if MAX_TAGS>0 else \
			"projects/{project}/clam/con/{con_sample}/unique.sorted.bam"
	output:
		"projects/{project}/clam/con/{con_sample}/realigned.sorted.bam"
	log:
		"projects/{project}/logs/clam/{con_sample}-em.log"
	params:
		outdir="projects/{project}/clam/con/{con_sample}",
		max_tags=MAX_TAGS,
		winsize=100
	shell:
		"CLAM realigner -i {input} -o {params.outdir} --winsize {params.winsize} --max-tags {params.max_tags} --unstranded >{log} 2>&1"
			

### peak calling

rule clam_callpeak:
	input:
		ip_prep="projects/{project}/clam/ip/{ip_sample}/unique.sorted.collapsed.bam" if MAX_TAGS>0 else \
			"projects/{project}/clam/ip/{ip_sample}/unique.sorted.bam",
		con_prep="projects/{project}/clam/con/{con_sample}/unique.sorted.collapsed.bam" if MAX_TAGS>0 else \
			"projects/{project}/clam/con/{con_sample}/unique.sorted.bam"
	output:
		"projects/{project}/clam/peaks-{ip_sample}-{con_sample}/narrow_peak.unique.bed"
	log:
		"projects/{project}/logs/clam/{ip_sample}-{con_sample}-callpeak.log"
	params:
		outdir="projects/{project}/clam/peaks-{ip_sample}-{con_sample}",
		gtf=config[GENOME]['gtf'],
		binsize=100,
		qval_cutoff=0.1,
		fold_change='0.1',
		threads=4
		
	shell:
		"""
		CLAM peakcaller -i {input.ip_prep}  -c {input.con_prep} \
			-p {params.threads} \
			-o {params.outdir} --gtf {params.gtf} --unique-only --unstranded --binsize {params.binsize} \
			--qval-cutoff {params.qval_cutoff} --fold-change {params.fold_change} >{log} 2>&1
		mv {output} {output}.all
		awk '$9<0.005 && $7>1.609' {output}.all > {output}
		"""
		


rule clam_callpeak_mread:
	input:		
		ip_prep=[
			"projects/{project}/clam/ip/{ip_sample}/unique.sorted.collapsed.bam" if MAX_TAGS>0 else "projects/{project}/clam/ip/{ip_sample}/unique.sorted.bam", 
			"projects/{project}/clam/ip/{ip_sample}/realigned.sorted.bam"],
		con_prep=[
			"projects/{project}/clam/con/{con_sample}/unique.sorted.collapsed.bam" if MAX_TAGS>0 else "projects/{project}/clam/con/{con_sample}/unique.sorted.bam", 
			"projects/{project}/clam/con/{con_sample}/realigned.sorted.bam"]
	output:
		"projects/{project}/clam/peaks-{ip_sample}-{con_sample}/narrow_peak.combined.bed"
	log:
		"projects/{project}/logs/clam/{ip_sample}-{con_sample}-callpeak_mread.log"
	params:
		outdir="projects/{project}/clam/peaks-{ip_sample}-{con_sample}",
		gtf=config[GENOME]['gtf'],
		binsize=100,
		qval_cutoff=0.1,
		fold_change='0.1',
		threads=4
	shell:
		"""
		CLAM peakcaller -i {input.ip_prep}  -c {input.con_prep} \
			-p {params.threads} \
			-o {params.outdir} --gtf {params.gtf} --unstranded --binsize {params.binsize} \
			--qval-cutoff {params.qval_cutoff} --fold-change {params.fold_change} >{log} 2>&1
		mv {output} {output}.all
		awk '$9<0.005 && $7>1.609' {output}.all > {output}
		"""
		
		

rule macs2:
	input:
		ip="projects/{project}/star/ip/{ip_sample}/Aligned.sortedByCoord.out.bam",
		con="projects/{project}/star/con/{con_sample}/Aligned.sortedByCoord.out.bam"
	output:
		"projects/{project}/macs2/{ip_sample}-{con_sample}/{ip_sample}-{con_sample}--nomodel_peaks.narrowPeak"
	log:
		"projects/{project}/logs/macs2/{ip_sample}-{con_sample}.log"
	params:
		genome_size="2.7e9",
		name_="{ip_sample}-{con_sample}",
		extsize=50,
		fold_range="9 30",
		q_cutoff=0.01,
		outdir="projects/{project}/macs2/{ip_sample}-{con_sample}"
	shell:
		"macs2 callpeak -g {params.genome_size} -t {input.ip} -c {input.con} -n {params.name_}" \
		"--nomodel -m {params.fold_range} -q {params.q_cutoff} --extsize {params.extsize} --outdir {params.outdir} >{log} 2>&1"


rule compare_peaks:
	input:
		macs2_peak="projects/{project}/macs2/{ip_sample}-{con_sample}/{ip_sample}-{con_sample}--nomodel_peaks.narrowPeak",
		clam_upeak="projects/{project}/clam/peaks-{ip_sample}-{con_sample}/narrow_peak.unique.bed",
		clam_mpeak="projects/{project}/clam/peaks-{ip_sample}-{con_sample}/narrow_peak.combined.bed"
	output:
		clam_rescued_peak="projects/{project}/clam/peaks-{ip_sample}-{con_sample}/narrow_peak.rescue.bed",
		clam_macs2_peak="projects/{project}/clam/peaks-{ip_sample}-{con_sample}/narrow_peak.macs2.bed",
		#macs2_clam_peak="projects/{project}/macs2/{ip_sample}-{con_sample}/narrow_peak.clam.bed",
		plot_fn="projects/{project}/clam/peaks-{ip_sample}-{con_sample}/peak_num.png"
	params:
		plot_script="scripts/report/plot_peak_num.R",
	shell:
		# removed 1.3.2018: no output was generated; need debug
		# bedtools intersect -b {input.clam_mpeak} -a {input.macs2_peak} > {output.macs2_clam_peak} 
		"""
		bedtools intersect -v -a {input.clam_mpeak} -b {input.clam_upeak} > {output.clam_rescued_peak}
		bedtools intersect -a {input.clam_mpeak} -b {input.macs2_peak} > {output.clam_macs2_peak}
		Rscript {params.plot_script} {input.clam_mpeak} {input.clam_upeak} {input.macs2_peak} {output.plot_fn}
		"""

### evaluations

rule homer_motif:
	input:
		peak_fn="projects/{project}/clam/peaks-{ip_sample}-{con_sample}/narrow_peak.unique.bed"
	output:
		"projects/{project}/homer/{ip_sample}-{con_sample}/clam_unique/homerResults.html"
	params:
		outdir="projects/{project}/homer/{ip_sample}-{con_sample}/clam_unique",
		motif_len='5,6,7',
		genome=GENOME,
		nthread=4,
		size=100,
		motif_num=10
	log:
		"projects/{project}/logs/homer/log.homer.{ip_sample}-{con_sample}.txt"
	shell:
		"findMotifsGenome.pl {input.peak_fn} {params.genome} {params.outdir} "\
		" -rna -len {params.motif_len} "\
		"-p {params.nthread} -size {params.size} -S {params.motif_num} >{log} 2>&1"


rule topology_dist:
	input:
		peak_fn="projects/{project}/clam/peaks-{ip_sample}-{con_sample}/narrow_peak.unique.bed"
	output:
		dist_img="projects/{project}/topology/{ip_sample}-{con_sample}/clam_unique/dist.png"
	log:
		"projects/{project}/logs/topology/{ip_sample}-{con_sample}.log"
	params:
		outdir="projects/{project}/topology/{ip_sample}-{con_sample}/clam_unique",
		count_script="scripts/topological_dist/Peak_distribution_on_utr_cds.py",
		plot_script="scripts/topological_dist/plot.R",
		genome=GENOME,
		binnum=50
	shell:
		"""
		python2 {params.count_script} {input} {params.genome} {params.binnum} >{params.outdir}/dist.data 2>{log}
		Rscript {params.plot_script} {params.outdir}/dist.data {output.dist_img}
		"""


rule repeat_comp:
	input:
		peak_fn="projects/{project}/clam/peaks-{ip_sample}-{con_sample}/narrow_peak.unique.bed"
	output:
		"projects/{project}/repeats/{ip_sample}-{con_sample}/clam_unique/dist.png"
	params:
		genome=GENOME,
		outdir="projects/{project}/repeats/{ip_sample}-{con_sample}/clam_unique",
		count_script="scripts/repeat_composition/RepeatsPie.py",
		plot_script="scripts/repeat_composition/plot.R"
	shell:
		"""
		python2 {params.count_script} {input.peak_fn} {params.genome} >{params.outdir}/dist.data
		Rscript {params.plot_script} {params.outdir}/dist.data {output}
		"""


rule homer_motif_rescue:
	input:
		peak_fn="projects/{project}/clam/peaks-{ip_sample}-{con_sample}/narrow_peak.rescue.bed"
	output:
		"projects/{project}/homer/{ip_sample}-{con_sample}/clam_rescue/homerResults.html"
	params:
		outdir="projects/{project}/homer/{ip_sample}-{con_sample}/clam_rescue",
		motif_len='5,6,7',
		genome=GENOME,
		nthread=4,
		size=100,
		motif_num=10
	log:
		"projects/{project}/logs/homer/log.homer.{ip_sample}-{con_sample}.rescue.txt"
	shell:
		"findMotifsGenome.pl {input.peak_fn} {params.genome} {params.outdir} "\
		" -rna -len {params.motif_len} "\
		"-p {params.nthread} -size {params.size} -S {params.motif_num} >{log} 2>&1"


rule topology_dist_rescue:
	input:
		peak_fn="projects/{project}/clam/peaks-{ip_sample}-{con_sample}/narrow_peak.rescue.bed"
	output:
		dist_img="projects/{project}/topology/{ip_sample}-{con_sample}/clam_rescue/dist.png"
	log:
		"projects/{project}/logs/topology/{ip_sample}-{con_sample}-rescue.log"
	params:
		outdir="projects/{project}/topology/{ip_sample}-{con_sample}/clam_rescue",
		count_script="scripts/topological_dist/Peak_distribution_on_utr_cds.py",
		plot_script="scripts/topological_dist/plot.R",
		genome=GENOME,
		binnum=20
	shell:
		"""
		python2 {params.count_script} {input} {params.genome} {params.binnum} >{params.outdir}/dist.data 2>{log}
		Rscript {params.plot_script} {params.outdir}/dist.data {output.dist_img}
		"""		


rule repeat_comp_rescue:
	input:
		peak_fn="projects/{project}/clam/peaks-{ip_sample}-{con_sample}/narrow_peak.rescue.bed"
	output:
		"projects/{project}/repeats/{ip_sample}-{con_sample}/clam_rescue/dist.png"
	params:
		genome=GENOME,
		outdir="projects/{project}/repeats/{ip_sample}-{con_sample}/clam_rescue",
		count_script="scripts/repeat_composition/RepeatsPie.py",
		plot_script="scripts/repeat_composition/plot.R"
	shell:
		"""
		python2 {params.count_script} {input.peak_fn} {params.genome} >{params.outdir}/dist.data
		Rscript {params.plot_script} {params.outdir}/dist.data {output}
		"""


rule homer_motif_macs2:
	input:
		peak_fn="projects/{project}/macs2/{ip_sample}-{con_sample}/{ip_sample}-{con_sample}--nomodel_peaks.narrowPeak"
	output:
		"projects/{project}/homer/{ip_sample}-{con_sample}/macs2/homerResults.html"
	params:
		outdir="projects/{project}/homer/{ip_sample}-{con_sample}/macs2",
		motif_len='5,6,7',
		genome=GENOME,
		nthread=4,
		size=100,
		motif_num=10
	log:
		"projects/{project}/logs/homer/log.homer.{ip_sample}-{con_sample}.macs2.txt"
	shell:
		"findMotifsGenome.pl {input.peak_fn} {params.genome} {params.outdir} "\
		"-len {params.motif_len} "\
		"-p {params.nthread} -size {params.size} -S {params.motif_num} >{log} 2>&1"


rule topology_dist_macs2:
	input:
		peak_fn="projects/{project}/macs2/{ip_sample}-{con_sample}/{ip_sample}-{con_sample}--nomodel_peaks.narrowPeak"
	output:
		dist_img="projects/{project}/topology/{ip_sample}-{con_sample}/macs2/dist.png"
	log:
		"projects/{project}/logs/topology/{ip_sample}-{con_sample}-macs2.log"
	params:
		outdir="projects/{project}/topology/{ip_sample}-{con_sample}/macs2",
		count_script="scripts/topological_dist/Peak_distribution_on_utr_cds.py",
		plot_script="scripts/topological_dist/plot.R",
		genome=GENOME
	shell:
		"""
		python2 {params.count_script} {input} {params.genome} > {params.outdir}/dist.data
		Rscript {params.plot_script} {params.outdir}/dist.data {output.dist_img}
		"""		

rule make_bw:
	input:
		#"projects/{project}/clam/peaks-{ip_sample}-{con_sample}/narrow_peak.combined.bed"
		"projects/{project}/clam/con/{con_sample}/realigned.sorted.bam",	
		"projects/{project}/clam/ip/{ip_sample}/realigned.sorted.bam"	
	output:
		"projects/{project}/bigwig/{ip_sample}-{con_sample}/foo.txt"
	params:
		ip_mbam="projects/{project}/clam/ip/{ip_sample}/realigned.sorted.bam",
		ip_ubam="projects/{project}/clam/ip/{ip_sample}/unique.sorted.collapsed.bam" if MAX_TAGS>0 else \
			"projects/{project}/clam/ip/{ip_sample}/unique.sorted.bam",
		ip_bw_dir="projects/{project}/bigwig/{ip_sample}-{con_sample}/{ip_sample}/",
		con_mbam="projects/{project}/clam/con/{con_sample}/realigned.sorted.bam",
		con_ubam="projects/{project}/clam/con/{con_sample}/unique.sorted.collapsed.bam" if MAX_TAGS>0 else \
			"projects/{project}/clam/con/{con_sample}/unique.sorted.bam",
		con_bw_dir="projects/{project}/bigwig/{ip_sample}-{con_sample}/{con_sample}/",
		bw_script="scripts/make_bw/make_bigwig.py"
	shell:
		"""
		mkdir -p {params.ip_bw_dir}
		mkdir -p {params.con_bw_dir}
		python2 {params.bw_script} {params.ip_ubam} {params.ip_mbam} {params.ip_bw_dir}
		python2 {params.bw_script} {params.con_ubam} {params.con_mbam} {params.con_bw_dir}
		echo "`date` done making bw" > {output}
		"""




### generating reports and cleaning up


rule report:
	input:
		# require ip mapping stats
		[ "projects/{project}/star/ip/{sample}/mapping_stats.txt".format(
			project=PROJECT, 
			sample=x[0])
			for x in COMPARISON_LIST
		],
		# require con mapping stats
		[ "projects/{project}/star/con/{sample}/mapping_stats.txt".format(
			project=PROJECT, 
			sample=x[1])
			for x in COMPARISON_LIST
		],
		# require peak comparison
		[ "projects/{project}/clam/peaks-{ip_sample}-{con_sample}/peak_num.png".format(
			project=PROJECT, 
			ip_sample=x[0],
			con_sample=x[1])
			for x in COMPARISON_LIST
		],
		# require unique homer
		[ "projects/{project}/homer/{ip_sample}-{con_sample}/clam_unique/homerResults.html".format(
			project=PROJECT, 
			ip_sample=x[0], 
			con_sample=x[1] )
			for x in COMPARISON_LIST
		],
		# require unique topology
		[ "projects/{project}/topology/{ip_sample}-{con_sample}/clam_unique/dist.png".format(
			project=PROJECT, 
			ip_sample=x[0], 
			con_sample=x[1] )
			for x in COMPARISON_LIST
		],
		# require unique repeats
		[ "projects/{project}/repeats/{ip_sample}-{con_sample}/clam_unique/dist.png".format(
			project=PROJECT, 
			ip_sample=x[0], 
			con_sample=x[1] )
			for x in COMPARISON_LIST
		],
		# require rescued homer
		[ "projects/{project}/homer/{ip_sample}-{con_sample}/clam_rescue/homerResults.html".format(
			project=PROJECT, 
			ip_sample=x[0], 
			con_sample=x[1] )
			for x in COMPARISON_LIST
		],
		# require rescued topology
		[ "projects/{project}/topology/{ip_sample}-{con_sample}/clam_rescue/dist.png".format(
			project=PROJECT, 
			ip_sample=x[0], 
			con_sample=x[1] )
			for x in COMPARISON_LIST
		],
		# require rescued repeats
		[ "projects/{project}/repeats/{ip_sample}-{con_sample}/clam_rescue/dist.png".format(
			project=PROJECT, 
			ip_sample=x[0], 
			con_sample=x[1] )
			for x in COMPARISON_LIST
		],
		# require macs2 homer motif
		[ "projects/{project}/homer/{ip_sample}-{con_sample}/macs2/homerResults.html".format(
			project=PROJECT, 
			ip_sample=x[0], 
			con_sample=x[1] )
			for x in COMPARISON_LIST
		],
		# require macs2 topology
		[ "projects/{project}/topology/{ip_sample}-{con_sample}/macs2/dist.png".format(
			project=PROJECT, 
			ip_sample=x[0], 
			con_sample=x[1] )
			for x in COMPARISON_LIST
		]
	output:
		"projects/{project}/reports/report_{project}.pdf".format(project=PROJECT)
	params:
		out_html="projects/{project}/reports/report_{project}.html".format(project=PROJECT)
	run:
		from scripts.report import generate_report
		import pdfkit
		pardir = os.getcwd()
		generate_report.generate_report(COMPARISON_LIST, pardir, params.out_html, PROJECT, include_mread_analysis=INCLUDE_MREAD_ANALYSIS)
		pdfkit.from_file(params.out_html, output[0])

rule archive:
	input:
		clam_peak = [
			"projects/{project}/clam/peaks-{ip_sample}-{con_sample}/narrow_peak.unique.bed".format( 
				project=PROJECT, ip_sample=x[0], con_sample=x[1])
				for x in COMPARISON_LIST
				],
		macs2_peak = [
				"projects/{project}/macs2/{ip_sample}-{con_sample}/{ip_sample}-{con_sample}--nomodel_peaks.narrowPeak".format(
				project=PROJECT, ip_sample=x[0], con_sample=x[1])
				for x in COMPARISON_LIST
				],
		bigwig = [ "projects/{project}/bigwig/{ip_sample}-{con_sample}/foo.txt".format(
				project=PROJECT,
				ip_sample=x[0],
				con_sample=x[1] )
				for x in COMPARISON_LIST
				],
		report = "projects/{project}/reports/report_{project}.pdf".format(project=PROJECT),
	output:
		"projects/{project}/archive/{project}.tar.gz".format(project=PROJECT)
	params:
		project=PROJECT
	shell:
		"""
		tar -czvf {output} projects/{params.project}/clam/peaks-* projects/{params.project}/macs2/* {input.report} projects/{params.project}/bigwig/*
		##rm -rf projects/{params.project}/sra/* projects/{params.project}/star/*/*/*.bam projects/{params.project}/reads/*
		echo "`date` done archiving" > projects/{params.project}/foo.txt
		"""
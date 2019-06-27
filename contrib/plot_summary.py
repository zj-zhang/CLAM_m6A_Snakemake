""" update PDF generator for m6A summary report
Author:ZCP
Date: May 14, 2019
"""
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
from peakComposition.peakPieCharts import intersect_gtf_regions

import re
import numpy as np
from subprocess import *
import time
from collections import defaultdict
from scipy.stats import zscore
import yaml
from PyPDF2 import PdfFileMerger, PdfFileReader

from itertools import combinations
from pybedtools import BedTool

try:
    PROJECT_DIR = sys.argv[1]
except:
    sys.exit('command insufficient project dir')

config = yaml.safe_load(open(PROJECT_DIR + '/config/config.yaml'))

cur_dir = os.path.dirname(os.path.realpath(__file__))

SNAKEMAKE_FILE_DIR = config['SCRIPTS_DIR']
GZIP=False
PAIRED_END = False if 'paired_end' not in config else config['paired_end']
INCLUDE_MREAD_ANALYSIS = True if 'include_mread_analysis' not in config else config['include_mread_analysis']

OUT_DIR = PROJECT_DIR + '/reports/'
TMP_DIR = PROJECT_DIR + '/tmp/'

GENOME = config['genome']

SAMPLE_TYPE_DICT = config['sample_type_dict']
COMPARISON_LIST = config['sample_comparison']
MAX_TAGS = config['clam']['max_tags']

def makedirs(_dir):
    try:
        os.stat(_dir)
    except:
        os.makedirs(_dir)

def read_file(filename):
    with open(filename) as fp:
        List = [x.strip() for x in fp if len(x.strip()) > 0]
        return List

def convert2log2pwm(pwm, scale = False):
    if scale:
        scale_size = 0.01
    else:
        scale_size = 0.00
    array = np.array([[float(_x) for _x in x.split('\t')] for x in pwm.strip().split('\n')[1:]])
    array += scale_size
    array /= np.sum(array, axis = 1).reshape((len(array),1))
    array /= 0.25
    return array

def star_mapping_stats():
    makedirs(OUT_DIR)
    fw = open(OUT_DIR + 'mapping_stats.counts.summary.txt', 'w')
    fw.write('sample\ttype\ttotal reads\tsplice junction reads\texon reads\tuniquely mapped reads (%)\treads mapped to multiple loci(%)\n')
    star_dir = PROJECT_DIR + '/star/'
    for sample in sorted(SAMPLE_TYPE_DICT.keys()):
        sample_type = SAMPLE_TYPE_DICT[sample]
        if sample_type == 'ip':
            sample_type = 'IP'
        if sample_type == 'con':
            sample_type = 'Inp'
        mapping_info = read_file(star_dir + sample + '/mapping_stats.txt')        
        input_reads = "{:,}".format(int(mapping_info[0].split()[-1]))
        splice_junction_reads = "{:,}".format(int(mapping_info[4].split()[-1]))
        exon_reads = "{:,}".format(int(mapping_info[5].split()[-1]))
        unique_mapped = mapping_info[1].split()[-1]
        multip_mapped = mapping_info[2].split()[-1]
        fw.write('\t'.join([sample, sample_type, input_reads, splice_junction_reads, exon_reads, unique_mapped, multip_mapped]) + '\n')
    fw.close()

def table_generate(table_inp, out_file):
    table_row = len(read_file(table_inp))
    print "library(gtable)"
    print "library(gridExtra)"
    print "library(grid)"
    print "data <- read.table('{}', header = T, sep = '\\t', check.names=FALSE)".format(table_inp)
    print "g <- tableGrob(data, rows = NULL)"
    print "g <- gtable_add_grob(g,"
    print "                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),"
    print "                     t = 2, b = nrow(g), l = 1, r = ncol(g))"
    print "g <- gtable_add_grob(g,"
    print "                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),"
    print "                     t = 1, l = 1, r = ncol(g))"
    print "pdf('{}', height = {}, width = 15)".format(out_file, round(table_row * 0.3))
    print "grid.draw(g)"
    print "dev.off()"

def peaks_num_stats():
    peaks_dir = PROJECT_DIR + '/clam/peaks/'
    fw = open(OUT_DIR + 'peaks.counts.summary.txt', 'w')
    fw.write('sample\tpeak_num\n')
    for samples in COMPARISON_LIST:
        peaks_num = len(read_file(peaks_dir + 'peaks-{}_{}/narrow_peak.unique.bed'.format(samples[0], samples[1])))
        fw.write('{}\t{}\n'.format('-'.join(samples[0].split('-')[:-1]), peaks_num))
    fw.close()

def load_theme():
    print 'library(ggplot2)'
    print 'library(ggthemes)'
    for line in read_file(cur_dir + '/theme_pub.R'):
        print line.strip()

def plot_peaks_num(out_file):
    load_theme()
    print "data <- read.table('{}', header = T , sep = '\\t')".format(OUT_DIR + 'peaks.counts.summary.txt')
    print "p1 <- ggplot(data, aes(x = sample, y = peak_num)) + geom_bar(stat='identity', aes(fill = peak_num)) + scale_fill_Publication() +theme_Publication() + "
    print "  theme(axis.text.x = element_text(angle = 45, hjust = 1),  legend.position = 'none',"
    print "        strip.text = element_text(size=15), axis.title.x = element_blank(), legend.title = element_blank()) + "
    print "  ylab('m6A Peak number') + scale_y_continuous(expand = c(0, 0)) + scale_fill_gradient2()"
    print "print(p1)"
    print "ggsave('{}', height = 5, width = 8)".format(out_file)
    print "dev.off()"

def get_motif_rScript(motif_files, title, index):
    print 'seq_matrix_list <- list()'
    for rank, motif_file in enumerate(motif_files):
        pvalue = read_file(motif_file)[0].split(':')[-1]
        pwm = '\n'.join(read_file(motif_file))
        pssm = np.array([[float(_x) for _x in x.split('\t')] for x in pwm.strip().split('\n')[1:]]).T
        pssm_flatten = pssm.flatten()
        seq_len = len(pssm[0])
        motif_strip_text = 'Rank ' + str(rank + 1) + '\np-value = {}'.format(pvalue)
        print 'seq_profile = c(' + ','.join([str(x) for x in pssm_flatten]) + ')'
        print 'seq_matrix = matrix(seq_profile, 4, {}, byrow = T)'.format(seq_len)
        print "rownames(seq_matrix) <- c('A', 'C', 'G', 'U')"
        print 'seq_matrix_list$`{}` <- seq_matrix'.format(motif_strip_text)
    print "library(ggseqlogo)"
    print "library(ggplot2)"
    print "p{} <- ggplot() + geom_logo(seq_matrix_list, method = 'prob') + ".format(index)
    if len(motif_files) > 1:
        print "   facet_wrap(~seq_group,  scales = 'free_x', ncol = {}) + ".format(len(motif_files))
    print "       theme_logo() + theme(axis.text.x = element_blank(), "
    print "       panel.spacing = unit(0.5, 'lines'),"
    print "       axis.text.y = element_blank(), "
    print "       axis.title.y = element_blank(),"
    print "       axis.title.x = element_blank(), "
    print "       plot.title = element_text(hjust = 0.5, size=30, face='bold'),"
    print "       legend.position = 'none') + ggtitle('{}')".format(title)

def get_dist_topology(outfn, index):
    print "data_fn='{}'".format(outfn)
    print "data = read.table(data_fn)"
    print "tot = nrow(data)"
    print "sep1 = tot/3+1"
    print "sep2 = tot/3*2+1"
    print "library(ggplot2)"
    print "colnames(data) = c('Peak_Tx', 'Peak_Tx_WithPeak', 'NumPeak')"
    print "data$coord = seq(1, nrow(data))"
    print "p{} <- ggplot(data, aes(x=coord, y=Peak_Tx_WithPeak, fill=NumPeak, colour=NumPeak)) + geom_point(colour='#386cb0', size=2) +".format(index)
    print "  geom_line() +"
    print "  guides(fill=F, colour=F) +"
    print "  theme_bw() +"
    print "  geom_vline(xintercept=c(sep1, sep2), linetype='dashed', size = 1.2) +"
    print "  annotate('text', label=c(\"5'UTR\", \"CDS\", \"3'UTR\"), size = 7, "
    print "           x=c( sep1/2, (sep1+sep2)/2, sep2+sep1/2 ), y=max(data[,2])*0.9) +  theme(axis.text.x = element_blank(), "
    print "                                                                                    panel.spacing = unit(0.5, 'lines'),"
    print "                                                                                    panel.grid.major = element_blank(), "
    print "                                                                                    panel.grid.minor = element_blank(),"
    print "                                                                                    axis.title.x = element_blank(), "
    print "                                                                                    text = element_text(size=15),"
    print "                                                                                    axis.ticks.x = element_blank(),"
    print "                                                                                    legend.position = 'none') +"
    print "  scale_x_continuous(expand=c(0,0)) + ylab('Number of m6A peaks per transcript')"

def get_intersection_pie(inpfn, index):
    inp_list = read_file(inpfn)[1:]
    labels = []
    names = []
    for _inp in inp_list:
        sp = _inp.strip().split('\t')
        labels.append(str(round(float(sp[-1]) * 100, 1)) +  '%')
        names.append(sp[0])
    labels[-1] = names[-1] + ' ' + labels[-1]
    label_info = ', '.join(["'{}'".format(label) for label in labels[::-1]])
    print "library(ggrepel)"
    print "data <- read.table('{}', header = T, sep = '\\t', quote=\"\\\"\")".format(inpfn)
    print "data$category <- factor(data$category, levels = c('3\\'UTR', '5\\'UTR', 'CDS', 'Other exon', 'Intron'))"
    print "data <- data[c(5, 4, 3, 2, 1),]"
    print "data$pos = (cumsum(c(0, data$norm_count)) + c(data$norm_count /2, .01))[1:nrow(data)]"
    print "data$x <- c(1.5, 1.2, 1, 1, 1)"
    print "data$label <- c({})".format(label_info)
    print "data$nudge_x <- c(0.3, 0, 0.0, 0.0, 0.0)"
    print "data$nudge_y <- c(-0.05, 0, 0, 0, 0)"
    print 'data$label_color <- c("black","white","white","white","white")'
    print 'data$color <- rev(c("#7fc97f","#e79c47","#ef3b2c","#386cb0","#662506"))'
    print "p{} <- ggplot(data, aes(\"\", norm_count, fill = category)) + ".format(index)
    print "  geom_bar(width = 1, size = 1, color = 'white', stat = 'identity') +"
    print '  coord_polar("y") + '
    print "  geom_text_repel(aes(x = x, y = pos, label = label), color = data$label_color,"
    print "                  nudge_x = data$nudge_x,"
    print "                  nudge_y = data$nudge_y,"
    print "                  segment.size = .7, "
    print "                  fontface = 'bold',"
    print "                  size = 8,"
    print "                  show.legend = FALSE) +"
    print "  scale_fill_manual(values = data$color) +"
    print "  theme_classic() +"
    print "  theme(axis.line = element_blank(),"
    print "        axis.text = element_blank(),"
    print "        axis.ticks = element_blank(),"
    print "        legend.key = element_rect(colour = NA),"
    print "        legend.position = c(0.5, 0.03),"
    print "        legend.direction = 'horizontal',"
    print "        legend.key.size= unit(0.2, 'cm'),"
    print "        legend.title = element_blank(),"
    print "        legend.spacing.x = unit(0.2, 'cm'),"
    print "        legend.key.width =unit(1,'line'),"
    print "        legend.text=element_text(size=20),"
    print "        axis.title.x = element_blank(),"
    print "        axis.title.y = element_blank(),"
    print "        plot.title = element_text(hjust = 0.5, color = '#666666')) "

def draw_peaks_dis_figures(samples):
    motif_files = []
    for x in xrange(0, 5):
        motif_file = PROJECT_DIR + '/homer/{}_{}/clam_unique/homerResults/motif{}.motif'.format(samples[0], samples[1], x + 1)
        motif_files.append(motif_file)
    peak = PROJECT_DIR + '/clam/peaks/peaks-{}_{}/narrow_peak.combined.bed'.format(samples[0], samples[1])
    dist_data = PROJECT_DIR + '/topology/{}_{}/clam_unique/dist.data'.format(samples[0], samples[1])
    title = '-'.join(samples[0].split('-')[:-1])
    gtf_bed_dir = '{}/peakComposition/{}'.format(cur_dir, GENOME)
    get_motif_rScript(motif_files, title, 1)
    get_dist_topology(dist_data, 2)
    outfn = '{}{}data.pie.txt'.format(TMP_DIR, title)
    #intersect_gtf_regions(peak, outfn, gtf_bed_dir)
    #get_intersection_pie(outfn, 3)
    pdf_output = OUT_DIR + '{}.pdf'.format(title)
    print 'library(ggpubr)'
    #print 'pdf(file="{}", width=7, height=11)'.format(pdf_output)
    print 'pdf(file="{}", width=7, height=6)'.format(pdf_output)
    print 'ggarrange(p1, '
    print '          p2,'
    #print '          p3,'
    #print '          heights = c(1, 1.3, 2.5),'
    #print '          ncol = 1, nrow = 3, align = \'v\')'
    print '          heights = c(1, 2.1),'
    print '          ncol = 1, nrow = 2)'
    print 'dev.off()'
    return pdf_output

def draw_venn_from_sample(out_dir):
    dict_replicate = defaultdict(lambda: [])
    dict_sample_bed_files = defaultdict(lambda: [])
    pdf_output_files = []
    for samples in COMPARISON_LIST:
        dict_replicate['-'.join(samples[0].split('-')[:-2])].append(samples)
    for sample in dict_replicate:
        peak_files = []
        samples = []
        for replicate in dict_replicate[sample]:
            peak_file = PROJECT_DIR + '/clam/peaks/peaks-{}_{}/narrow_peak.unique.bed'.format(replicate[0], replicate[1])
            peak_files.append(peak_file)
            samples.append('-'.join(replicate[0].split('-')[:-1]))
            dict_sample_bed_files[sample].append(peak_file)
        if len(samples) == 2:
            venn_pairwise(peak_files, samples, out_dir + sample + '.venn.pdf')
            pdf_output_files.append(out_dir + sample + '.venn.pdf')
        elif len(samples) == 3:
            venn_triple(peak_files, samples, out_dir + sample + '.venn.pdf')
            pdf_output_files.append(out_dir + sample + '.venn.pdf')
        elif len(samples) == 4:
            vennPlot_quad(peak_files, samples, out_dir + sample + '.venn.pdf')
            pdf_output_files.append(out_dir + sample + '.venn.pdf')
        elif len(samples) == 5:
            vennPlot_quintuple(peak_files, samples, out_dir + sample + '.venn.pdf')
            pdf_output_files.append(out_dir + sample + '.venn.pdf')
        else:
            print 'cannot draw sample {}, len {}'.format(sample, len(samples))    

    samples = dict_replicate.keys()
    samples_len = len(samples)
    for element in combinations([i for i in range(samples_len)], 2):
        sample_1 = samples[element[0]]
        sample_2 = samples[element[1]]
        merged_sample_1 = merged_bed(dict_sample_bed_files[sample_1], '/tmp/peak.sample1.bed')
        merged_sample_2 = merged_bed(dict_sample_bed_files[sample_2], '/tmp/peak.sample2.bed')
        merged_samples = merged_bed(['/tmp/peak.sample1.bed', '/tmp/peak.sample2.bed'], '/tmp/merged.bed')
        venn_pairwise(['/tmp/peak.sample1.bed', '/tmp/peak.sample2.bed'], [sample_1, sample_2], out_dir + sample_1 + '_' + sample_2 + '.venn.pdf')
        pdf_output_files.append(out_dir + sample_1 + '_' + sample_2 + '.venn.pdf')
    return pdf_output_files

def get_pdf_from_sample():
    rscript_out = TMP_DIR  + 'generate_figure.R'
    pdf_output_files = []
    sys.stdout = open(rscript_out, "w+")
    star_mapping_stats()
    table_generate(OUT_DIR + 'mapping_stats.counts.summary.txt', OUT_DIR + 'mapping.stats.table.pdf')
    pdf_output_files.append(OUT_DIR + 'mapping.stats.table.pdf')
    peaks_num_stats()
    plot_peaks_num(OUT_DIR + 'peaks.stats.pdf')
    pdf_output_files.append(OUT_DIR + 'peaks.stats.pdf')
    dict_comparison = {}
    for samples in COMPARISON_LIST:
        sample_name = '-'.join(samples[0].split('-')[:-1])
        dict_comparison[sample_name] = samples
    for sample_name in sorted(dict_comparison.keys()):
        samples = dict_comparison[sample_name]
        pdf_output = draw_peaks_dis_figures(samples)
        pdf_output_files.append(pdf_output)
    pdf_output_files.extend(draw_venn_from_sample(OUT_DIR))
    sys.stdout.flush()
    sys.stdout.close()
    call('Rscript {}'.format(rscript_out), shell=True)
    sys.stdout = sys.__stdout__
    merger = PdfFileMerger()
    for filename in pdf_output_files:
        merger.append(PdfFileReader(file(filename, 'rb')))
    merger.write(OUT_DIR + "report_summary.pdf")

def make_png_output():
    for sample in SAMPLE_MAPPING_DICT:
        for replicate in [1, 2, 3]:
            sample_name = sample + '-' + str(replicate)
            title = SAMPLE_MAPPING_DICT[sample_name[:-2]] + sample_name[-2:]
            print title
            OUT_DIR = LOCAL_PROJECT_DIR + 'summary/'
            pdf_output = OUT_DIR + '{}.pdf'.format(title)
            cmd = 'convert -density 150 {} -quality 100 {}{}.png'.format(pdf_output, OUT_DIR, title)
            os.system(cmd)

def copy_file(src, dist):
    call('cp {} {}'.format(src, dist), shell = True)

def intersect_bed(bed1, bed2, out):
    cmd = 'bedtools intersect -a ' + bed1 +  ' -b ' + bed2 + ' -s | awk \'{if(($3 - $2) > 50) print($_)}\' | sort | uniq > ' + out
    call(cmd, shell = True)

def merged_bed(peak_files, out_file):
    tmp_out = '/tmp/{}.tmp'.format(out_file.split('/')[-1])
    copy_file(peak_files[0], tmp_out)
    for peaks in peak_files[1:]:
        intersect_bed(tmp_out, peaks, out_file)
        copy_file(out_file, tmp_out)

def color_scheme(col_num):
    return ', '.join(['"dodgerblue"', '"goldenrod1"', '"darkorange1"', '"seagreen3"', '"orchid3"'][0:col_num])

def venn_intersect_call(peak_files, samples):
    for index, peaks in enumerate(peak_files):
        print 'area{} = {},'.format(index + 1, len(read_file(peaks)))
    file_num = len(peak_files)
    for x in xrange(1, len(peak_files)):
        for element in combinations([i + 1 for i in range(file_num)], x + 1):
            merged_peaks = [peak_files[_e - 1] for _e in element]
            merged_bed(merged_peaks, '/tmp/peaks.out')
            num_peaks = len(read_file('/tmp/peaks.out'))
            if file_num != 2:
                print 'n{} = {},'.format(''.join([str(x) for x in element]), num_peaks)
            else:
                print 'cross.area = {},'.format(num_peaks)
    print 'category = c({}),'.format(', '.join(['"{}"'.format(sample) for sample in samples]))
    print 'fill = c({}),'.format(color_scheme(file_num))
    print 'cat.col = c({}),'.format(color_scheme(file_num))

def venn_pairwise(peak_files, samples, out_file):
    print 'grid.newpage()'
    print 'library(VennDiagram)'
    print 'pdf("{}", height = 11, width = 11)'.format(out_file)
    print 'venn.plot <- draw.pairwise.venn('
    venn_intersect_call(peak_files, samples)
    print '    lty = "blank",'
    print '    cex = 2,'
    print '    cat.cex = 1.5,'
    print '    cat.pos = c(-10, 0),'
    print '    ext.line.lty = "dashed"'
    print '    )'
    print 'grid.draw(venn.plot)'
    print 'dev.off()'

def venn_triple(peak_files, samples, out_file):
    print 'grid.newpage()'
    print 'pdf("{}", height = 11, width = 11)'.format(out_file)
    print 'library(VennDiagram)'
    print 'venn.plot <- draw.triple.venn('
    venn_intersect_call(peak_files, samples)
    print 'lty = "blank",'
    print 'cex = 2,'
    print 'cat.cex = 1.5,'
    print 'margin = 0.05'
    print ')'
    print 'grid.draw(venn.plot)'
    print 'dev.off()'

def vennPlot_quad(peak_files, samples, out_file):
    print 'grid.newpage()'
    print 'pdf("{}", height = 11, width = 11)'.format(out_file)
    print 'library(VennDiagram)'
    print 'venn.plot <- draw.quad.venn('
    venn_intersect_call(peak_files, samples)
    print 'lty = "blank",'
    print 'cex = 2,'
    print 'cat.cex = 1.5,'
    print 'margin = 0.05'
    print ')'
    print 'grid.draw(venn.plot)'
    print 'dev.off()'

def vennPlot_quintuple(peak_files, samples, out_file):
    print 'grid.newpage()'
    print 'pdf("{}", height = 11, width = 11)'.format(out_file)
    print 'library(VennDiagram)'
    print 'venn.plot <- draw.quintuple.venn('
    venn_intersect_call(peak_files, samples)
    print 'cat.cex = 1.5,'
    print 'margin = 0.5,'
    print 'cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, '
    print '        1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),'
    print 'ind = TRUE)'
    print 'grid.draw(venn.plot)'
    print 'dev.off()'

def get_gene_expression(gtf_path):
    transcript2gene_id = {}
    gtf_annotation = open('{}gtf/gencode.vM13.annotation.gtf'.format(LOCAL_PROJECT_DIR))
    transcript2gene_id = {}
    for gtf in gtf_annotation:
        if gtf.startswith('#'):
            continue
        sp = gtf.strip().split('\t')
        gene_id = re.sub('.*gene_id "|\".*', '', sp[8])
        if sp[2] == 'transcript':
            transcript_id = re.sub('.*transcript_id "|\".*', '', sp[8])
            gene_id = gene_id.split('.')[0]
            transcript_id = transcript_id.split('.')[0]
            transcript2gene_id[transcript_id] = gene_id
    gtf_annotation.close()
    sample_name_list = []
    transcript_exp_dict = defaultdict(lambda:[0] * len(SAMPLE_MAPPING_DICT) * 3)
    sample_index = -1
    for sample in SAMPLE_KEY_LIST:
        for replicate in [1, 2, 3]:
            sample_index += 1
            print 'getting gene exp profile from sample {}-{}'.format(sample, replicate)
            kallisto_input = read_file(LOCAL_PROJECT_DIR + 'kallisto/{}-{}-Inp/abundance.tsv'.format(sample, replicate))
            sample_name_list.append(SAMPLE_MAPPING_DICT[sample] + '-{}'.format(replicate))
            for _exp in kallisto_input[1:]:
                sp = _exp.strip().split('\t')
                transcript_exp_dict[sp[0]][sample_index] = float(sp[4])
    dict_gene_exp = defaultdict(lambda:np.array([0.0 for x in xrange(len(sample_name_list))]))
    for transcript_id in transcript_exp_dict:
        transcript_exp = transcript_exp_dict[transcript_id]
        if transcript_id not in transcript2gene_id:
            continue
        dict_gene_exp[transcript2gene_id[transcript_id]] += np.array([float(x) for x in transcript_exp])
    fw = open(LOCAL_PROJECT_DIR + 'summary/gene_exp.txt', 'w')
    fw.write('gene_id\t{}\n'.format('\t'.join(sample_name_list)))
    for gene_id in dict_gene_exp:
        fw.write('{}\t{}\n'.format(gene_id, '\t'.join([str(x) for x in dict_gene_exp[gene_id]])))
    fw.close()

def get_top_gene_selected():
    gene_exp = read_file(LOCAL_PROJECT_DIR + 'summary/gene_exp.txt')
    fw = open(LOCAL_PROJECT_DIR + 'tmp/gene_exp_top.txt', 'w')
    header =gene_exp[0]
    fw.write(header + '\n')
    count = 0
    for line in gene_exp[1:]:
            sp = line.strip().split('\t')
            gene_list = [float(x) for x in sp[1:]]
            if np.sum(gene_list) == 0:
                    continue
            gene_list = zscore(gene_list)
            if np.max(gene_list) - np.min(gene_list) > 3.7:
                    fw.write(sp[0] + '\t' + '\t'.join([str(_value) for _value in gene_list]) + '\n')
                    count += 1
    fw.close()
    print count

get_pdf_from_sample()

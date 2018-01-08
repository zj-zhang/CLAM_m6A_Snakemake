import sys
import pysam

aligner_out = sys.argv[1]
is_stranded = sys.argv[2]
is_stranded = True if is_stranded.lower().startswith('t') else False ## default: is_stranded=F
strand = sys.argv[3]
is_multi = False if len(sys.argv)>4 and sys.argv[4].lower().startswith('f') else True ## default: is_multi=T
is_reverse = True if strand=='-' else False

bam = pysam.Samfile(aligner_out, 'rb')
chr_list = chr_list=[x['SN'] for x in bam.header['SQ']]
for read in bam:
	this_score = str(read.opt('AS')) if is_multi else '1.0'
	if this_score=='0.0':
		continue
	this_strand = '-' if read.is_reverse else '+'
	this_qlen = sum( [ x[1] for x in read.cigar if x[0]==0 ] )
	line = "\t".join(
		[chr_list[read.rname], 
		str(read.positions[0]), 
		str(read.positions[0]+this_qlen), 
		this_strand,
		this_score])
	if is_stranded:
		if read.is_reverse==is_reverse:
			print >> sys.stdout, line
		else:
			continue
	else:
		print >> sys.stdout, line
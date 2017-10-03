#!/usr/bin/env python2
import sys
import os

cur_dir = os.path.dirname(os.path.realpath(__file__))
binnum=50
peakfile=sys.argv[1]
genome = sys.argv[2] if len(sys.argv)>2 else 'hg19'
filename=peakfile.split('/')[-1]

def exon_cds (start, end, cdsstart,cdsend):
            for i in range(len(start)):
                if (start[i]<=cdsstart<end[i]):
                    m=i
                elif(cdsstart==end[i]):
                    #print ['cdsstart==end[i]', start, end, cdsstart,cdsend]
                    m=i+1
                if (start[i]<=cdsend<=end[i]):
                    n=i
                    break
            #print['exon_cds', start, end, cdsstart,cdsend]
            changedstart=[cdsstart]+start[m+1:n+1]
            changedend=end[m:n]+[cdsend]
            utr5start=start[:m+1]
            utr5end=end[:m]+[cdsstart]
            utr3start=[cdsend]+start[n+1:]
            utr3end=end[n:]
            cdslength=0
            utr5length=0
            utr3length=0
            for i in range(len(changedstart)):
                cdslength+=changedend[i]-changedstart[i]
            for i in range(len(utr5start)):
                utr5length+=utr5end[i]-utr5start[i]
            for i in range(len(utr3start)):
                utr3length+=utr3end[i]-utr3start[i]
            #if (cdsstart==152789537):
                #print [changedstart, changedend,cdslength, utr5start, utr5end, utr5length,utr3start, utr3end, utr3length]
            return  [changedstart, changedend,cdslength, utr5start, utr5end, utr5length,utr3start, utr3end, utr3length]

def bin_tx (mid, exonstart, exonend, strand, binnum):
    txlen=0
    pos=0
    for i in range(len(exonstart)):
        if (exonstart[i]<=mid<exonend[i]):
            pos=txlen+mid-exonstart[i]
        txlen+=exonend[i]-exonstart[i]
    bin=pos*binnum/(txlen)
    if (strand=='-'):
        bin=binnum-bin-1
    elif (strand=='+'):
        pass
    else:
        print 'strand error'
    return bin

longtx={}
input=open(os.path.join(cur_dir,'%s_ensGene.txt'%genome))
for line in input:
    line=line[:-1]
    a=line.split('\t')
    tx=a[1]
    chr=a[2]
    strand=a[3]
    exonstart=a[9].split(',')[:-1]
    exonend=a[10].split(',')[:-1]
    ensg=a[12]
    length=0
    for i in range(len(exonstart)):
        length+=int(exonend[i])-int(exonstart[i])
    if (not longtx.has_key(ensg)):
        longtx[ensg]=[tx, length, a]
    elif(length> longtx[ensg][1]):
        longtx[ensg]=[tx, length, a]
    else:
        pass
input.close()

peak={}
for line in open(peakfile):
    line=line[:-1]
    a=line.split('\t')
    """
    if(a[0]=='allid'):
        continue
    chr, peakstart, peakend=a[0].split(':')
    peakstart=int(peakstart)
    peakend=int(peakend)
    """
    if(a[0]=='Chr'):
        continue
    chr=a[0]
    peakstart=int(a[1])
    peakend=int(a[2])
    mid=(peakstart+peakend)/2
    if(peak.has_key(chr)):
        peak[chr].append(mid)
    else:
        peak[chr]=[mid]
        
coding={}
for line in open(os.path.join(cur_dir,'%s_ensg_genename.txt'%genome)):
    line=line[:-1]
    a=line.split()
    if a[0].startswith('#'):
        continue
    coding[a[0]]=a[1]


peakbin=[0]*binnum
peakbin5=[0]*binnum
peakbin3=[0]*binnum
peakbincds=[0]*binnum
txnum=0
tx_with_peak={}  #gene have at least one peak
#for line in open('/u/nobackup/yxing/jinkwang_backup/m6A/level_result/hg19_ensGene.txt'):
for line in open(os.path.join(cur_dir,'%s_ensGene.txt'%genome)):
    line=line[:-1]
    a=line.split('\t')
    if (a[0]=='#bin'):
        continue
    tx=a[1]
    chr=a[2]
    strand=a[3]
    exonstart=a[9].split(',')[:-1]
    exonend=a[10].split(',')[:-1]
    genename=a[12]
    cdsstart=int(a[6])
    cdsend=int(a[7])
    if (longtx[genename][0]!=tx or cdsstart==cdsend):
        continue
    #if (not coding.has_key(genename) or coding[genename]!='protein_coding'):
    if (not coding.has_key(tx) or coding[tx]!='protein_coding'):
        continue
    if (not peak.has_key(chr)):
        continue
    for i in range(len(exonstart)):
        exonstart[i]=int(exonstart[i])
        exonend[i]=int(exonend[i])
    change=exon_cds(exonstart, exonend, cdsstart, cdsend)
    if (change[2]==0 or change[5]==0 or change[8]==0):
        continue
    txnum+=1
    for i in range(len(exonstart)):
        for mid in peak[chr]:
            if (exonstart[i]<=mid<exonend[i]):
                bin=bin_tx(mid, exonstart, exonend, strand, binnum)
                peakbin[bin]+=1
                if (strand=='+'):
                    if (mid<cdsstart):
                        bin5=bin_tx(mid, change[3], change[4], strand, binnum)
                        peakbin5[bin5]+=1
                        #if (bin5==0):
                            #print [chr, mid, strand]
                    elif(mid>=cdsend):
                        bin3=bin_tx(mid, change[6], change[7], strand, binnum)
                        peakbin3[bin3]+=1
                    elif(cdsstart<=mid<cdsend):
                        bincds=bin_tx(mid, change[0], change[1], strand, binnum)
                        peakbincds[bincds]+=1
                        #if(bincds==0):
                        #   print [chr, mid, strand]
                    else:
                        print 'utr error'
                elif(strand=='-'):
                    if (mid<cdsstart):
                        bin3=bin_tx(mid, change[3], change[4], strand, binnum)
                        peakbin3[bin3]+=1
                    elif(mid>=cdsend):
                        bin5=bin_tx(mid, change[6], change[7], strand, binnum)
                        peakbin5[bin5]+=1
                        #if (bin5==0):
                            #print [chr, mid, strand]
                    elif(cdsstart<=mid<cdsend):
                        bincds=bin_tx(mid, change[0], change[1], strand, binnum)
                        peakbincds[bincds]+=1
                    else:
                        print 'utr error'
                else:
                    print "strand error"
                if (tx_with_peak.has_key(genename)):
                    tx_with_peak[genename]+=1
                else:
                    tx_with_peak[genename]=1
                    
print >>sys.stderr, ['type', 'genenum', txnum, 'gene_with_peak', len(tx_with_peak.keys()), 'peaknum', len([x for chr in peak for x in peak[chr]])]



output = sys.stdout
for i in range(len(peakbin)):
    output.write(str(peakbin5[i]*1.0/txnum)+'\t'+str(peakbin5[i]*1.0/len(tx_with_peak.keys()))+'\t'+str(peakbin5[i])+'\n')
for i in range(len(peakbin)):
    output.write(str(peakbincds[i]*1.0/txnum)+'\t'+str(peakbincds[i]*1.0/len(tx_with_peak.keys()))+'\t'+str(peakbincds[i])+'\n')
for i in range(len(peakbin)):
    output.write(str(peakbin3[i]*1.0/txnum)+'\t'+str(peakbin3[i]*1.0/len(tx_with_peak.keys()))+'\t'+str(peakbin3[i])+'\n')


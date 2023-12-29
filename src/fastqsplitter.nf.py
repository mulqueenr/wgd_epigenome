import gzip
from Bio import SeqIO
import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pandas as pd
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--fq1')
parser.add_argument('--fq2')
parser.add_argument('--idx3')
parser.add_argument('--idx4')
parser.add_argument('--platform',default="echo")
parser.add_argument('--plate',default="2")


args = parser.parse_args()

#fq1="Undetermined_S0_L001_R1_001.chunk050.fastq.gz"
#fq2="Undetermined_S0_L001_R2_001.chunk050.fastq.gz"
#idx3="Undetermined_S0_L001_I1_001.chunk050.fastq.gz"
#idx4="Undetermined_S0_L001_I2_001.chunk050.fastq.gz"
#platform="echo"
#plate="1,2,3,4"
fq1=args.fq1
fq2=args.fq2
idx3=args.idx3
idx4=args.idx4
platform=args.platform

#accept comma separated list for plate args
if ',' in args.plate:
        plate=[args.plate.split(',')]
else:
        plate=[args.plate]

#set up plate
if platform == "icell8":
        n7="/volumes/seq/projects/wgd/ref/wafer_v2_N7_72.txt" #copied to home directory as backup
        s5="/volumes/seq/projects/wgd/ref/wafer_v2_S5_72.txt" #copied to home directory as backup
        n7=pd.read_table(n7)
        s5=pd.read_table(s5)
        plate_i7=n7["RC_N7"]
        plate_i5=s5["RC_S5_Hiseq4000_nextseq"]
        plate_i7=[i.strip() for i in plate_i7]
        plate_i5=[i.strip() for i in plate_i5]
elif platform == "echo":
        idx="/volumes/seq/projects/wgd/ref/echo_idx.tsv"
        idx=pd.read_table(idx)
        plate_i7=[]
        plate_i5=[]
        for plate_number in plate:
                if plate_number == "1":
                        plate_i7.append([i.strip() for i in idx.loc[(idx['idx'] == 'i7') & (idx['set'] == 1)]['idx_seq']])
                        plate_i5.append([i.strip() for i in idx.loc[(idx['idx'] == 'i5') & (idx['set'] == 1)]['idx_seq']])
                if plate_number == "2":
                        plate_i7.append([i.strip() for i in idx.loc[(idx['idx'] == 'i7') & (idx['set'] == 2)]['idx_seq']])
                        plate_i5.append([i.strip() for i in idx.loc[(idx['idx'] == 'i5') & (idx['set'] == 1)]['idx_seq']])
                if plate_number == "3":
                        plate_i7.append([i.strip() for i in idx.loc[(idx['idx'] == 'i7') & (idx['set'] == 1)]['idx_seq']])
                        plate_i5.append([i.strip() for i in idx.loc[(idx['idx'] == 'i5') & (idx['set'] == 2)]['idx_seq']])
                if plate_number == "4":
                        plate_i7.append([i.strip() for i in idx.loc[(idx['idx'] == 'i7') & (idx['set'] == 2)]['idx_seq']])
                        plate_i5.append([i.strip() for i in idx.loc[(idx['idx'] == 'i5') & (idx['set'] == 2)]['idx_seq']])
        plate_i7=[k for i in plate_i7 for k in i]
        plate_i5=[k for i in plate_i5 for k in i]
else:
        print("Please use a valid platform and plate combination for demultiplexing.")

i=0
#open fastq files, correct barcode read names then out fastq 1 and 2  with new read name
with gzip.open(fq1, "rt") as handle1:
        with gzip.open(fq2, "rt") as handle2:
                with gzip.open(idx3, "rt") as handle3:
                        with gzip.open(idx4, "rt") as handle4:
                                with open(fq1[:-9]+".barc.fastq", "w") as outfile_fq1:
                                        with open(fq2[:-9]+".barc.fastq", "w") as outfile_fq2:
                                                for (title1, seq1, qual1), (title2, seq2, qual2), \
                                                (title3,seq3,qual3), \
                                                (title4,seq4,qual4) in \
                                                zip(FastqGeneralIterator(handle1), \
                                                        FastqGeneralIterator(handle2),\
                                                        FastqGeneralIterator(handle3),\
                                                        FastqGeneralIterator(handle4)):
                                                        if seq3[:8] in plate_i7 and seq4[:8] in plate_i5:
                                                                i+=1
                                                                readname=seq3[:8]+seq4[:8]
                                                                outfile_fq1.write("@%s:%s\n%s\n+\n%s\n" % (readname, i, seq1, qual1))
                                                                outfile_fq2.write("@%s:%s\n%s\n+\n%s\n" % (readname, i, seq2, qual2))




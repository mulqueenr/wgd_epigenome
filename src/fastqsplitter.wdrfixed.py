import gzip
from Bio import SeqIO
import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
import pandas as pd
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--fq1')
parser.add_argument('--fq2')
parser.add_argument('--idx3')
parser.add_argument('--idx4')
parser.add_argument('--samples')
parser.add_argument('--sample_layout')
parser.add_argument('--dna_chip_primer',default="MD_N708_out")
parser.add_argument('--rna_chip_primer',default="DDR_PCR_p5_UDI8_V2")

args = parser.parse_args()


"""
Break down of indexes in a Nextseq sequencing run of --use-bases-mask=Y46,I34,I8,Y50
zless Undetermined_S0_I1_001.fastq.gz | grep -E "^[A|T|C|G|N]+$" | cut -c1-8 | sort | uniq -c | sort -k1,1n
        #DNA MD_atacN7##_in REVCOMP (8bp) #chip column
        #RNA ArcDR_RNA_RPI_atacN714 REVCOMP (8bp) #chip column
zless Undetermined_S0_I1_001.fastq.gz | grep -E "^[A|T|C|G|N]+$" | cut -c9-27 | sort | uniq -c | sort -k1,1n
        #DNA i7 Constant region REVCOMP (18bp)
        #RNA Ignore
zless Undetermined_S0_I1_001.fastq.gz | grep -E "^[A|T|C|G|N]+$" | cut -c27- | sort | uniq -c | sort -k1,1n
        #DNA i7 MD_N##_out (8bp) #chip pool
        #RNA Ignore
zless Undetermined_S0_I2_001.fastq.gz | grep -E "^[A|T|C|G|N]+$" | cut -c1-8 | sort | uniq -c | sort -k1,1n
        #DNA i5 S5##/atac/Nflex (8bp) #chip row
        #RNA DDR_PCR_p5_UDI#_v2 (8bp) #chip pool
zless Undetermined_S0_R1_001.fastq.gz | grep -E "^[A|T|C|G|N]+$" | cut -c1-8 | sort | uniq -c | sort -k1,1n
        #DNA Ignore
        #RNA atacS## i5 (8bp) #chip row
"""


# TO ADD IN FUTURE USE CASES #If no samplesheet given, just try everything,
# TO ADD IN FUTURE USE CASES #Set up for proper reverse complementing for different sequencers
#DNA                
DNA_i7_chip_col="/ref/wafer_v2_multichip_N7_72.txt" #REVCOMP of 1-8bp in I1
DNA_i7_chip_col=pd.read_table(DNA_i7_chip_col)
DNA_i7_chip_col.columns=[x+"_DNA" for x in DNA_i7_chip_col.columns]

DNA_i7_constant="CTGAGTCGGAGACACGCA" #REVCOMP of 9-27bp in I1
DNA_chip_pool="/ref/wafer_v2_multichip_N7chipbc_11.txt" #REVCOMP of 27-34bp in I1
DNA_chip_pool=pd.read_table(DNA_chip_pool)
DNA_chip_pool.columns=[x+"_DNA" for x in DNA_chip_pool.columns]

DNA_i5_chip_row="/ref/wafer_v2_S5_72.txt" #1-8bp in I2
DNA_i5_chip_row=pd.read_table(DNA_i5_chip_row)
DNA_i5_chip_row.columns=[x+"_DNA" for x in DNA_i5_chip_row.columns]

#RNA
RNA_i7_chip_col="/ref/arcdr_RNA_i7_72.tsv" #REVCOMP of 1-8bp in I1
RNA_i7_chip_col=pd.read_table(RNA_i7_chip_col) 
RNA_i7_chip_col.columns=[x+"_RNA" for x in RNA_i7_chip_col.columns]

RNA_chip_pool="/ref/wafer_v2_multichip_i5_RNA_chipbc_9.txt" #1-8bp in I2
RNA_chip_pool=pd.read_table(RNA_chip_pool)
RNA_chip_pool.columns=[x+"_RNA" for x in RNA_chip_pool.columns]

RNA_i5_chip_row="/ref/arcdr_RNA_i5_72.tsv" #1-8bp in R2
RNA_i5_chip_row=pd.read_table(RNA_i5_chip_row)
RNA_i5_chip_row.columns=[x+"_RNA" for x in RNA_i5_chip_row.columns]

#Make full DataFrame
sample=pd.read_table(args.samples)

#Merge sample layout to metadata
sample_layout=pd.read_table(args.sample_layout,sep="\t")
sample=pd.merge(sample,sample_layout,left_on='Source well',right_on='Well')
sample=pd.merge(sample,DNA_i5_chip_row,left_on='Row',right_on='Row_DNA')
sample=pd.merge(sample,DNA_i7_chip_col,left_on='Col',right_on='Col_DNA',suffixes=("_i5","_i7"))
sample=pd.merge(sample,RNA_i5_chip_row,left_on='Row',right_on='Row_RNA')
sample=pd.merge(sample,RNA_i7_chip_col,left_on='Col',right_on='Col_RNA',suffixes=("_i5","_i7"))
sample['cellID']=["C_"+'_'.join(i) for i in zip(sample['Col'].map(str),sample['Row'].map(str))] #define cells by Col x Row

#Set chip indexes
sample['idx_chip_DNA']=DNA_chip_pool.loc[DNA_chip_pool['idx_name_DNA']==args.dna_chip_primer]['pool_idx_DNA'].values[0]
sample['idx_chip_RNA']=RNA_chip_pool.loc[RNA_chip_pool['idx_name_RNA']==args.rna_chip_primer]['pool_idx_RNA'].values[0]

#Perform Reverse Complementing
sample['col_idx_DNA']=[str(Seq(x).reverse_complement()) for x in sample['col_idx_DNA']] #REVCOMP of 1-8bp in I1
sample['row_idx_DNA']=[str(Seq(x).reverse_complement()) for x in sample['row_idx_DNA']] #REVCOMP of 1-8bp in I2
sample['idx_chip_DNA']=[str(Seq(x).reverse_complement()) for x in sample['idx_chip_DNA']] #REVCOMP of 26-34bp in I1
sample['col_idx_RNA']=[str(Seq(x).reverse_complement()) for x in sample['col_idx_RNA']] #REVCOMP of 1-8bp in I1

fq1=args.fq1
fq2=args.fq2
idx3=args.idx3
idx4=args.idx4

i=0

#open fastq files, correct barcode read names then out fastq 1 and 2  with new read name
with gzip.open(fq1, "rt") as handle1, \
     gzip.open(fq2, "rt") as handle2, \
     gzip.open(idx3, "rt") as handle3, \
     gzip.open(idx4, "rt") as handle4, \
     open(fq1[:-9]+".DNA.barc.fastq", "w") as outfile_DNA_fq1, \
     open(fq2[:-9]+".DNA.barc.fastq", "w") as outfile_DNA_fq2, \
     open(fq1[:-9]+".RNA.barc.fastq", "w") as outfile_RNA_fq1, \
     open(fq2[:-9]+".RNA.barc.fastq", "w") as outfile_RNA_fq2:
        for (title1, seq1, qual1), \
        (title2, seq2, qual2), \
        (title3,seq3,qual3), \
        (title4,seq4,qual4) in \
        zip(FastqGeneralIterator(handle1), \
                FastqGeneralIterator(handle2),\
                FastqGeneralIterator(handle3),\
                FastqGeneralIterator(handle4)):
                if any(sample.idx_chip_DNA==seq3[26:34]): #treat read as DNA
                        dna_col=sample.index[sample.col_idx_DNA == seq3[:8]].tolist()
                        dna_row=sample.index[sample.row_idx_DNA == seq4[:8]].tolist()
                        sample_idx=list(set(dna_col).intersection(dna_row))
                        if len(sample_idx)==1:
                                sample_idx=sample_idx[0]
                                read_name="%s:%s_%s_%s" % (sample.loc[sample_idx,'cellID'], sample.loc[sample_idx,'idx_chip_DNA'], sample.loc[sample_idx,'col_idx_DNA'],sample.loc[sample_idx,'row_idx_DNA'])
                                outfile_DNA_fq1.write("@%s:%s\n%s\n+\n%s\n" % (read_name, i, seq1, qual1))
                                outfile_DNA_fq2.write("@%s:%s\n%s\n+\n%s\n" % (read_name, i, seq2, qual2))
                if any(sample.idx_chip_RNA == seq4[:8]): #treat read as RNA
                        rna_col=sample.index[sample.col_idx_RNA == seq3[:8]].tolist()
                        rna_row=sample.index[sample.row_idx_RNA == seq1[:8]].tolist()
                        sample_idx=list(set(rna_col).intersection(rna_row))
                        if len(sample_idx)==1:
                                sample_idx=sample_idx[0]
                                read_name="%s:%s_%s_%s" % (sample.loc[sample_idx,'cellID'], sample.loc[sample_idx,'idx_chip_RNA'], sample.loc[sample_idx,'col_idx_RNA'],sample.loc[sample_idx,'row_idx_RNA'])
                                outfile_RNA_fq1.write("@%s:%s\n%s\n+\n%s\n" % (read_name, i, seq1, qual1))
                                outfile_RNA_fq2.write("@%s:%s\n%s\n+\n%s\n" % (read_name, i, seq2, qual2))



sample.to_csv('metadata.csv')


/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/wgd/240505_RNADNA_WGD_t0/240322_ARCDR_FIXED_WGDt0_141014C_scan2/141014-manual3_reprocess_WellList.TXT



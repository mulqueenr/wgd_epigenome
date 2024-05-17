
#set up environment variables 
export projDir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/wgd"
export sif="/home/rmulqueen/singularity/scmetR.sif"

module load singularity 
singularity shell \
--bind ${projDir} \
$sif

library(optparse)
library(Biostrings)
library(reshape2)
library(ggplot2)
option_list = list(
  make_option(c("-f", "--filter_file"), type="character", default=NULL, 
              help="*FilterFile.csv filter file from Takara icell8."),
  make_option(c("-d", "--data_file"), type="character", default=NULL, 
              help="*WellList.TXT data file from Takara icell8."),
  make_option(c("-s", "--sample_layout"), type="character", default=NULL, 
              help="CSV file with sample names per well during icell8 dispense."),
  make_option(c("-b", "--wafar_bc"), type="character", default="/volumes/seq/projects/wgd/src/Wafar_cell.csv", 
              help="Barcode file from standard wafergen/icell8 barcodes setup."),
  make_option(c("-c", "--wdr_bc"), type="character", default="/volumes/seq/projects/wgd/src/WDR_FIX_cellBC.txt", 
              help="Barcode file from WDR barcodes setup.")

);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

data_file=read.delim("/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/wgd/240505_RNADNA_WGD_t0/240322_ARCDR_FIXED_WGDt0_141014C_scan2/141014_WellList.TXT",header = T)
sample_layout<-read.delim("/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/wgd/240505_RNADNA_WGD_t0/sample_layout.txt",sep="\t",header=T)

#Add RNA Indexes
rna_i5_in=read.delim("/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/wgd/ref/arcdr_RNA_i5_72.tsv",header = F)
rna_i5_in<-rna_i5_in[,c("V1","V3")];colnames(rna_i5_in)<-c("rna_i5in_idx_name","rna_i5in_idx")
rna_i5_in$Row<-1:72

rna_i7=read.delim("/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/wgd/ref/arcdr_RNA_i7_72.tsv",header = F)
rna_i7<-rna_i7[,c("V1","V3")];colnames(rna_i7)<-c("rna_i7_idx_name","rna_i7_idx")
rna_i7$Col<-1:72

data_file<-merge(data_file,rna_i5_in,by="Row")
data_file<-merge(data_file,rna_i7,by="Col")

data_file$rna_i5_out_name<-"DDR_PCR_P5_UDI9_V2"
data_file$rna_i5_out_idx<-"TCTGTTGG"

#Add DNA Indexes
dna_i7_out=read.delim("/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/wgd/ref/wafer_v2_multichip_N7chipbc_11.txt",header = F)
dna_i7_out<-dna_i7_out[,c("V2","V4")];colnames(dna_i7_out)<-c("dna_i7out_idx_name","dna_i7out_idx")

dna_i7_in=read.delim("/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/wgd/ref/wafer_v2_multichip_N7_72.txt",header = F)
dna_i7_in<-dna_i7_in[,c("V2","V4")];colnames(dna_i7_in)<-c("dna_i7in_idx_name","dna_i7in_idx")
dna_i7_in$Col<-1:72

dna_i5=read.delim("/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/wgd/ref/wafer_v2_S5_72.txt",header = T)
dna_i5<-dna_i5[,c("Name","RC_S5_Hiseq4000_nextseq")];colnames(dna_i5)<-c("dna_i5_idx_name","dna_i5_idx")
dna_i5$Row<-1:72

data_file<-merge(data_file,dna_i5,by="Row")
data_file<-merge(data_file,dna_i7_in,by="Col")

data_file$dna_i7out_idx_name<-"MD_OUT_N709"
data_file$dna_i7in_idx<-"AGCGTAGC"

#as.character(complement(DNAStringSet(reverse(data_file$DNA_S5))))
sample_layout<-read.delim("/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/wgd/240505_RNADNA_WGD_t0/sample_layout.txt",header=T,sep="\t",col.names=c("SampleWell","Sample_name"))
data_file<-merge(data_file,sample_layout,by="SampleWell")
data_file$Sample_name<-gsub(" ","_",data_file$Sample_name)

data_file$CellID<-paste0("Cell_",1:nrow(data_file))

#RNA Samplesheet:


rna_out<-data.frame(
  Lane="*",
  Sample_ID=data_file$CellID,
  Sample_Name=data_file$Sample_name,
  Sample_Well=paste0("Row",data_file$Row,"_","Col",data_file$Col),
  index=
  index2=
  )

dna_out<-data.frame(
  Lane="*",
  Sample_ID=data_file$CellID,
  Sample_Name=data_file$Sample_name,
  Sample_Well=paste0("Row",data_file$Row,"_","Col",data_file$Col),
  index=data_file$dna_i7in_idx,
  index2=paste0(data_file$dna_i7out_idx,data_file$dna_i5_idx)
  )
rna<-data_file
  echo "[Header],,,,,,,,
  Investigator Name,Ryan,,,,,,,
  Date,\$(date +'%d/%m/%Y'),,,,,,,
  Workflow,GenerateFASTQ,,,,,,,
  [Data],,,,,,,,
  Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,Sample_Project,Description" > SampleSheet.csv

#bclfiles in 
#/volumes/seq/flowcells/MDA/nextseq2000/2024/20240426_Ryan_231_wgdT0_DNA_RNA_ESO_P17_Chip1_DNA_RNA_reseq_P18_Chip2_DNA/240425_VH00219_582_AACFVYNHV

#DNA by cellTAGATCGC
#ArcDR_ESO_pre_P48_Chip2_DNA_C1  MD_atacN714_in(') - MD_OUT_N709_S501(')(CCGTTTGT - GCTACGCTGCGATCTA)

#RNA by column
#ArcDR_ESO_preTX_P48_RNA_Chip2_Column1   MD_atacN714_in-DDR_PCR_P5_UDI9_V2(CCGTTTGTAT-TCTGTTGG)

#Cell,DNA_Barcode,RNA,RNA_Barcode,DNA_Cell_Name,row,col
#C1,N701-S501(TAAGGCGA-GCGATCTA),WDR_RNA_N702_2-WDR_RNA_TS503_2,CGTACTAGTATCCTCT,TAAGGCGAGCGATCTA,0,0



Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,Sample_Project,Description
,ArcDR_ESO_preTX_P17_C1,ArcDR_ESO_preTX_P17_C1,plate1,A1,CCGTTTGT,TAAGGCGAGCGATCTA,ArcDR_ESO_preTX_P17_C1,ArcDR_ESO_preTX_P17_C1

fq_out="/volumes/seq/projects/wgd/240322_ARCDR_FIXED_WGDt0/fastq"
runDir="/volumes/seq/flowcells/MDA/nextseq2000/2024/20240426_Ryan_231_wgdT0_DNA_RNA_ESO_P17_Chip1_DNA_RNA_reseq_P18_Chip2_DNA/240425_VH00219_582_AACFVYNHV"


#Running RNA
bcl2fastq \
--use-bases-mask=Y46,I8n18,I16,Y50 \
--create-fastq-for-index-reads \
--barcode-mismatches 0 \
-r 60 -p 60 -w 60 \
--no-lane-splitting \
-R $runDir \
-o $fq_out

zcat Undetermined_S0_I1_001.fastq.gz | grep -E '^[ATGC]+$' | sort | uniq -c | sort -k1,1n

#Running DNA
bcl2fastq \
--use-bases-mask=Y46,I34,I8,Y50 \
--create-fastq-for-index-reads \
--barcode-mismatches 0 \
-r 60 -p 60 -w 60 \
--no-lane-splitting \
-R $runDir \
-o $fq_out



//<Read Number="1" NumCycles="46" IsIndexedRead="N" IsReverseComplement="N"/>
//<Read Number="2" NumCycles="34" IsIndexedRead="Y" IsReverseComplement="N"/> #pool barcode+icell8columnidx
//<Read Number="3" NumCycles="16" IsIndexedRead="Y" IsReverseComplement="Y"/> #icell8rowidx
//<Read Number="4" NumCycles="50" IsIndexedRead="N" IsReverseComplement="N"/>


#--sample-sheet=$WPATH/SampleSheet.csv
#dat$RNA_BC is true RNA barcodes for sample demultiplexing
#RNA index2 (i7) MD_atacN714-792[column] [10bp]
#RNA index1 (i5) DDR_PCR_P5_UDI8_V2 [pool] [8bp]

#DNA:GTGATAGC


#set up environment variables 
export projDir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/wgd"
export sif="/home/rmulqueen/singularity/scmetR.sif"

module load singularity 
singularity shell \
--bind ${projDir} \
$sif

/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/wgd/rundir/240425_VH00219_582_AACFVYNHV



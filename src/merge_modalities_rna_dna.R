#bsub -Is -W 6:00 -q medium -n 10 -M 50 -R rusage[mem=16] /bin/bash 
#module load singularity
#proj_dir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/wgd"
#sif="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/wgd/src/copykit.sif"
#singularity shell \
#--bind ${proj_dir}/ref:/ref \
#--bind ${proj_dir}/src:/src \
#--bind /rsrch4/scratch/genetics/rmulqueen/ \
#--bind $proj_dir \
#$sif

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
BiocParallel::bpparam()
library(optparse)

option_list = list(
  make_option(c("-i", "--input_dir"), type="character", default=".", 
              help="Directory location with necessary files.", metavar="character"),
  make_option(c("-o", "--outname"), type="character", default="WGD", 
              help="Outname prefix for files.", metavar="character"),
  make_option(c("-d", "--dna_cnv_metadata"), type="character", default="WGD.scCNA.rds", 
              help="DNA Object with CNV calls.", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
setwd(opt$input_dir)

#read in RNA probes
rna<-ReadMtx(
  mtx="probe_count_matrix.mtx.gz",
  cells="barcodes.tsv.gz",
  features="features.tsv.gz",
  cell.column = 1,feature.column = 2,
  skip.cell=1,skip.feature=1,
  cell.sep = "\t",feature.sep = "\t",
)
#create object
dat <- CreateSeuratObject(counts = rna, project = opt$outname, min.cells = 5, min.features = 500)

#add dna clonality info
dna<-read.csv(opt$dna_cnv_metadata)
row.names(dna)<-dna$cellID

#combine to seurat object
dat<-AddMetaData(dat,dna)

dat <- PercentageFeatureSet(dat, pattern = "^MT-", col.name = "percent.mt")
dat <- subset(dat, subset = nFeature_RNA < 20000 & percent.mt < 30 & reads_assigned_bins > 50000)

dat <- NormalizeData(dat)
dat <- ScaleData(dat)
dat <- FindVariableFeatures(dat,nfeatures=3000)
dat <- RunPCA(dat, verbose = FALSE, features=VariableFeatures(dat))
dat <- FindNeighbors(dat, dims = 1:30)
dat <- FindClusters(dat, resolution = 0.4, verbose = FALSE)
dat <- RunUMAP(dat, dims = 1:30)

# Now, we can include the key in the feature name, which overrides the default assay
p1 <- DimPlot(dat, group.by="Sample_y") + ggtitle("Treatment")
p2 <- FeaturePlot(dat, "Signal1",order=T) + ggtitle("DAPI Signal")

p3 <- DimPlot(dat, group.by="superclones") + ggtitle("DNA Superclones")
p4 <- DimPlot(dat, group.by="subclones") + ggtitle("DNA Subclones")

p5 <- DimPlot(dat, group.by="seurat_clusters") + ggtitle("Seurat Clusters")
p6 <- FeaturePlot(dat, "Circularity1",order=T) + ggtitle("DAPI Circularity")

p7 <- FeaturePlot(dat, "reads_assigned_bins") + ggtitle("DNA Read Count")
p8 <- FeaturePlot(dat, "nCount_RNA") + ggtitle("RNA Read Count")

plt<-(p1|p2)/(p3|p4)/(p5|p6)/(p7|p8)
ggsave(plt,file=paste0(opt$outname,".qc_metrics.umap.pdf"),height=16,width=10)

saveRDS(dat,file=paste0(opt$outname,".SeuratObj.Rds"))

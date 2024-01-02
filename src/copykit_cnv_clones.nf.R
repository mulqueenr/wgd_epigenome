library(copykit)
library(BiocParallel)
library(EnsDb.Hsapiens.v86)


args = commandArgs(trailingOnly=TRUE)
setwd(args[1])
meta_out=args[2]
register(MulticoreParam(progressbar = T, workers = args[3]), default = T)
BiocParallel::bpparam()   
data <- runVarbin(args[1],
                 remove_Y = TRUE,
                 genome="hg38",
                 is_paired_end=TRUE)



# Mark euploid cells if they exist
data <- findAneuploidCells(data)

# Mark low-quality cells for filtering
data <- findOutliers(data)

# Visualize cells labeled by filter and aneuploid status
pdf(paste(meta_out,"outlier_qc.heatmap.pdf",sep="/"))
plotHeatmap(data, label = c('outlier', 'is_aneuploid'), row_split = 'outlier')
dev.off()

# Remove cells marked as low-quality and/or aneuploid from the copykit object
#data <- data[,SummarizedExperiment::colData(data)$outlier == FALSE]
#data <- data[,SummarizedExperiment::colData(data)$is_aneuploid == TRUE]


# kNN smooth profiles
data <- knnSmooth(data)

data <- runUmap(data)
data<-findSuggestedK(data)
# Create a umap embedding 

# Find clusters of similar copy number profiles and plot the results
# If no k_subclones value is provided, automatically detect it from findSuggestedK()
data  <- findClusters(data,k_subclones=metadata(data)$suggestedK)#output from k_clones

pdf(paste(meta_out,"subclone.umap.pdf",sep="/"))
plotUmap(data, label = 'subclones')
dev.off()

# Calculate consensus profiles for each subclone, 
# and order cells by cluster for visualization with plotHeatmap
data <- calcConsensus(data)
data <- runConsensusPhylo(data)

# Plot a copy number heatmap with clustering annotation
pdf(paste(meta_out,"subclone.heatmap.pdf",sep="/"))
plotHeatmap(data, label = 'subclones',order='hclust')
dev.off()

saveRDS(data,file=paste(meta_out,"scCNA.rds",sep="/"))
clone_out<-data.frame(bam=paste0(getwd(),"/",row.names(data@colData),".bam"),clone=data@colData$subclones)
for (i in unique(clone_out$clone)){
	tmp<-clone_out[clone_out$clone==i,]
	write.table(tmp$bam,file=paste0("clone_",i,".bam_list.txt"),row.names=F,col.names=F,quote=F)
}

met<-as.data.frame(data@colData)
write.table(met,file=paste(meta_out,"meta.tsv",sep="/"),quote=F,sep="\t",col.names=T,row.names=T)
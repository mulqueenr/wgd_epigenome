Make metadata
```python
#meta file generator
import gzip
from Bio import SeqIO
import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pandas as pd

idx="/volumes/seq/projects/wgd/ref/echo_idx.tsv"
meta="/volumes/seq/projects/wgd/231228_RM_WGD_ACT/meta.tsv"

idx=pd.read_table(idx)
met=pd.read_table(meta)

apriori_meta=list()
for index, row in met.iterrows():
        sample=row["sample"]
        i7=row["sample"][:8]
        i5=row["sample"][8:16]
        i7_well=idx.loc[(idx['idx'] == 'i7') & (idx['idx_seq'] == i7)]['well'].values[0]
        i5_well=idx.loc[(idx['idx'] == 'i5') & (idx['idx_seq'] == i5)]['well'].values[0]
        i7_set=idx.loc[(idx['idx'] == 'i7') & (idx['idx_seq'] == i7)]['set'].values[0]
        i5_set=idx.loc[(idx['idx'] == 'i5') & (idx['idx_seq'] == i5)]['set'].values[0]
        if i7_set==1:
                if i5_set==1:
                        plate=1
                else:
                        plate=2
        else:
                if i5_set==1:
                        plate=3
                else:
                        plate=4
        tmp=[sample,plate,i7_well,i5_well]
        apriori_meta.append(tmp)

df = pd.DataFrame(apriori_meta,columns=['sample','plate','row','column'])
df['ploidy_gate'] = ["diploid" if int(row['column']) <= 12 else "aneuploid" for index,row in df.iterrows()]
df['genotype'] = ["WT" if int(row['plate']) <= 2 else "p53ko" for index,row in df.iterrows()]

treatment=list()
for index, row in df.iterrows():
        if row['column'] in ['1','2','13','14']:
                tmp="mon_mpi"
        elif row['column'] in ['3','4','15','16']:
                tmp="dcb"
        elif row['column'] in ['5','6','17','18']:
                tmp="dcb"
        elif row['column'] in ['7','8','19','20']:
                tmp="sp"
        elif row['column'] in ['9','10','21','22']:
                tmp="ro"
        elif row['column'] in ['11','12','23','24']:
                tmp="veh"
        treatment.append(tmp)

df['treatment']=treatment
df.to_csv("/volumes/seq/projects/wgd/231228_RM_WGD_ACT/meta_apriori.csv",index=False)
```

Plot Heatmap

```R
library(copykit)
library(ggplot2)
library(reshape2)
library(HMMcopy)
setwd("/volumes/seq/projects/wgd/231228_RM_WGD_ACT")
meta_apriori<-read.csv("meta_apriori.csv")
meta<-read.csv("meta.tsv",sep="\t")
dat<-readRDS("scCNA.rds")

row.names(meta_apriori)<-meta_apriori$sample
dat$plate<-meta_apriori[row.names(dat@colData),]$plate
dat$plate_row<-meta_apriori[row.names(dat@colData),]$row
dat$plate_column<-meta_apriori[row.names(dat@colData),]$column
dat$ploidy_gate<-meta_apriori[row.names(dat@colData),]$ploidy_gate
dat$treatment<-meta_apriori[row.names(dat@colData),]$treatment
dat$treatment_gate<-paste(dat$treatment,dat$ploidy_gate,sep="_")

# Plot a copy number heatmap with clustering annotation
pdf("subclone.heatmap.expanded.pdf")
plotHeatmap(dat, label = c('treatment_gate','subclones',"plate","ploidy_gate","treatment"),group = 'treatment_gate',order="hclust",row_split="treatment")
dev.off()


dat <- calcInteger(dat, method = 'scquantum', assay = 'smoothed_bincounts')

pdf("/volumes/seq/projects/wgd/231228_RM_WGD_ACT/subclone.ploidy.pdf")
plotMetrics(dat, metric = 'ploidy', label = 'ploidy_score')
dev.off()

saveRDS(dat,"scCNA.rds")

#cell metadata
out_meta<-as.data.frame(dat@colData[c("sample","plate","plate_row","plate_column","ploidy_gate","treatment","reads_assigned_bins")])
colnames(out_meta)<-c("cell_id","plate","plate_row","plate_column","ploidy_gate","treatment","reads_assigned_bins")

#per bin info
range_dat<-data.frame(chr=dat@rowRanges@seqnames,
        start=dat@rowRanges@ranges@start,
        end=dat@rowRanges@ranges@start+dat@rowRanges@ranges@width,
        gc=dat@rowRanges$gc_content)
range_dat$bin_id<-1:nrow(range_dat)

out_dat<-cbind(range_dat,dat@assays@data$bincounts)
out_dat_counts<-melt(dat@assays@data$bincounts)
out_dat_counts$bin_id<-1:nrow(range_dat)
out_dat<-merge(range_dat,out_dat_counts,by="bin_id")

out_dat_state<-melt(dat@assays@data$integer)
out_dat_state$bin_id<-1:nrow(range_dat)
out_dat<-cbind(out_dat,out_dat_state$value)

colnames(out_dat)<-c("bin_id","chr","start","end","gc","cell_id","true_reads_norm","state")
out_dat<-merge(out_dat,out_meta,by="cell_id")
#normalize reads instead of raw bincounts?? seems like they are normalized

cn_s=out_dat[out_dat$ploidy_gate=="aneuploid",]
cn_g1=out_dat[out_dat$ploidy_gate=="diploid",]

write.table(cn_s,file="s_phase_cells_hmmcopy_trimmed.csv",sep=",")
system("gzip s_phase_cells_hmmcopy_trimmed.csv")

write.table(cn_g1,file="g1_phase_cells_hmmcopy_trimmed.csv",sep=",")
system("gzip g1_phase_cells_hmmcopy_trimmed.csv")

```
Output data for PERT
```R
library(copykit)
library(ggplot2)
library(reshape2)
library(HMMcopy)
library(SCOPE)
library(WGSmapp)
library(BSgenome.Hsapiens.UCSC.hg38)

setwd("/volumes/seq/projects/wgd/231228_RM_WGD_ACT")
dat<-readRDS("scCNA.rds")

cell_dir="/volumes/seq/projects/wgd/231228_RM_WGD_ACT/data/cells"
dat_hmm<-runCountReads(cell_dir,
        remove_Y=TRUE,
        resolution="500kb",
        is_paired_end=TRUE)


set_up_ref<-function(bins){ #modified version of SCOPE's get_bam_bed function
  genome <- BSgenome.Hsapiens.UCSC.hg38
  ref <- bins[which(as.character(seqnames(bins)) %in% paste0("chr", c(seq_len(22), "X", "Y")))] #autosomes and X Y

  #Compute mappability for each reference bin.
  mapp_gref<-mapp_hg38 #this is packaged with SCOPE, mappability across bins
  mapp <- rep(1, length(ref))
  #seqlevelsStyle(ref) <- "UCSC"
  for (chr in as.character(unique(seqnames(ref)))) {
      message("Getting mappability for ", chr, sep = "")
      chr.index <- which(as.matrix(seqnames(ref)) == chr)
      ref.chr <- ref[which(as.character(seqnames(ref)) == chr)]
      mapp.chr <- rep(1, length(ref.chr))
      overlap <- as.matrix(findOverlaps(ref.chr, mapp_gref))
      for (i in unique(overlap[, 1])) {
          index.temp <- overlap[which(overlap[, 1] == i), 2]
          overlap.sub <- findOverlaps(ref.chr[i], mapp_gref[index.temp])
          overlap.intersect <- pintersect(ref.chr[i][queryHits(overlap.sub)],mapp_gref[index.temp][subjectHits(overlap.sub)])
          mapp.chr[i] <- sum((mapp_gref$score[index.temp]) * (width(overlap.intersect)))/sum(width(overlap.intersect))
      }
      mapp[chr.index] <- mapp.chr
  }

  #Compute GC for each bin, also from SCOPE
  gc <- rep(NA, length(ref))
  for (chr in unique(seqnames(ref))) {
      message("Getting GC content for chr ", chr, sep = "")
      chr.index <- which(as.matrix(seqnames(ref)) == chr)
      ref.chr <- IRanges(start = start(ref)[chr.index], end = end(ref)[chr.index])
      if (chr == "X" | chr == "x" | chr == "chrX" | chr == "chrx") {
          chrtemp <- "chrX"
      }
      else if (chr == "Y" | chr == "y" | chr == "chrY" | chr == "chry") {
          chrtemp <- "chrY"
      }
      else {
          chrtemp <- as.numeric(mapSeqlevels(as.character(chr), 
              "NCBI")[1])
      }
      if (length(chrtemp) == 0) 
      message("Chromosome cannot be found in NCBI database. ")
      chrm <- unmasked(genome[[chrtemp]])
      seqs <- Views(chrm, ref.chr)
      af <- alphabetFrequency(seqs, baseOnly = TRUE, as.prob = TRUE)
      gc[chr.index] <- round((af[, "G"] + af[, "C"]) * 100, 2)
  }

  ref@elementMetadata$gc<-gc
  ref@elementMetadata$mapp<-mapp
  return(ref)
}

genome <- BSgenome.Hsapiens.UCSC.hg38
ref<-set_up_ref(bins=dat_hmm@rowRanges) #bins is granges of windows to use

#HMM Correction
hmmcopy_sample<-function(x){
  count<-cbind(as.data.frame(ref)[c("seqnames","start","end","gc_content","mapp")],
        dat_hmm@assays@data$bincounts[,x])
  samp<-colnames(dat_hmm@assays@data$bincounts)[x]
  colnames(count)<-c("chr","start","end","gc","map","reads")
  count<-count[c("chr","start","end","reads","gc","map")]
  count<-data.table(count)
  count<-correctReadcount(count)
  count$chr<-as.character(count$chr)
  count<-count[count$chr!="chrY",]
  seg<-HMMsegment(count)
  count$state<-seg$state
  count$sample<-samp
  count$state<-as.character(count$state)
  count$chr<-factor(count$chr,levels=paste0("chr",c(1:22,"X")))
  count<-count[order(count$chr,count$start),]
  count$row_order<-1:nrow(count)
  return(count)
}

hmm_y<-lapply(1:ncol(dat_hmm@assays@data$bincounts),function(x) hmmcopy_sample(x))
out_hmm<-do.call("rbind",hmm_y)

meta<-as.data.frame(dat@colData)[c("reads_assigned_bins","plate","plate_row","plate_column","ploidy_gate","treatment")]
meta$sample<-row.names(meta)

out_hmm<-merge(out_hmm,meta,by="sample")

cn_s=out_hmm[out_hmm$ploidy_gate=="aneuploid",]
cn_g1=out_hmm[out_hmm$ploidy_gate=="diploid",]

write.table(cn_s,file="~/s_phase_cells_hmmcopy_trimmed.csv",sep=",")
system("gzip -f ~/s_phase_cells_hmmcopy_trimmed.csv")

write.table(cn_g1,file="~/g1_phase_cells_hmmcopy_trimmed.csv",sep=",")
system("gzip -f ~/g1_phase_cells_hmmcopy_trimmed.csv")
```

Set up PERT
https://github.com/shahcompbio/scdna_replication_tools
```bash
cd ~/singularity
singularity pull docker://adamcweiner/scdna_replication_tools:main
singularity shell ~/singularity/scdna_replication_tools_main.sif
#also clone repo to scdna_replication_tools-main.zip
```

https://github.com/shahcompbio/scdna_replication_tools/blob/main/notebooks/inference_tutorial_pt2.ipynb
```python
singularity shell ~/singularity/scdna_replication_tools_main.sif

import sys
sys.path.append('/volumes/USR2/Ryan/tools/scdna_replication_tools') #add this path as a module
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scdna_replication_tools.infer_scRT import scRT
from scdna_replication_tools.plot_pert_output import plot_model_results
from scdna_replication_tools.compute_ccc_features import compute_ccc_features
from scdna_replication_tools.predict_cycle_phase import predict_cycle_phase

from scgenome import refgenome



# load the simulated data from the paper
# this corresponds to the data seen in simulated dataset P8.2 which has subclonal and cell-specific CNAs
cn_s = pd.read_csv('/volumes/USR2/Ryan/tools/scdna_replication_tools/data/P8.2/s_phase_cells_hmmcopy_trimmed.csv.gz', dtype={'chr': str})
cn_g1 = pd.read_csv('/volumes/USR2/Ryan/tools/scdna_replication_tools/data/P8.2/g1_phase_cells_hmmcopy_trimmed.csv.gz', dtype={'chr': str})

#load data
cn_s = pd.read_csv('s_phase_cells_hmmcopy_trimmed.csv.gz', dtype={'chr': str})
cn_g1 = pd.read_csv('g1_phase_cells_hmmcopy_trimmed.csv.gz', dtype={'chr': str})

# add the replication columns for the G1-phase cells
cn_g1['true_rep'] = 0.0
cn_g1['true_p_rep'] = 0.0
cn_g1['true_t'] = 1.0
cn_g1['true_G1_state']=1.0

cn_s['true_rep'] = 0.0
cn_s['true_p_rep'] = 0.0
cn_s['true_t'] = 1.0
cn_s['true_G1_state']=0.0

cn_g1['true_reads_norm']=(cn_g1['reads']/cn_g1['reads_assigned_bins'])*1000000
cn_g1['clone_id']=cn_g1['treatment']
cn_g1['cell_id']=cn_g1['sample']
cn_g1['chr']=[i[3:] for i in cn_g1['chr']]
cn_s['true_reads_norm']=(cn_s['reads']/cn_s['reads_assigned_bins'])*1000000
cn_s['clone_id']=cn_s['treatment']
cn_s['cell_id']=cn_s['sample']
cn_s['chr']=[i[3:] for i in cn_s['chr']]

# plot the true simulated data
plot_model_results(
    cn_s, cn_g1, 
    clone_col='clone_id', second_sort_col='true_t', 
    input_cn_col='state', output_cn_col='true_G1_state', 
    output_rep_col='true_rep', rpm_col='true_reads_norm',
    rpm_title='Read depth', input_cn_title='HMMcopy states',
    output_cn_title='True somatic CN states', rep_title='True replication states'
)
plt.show()
plt.savefig("pert_g2m_s.png")

# note the true cell cycle phase of each cell and merge all the cells into one dataframe
cn_s['true_phase'] = 'S'
cn_g1['true_phase'] = 'G1/2'

cn = pd.concat([cn_s, cn_g1], ignore_index=True)
cn.head()
# compute per-cell features for all cells
cn, cell_features = compute_ccc_features(
    cn, 
    rpm_col='true_reads_norm', 
    cn_col='state',
    clone_col='treatment', 
    madn_col='madn', 
    lrs_col='lrs',
    num_reads_col='reads_assigned_bins' 
)
cell_features.head()
cn.head()

# merge the true_phase information back into the cell_features dataframe
cell_features = cell_features.merge(cn[['cell_id', 'true_phase']].drop_duplicates(), on='cell_id')
cell_features.head()

# split the data so cells with the 400 highest corrected breakpoints are in the G1/2-phase dataset
# sort the cells by corrected breakpoints
cell_features = cell_features.sort_values('corrected_breakpoints', ascending=False)
# get the cell_ids for the top 400 cells and the threshold for the 401st cell
top_400_cell_ids = list(cell_features['cell_id'].iloc[:405].values)
thresh = cell_features['corrected_breakpoints'].iloc[405]

# look at distribution of corrected breakpoints to set a threshold for identifying the high-confidence G1/2-phase cells
sns.histplot(cell_features['corrected_breakpoints'], kde=True, bins=100)
# plot a vertical line at the threshold
plt.axvline(thresh, color='red', linestyle='--')
plt.show()
plt.savefig("pert_g2m_s2.png")

# # split the data into G1/2-phase and S-phase
cn_s_input = cn[cn['true_G1_state'] == 0.0].reset_index(drop=True)
cn_g_input = cn[cn['true_G1_state'] == 1.0].reset_index(drop=True)
cn_g_input.cell_id.unique().shape

cn_s_input['library_id']=cn_s_input['treatment']+"_"+cn_s_input['ploidy_gate']
cn_g_input['library_id']=cn_g_input['treatment']+"_"+cn_g_input['ploidy_gate']

# temporarily remove columns that don't get used by PERT
temp_cn_s = cn_s_input[['cell_id', 'chr', 'start', 'end', 'gc', 'state','copy', 'library_id', 'true_reads_norm']]
temp_cn_g1 = cn_g_input[['cell_id', 'chr', 'start', 'end', 'gc', 'state', 'copy', 'library_id', 'true_reads_norm']]

print('number of cells in S-phase: ', len(temp_cn_s['cell_id'].unique()))
print('number of cells in G1-phase: ', len(temp_cn_g1['cell_id'].unique()))
print('number of loci in S-phase: ', len(temp_cn_s[['chr', 'start']].drop_duplicates()))
print('number of loci in G1-phase: ', len(temp_cn_g1[['chr', 'start']].drop_duplicates()))


# create scRT object with input columns denoted
scrt = scRT(temp_cn_s, temp_cn_g1, 
        input_col='true_reads_norm', 
        clone_col='library_id', 
        assign_col='copy', 
        rt_prior_col=None,
        cn_state_col='state', 
        gc_col='gc', 
        cn_prior_method='g1_composite', 
        max_iter=300)

# run inference using PERT
cn_s_with_scrt, supp_s_output, cn_g_with_scrt, supp_g_output = scrt.infer(level='pyro')
cn_s_with_scrt.head()

cn_g_with_scrt.head()

# plot the loss curve for 'loss_s' and 'loss_g' param values in the supp_s_output and supp_g_output dataframes
fig, ax = plt.subplots(1, 3, figsize=(12, 4))
sns.lineplot(data=supp_s_output.query("param=='loss_g'"), x='level', y='value', ax=ax[0])
sns.lineplot(data=supp_s_output.query("param=='loss_s'"), x='level', y='value', ax=ax[1])
sns.lineplot(data=supp_g_output.query("param=='loss_s'"), x='level', y='value', ax=ax[2])

for i in range(3):
    ax[i].set_xlabel('Iteration')
    ax[i].set_ylabel('Loss')
    ax[i].set_title('Loss curve PERT Step {}'.format(i+1))

plt.show()
plt.savefig("tutorial2_true_sim_losscurve.png")

cn_s_with_scrt['clone_id']=cn_s['treatment']
cn_g_with_scrt['clone_id']=cn_g1['treatment']

# plot the results of the inference
plot_model_results(
    cn_s_with_scrt, cn_g_with_scrt, clone_col='library_id', second_sort_col='model_tau',
    input_cn_col='state', output_cn_col='model_cn_state', output_rep_col='model_rep_state',
    rpm_col='true_reads_norm', rpm_title='Read depth', input_cn_title='HMMcopy states',
    output_cn_title='PERT CN states', rep_title='PERT replication states',
    top_title_prefix='Input set: unknown', bottom_title_prefix='Input set: high-confidence G1/2-phase'
)
plt.show()
plt.savefig("tutorial2_true_sim_losscurve.png")

# concatenate the two dataframes containing PERT output
cn_out = pd.concat([cn_s_with_scrt, cn_g_with_scrt], ignore_index=True)

# predict the cycle phase for each clone based on the PERT output with default parameters
cn_s_out, cn_g_out, cn_lq_out = predict_cycle_phase(
    cn_out, frac_rt_col='cell_frac_rep', rep_state_col='model_rep_state', 
    cn_state_col='model_cn_state', rpm_col='true_reads_norm'
)

cn_s_out.cell_id.unique().shape, cn_g_out.cell_id.unique().shape, cn_lq_out.cell_id.unique().shape

# plot the results of the inference, now sorted by the PERT predicted phase
plot_model_results(
    cn_s_out, cn_g_out, clone_col='library_id', second_sort_col='model_tau',
    input_cn_col='state', output_cn_col='model_cn_state', output_rep_col='model_rep_state',
    rpm_col='true_reads_norm', rpm_title='Read depth', input_cn_title='HMMcopy states',
    output_cn_title='PERT xCN states', rep_title='PERT replication states',
    top_title_prefix='PERT: S-phase', bottom_title_prefix='PERT: G1/2-phase'
)
plt.show()
plt.savefig("tutorial2_true_sim_inferreddat.png")

# show what the PERT predicted low-quality cells look like in contrast to the predicted S-phase cells
plot_model_results(
    cn_s_out, cn_lq_out, clone_col='library_id', second_sort_col='model_tau',
    input_cn_col='state', output_cn_col='model_cn_state', output_rep_col='model_rep_state',
    rpm_col='true_reads_norm', rpm_title='Read depth', input_cn_title='HMMcopy states',
    output_cn_title='PERT CN states', rep_title='PERT replication states',
    top_title_prefix='PERT: S-phase', bottom_title_prefix='PERT: LQ'
)
plt.show()
plt.savefig("tutorial2_true_sim_inferreddat_lowquality.png")


# merge the true time in S-phase with the PERT cell_frac_rep output for all cells
output_cell_times = pd.concat([
    cn_s_out[['cell_id', 'cell_frac_rep']].drop_duplicates(), 
    cn_g_out[['cell_id', 'cell_frac_rep']].drop_duplicates(), 
    cn_lq_out[['cell_id', 'cell_frac_rep']].drop_duplicates()
    ], ignore_index=True
)

true_cell_times = pd.concat([
    cn_s[['cell_id', 'true_t', 'true_phase']].drop_duplicates(), 
    cn_g1[['cell_id', 'true_t', 'true_phase']].drop_duplicates(), 
    ], ignore_index=True
)

cell_times = pd.merge(true_cell_times, output_cell_times, on='cell_id')

print(cell_times.shape)

cell_times.head()

sns.jointplot(data=cell_times, x='true_t', y='cell_frac_rep', hue='true_phase', alpha=0.1)
plt.xlabel('True time in S-phase')
plt.ylabel('PERT predicted time in S-phase')
plt.show()
plt.savefig("tutorial2_true_sim_inferreddat_dotplot.png")




```
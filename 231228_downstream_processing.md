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

```R
library(copykit)
library(ggplot2)

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
plotHeatmap(dat, label = c('treatment_gate','subclones',"plate","ploidy_gate","treatment"),group = 'treatment_gate',order="hclust")
dev.off()
```
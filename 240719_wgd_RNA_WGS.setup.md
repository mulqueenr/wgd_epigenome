Transfer data from geo onto seadragon

```bash
bsub -Is -W 6:00 -q interactive -n 1 -M 16 -R rusage[mem=16] /bin/bash #get interactive transfer node this has internet access for environment set up
```

Use SFTP to get run data
```bash
proj_dir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/wgd"
mkdir -p ${proj_dir}/rundir/
cd ${proj_dir}/rundir/
sftp mulqueen@qcprpgeo.mdanderson.edu
get -r /volumes/seq/flowcells/MDA/nextseq2000/2024/20240720_spDNA_BCIS41T_CHip1_TArgetDR_test2_55_61degree_Ryan_ArcDR_2nd_ECIS61T_BCIS120T_cDNA #grab full rundir
```

Set up directory for run analysis
```bash
mkdir -p ${proj_dir}/240515_RNADNA_WGD_t18
cp ../240505_RNADNA_WGD_t0/sample_layout.txt . #same sample layout as before
cp ../240505_RNADNA_WGD_t0/RNA_WGS.nextflow.submit.lsf . 
#modify RNA_WGS.nextflow.submit.lsf with new variables
```

Run analysis through seadragon

```bash
#BSUB -W 24:00
#BSUB -q e80medium
#BSUB -n 70
#BSUB -M 800
#BSUB -R rusage[mem=800]
#BSUB -o /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/wgd/240515_RNADNA_WGD_t18/%J_processing.log
#BSUB -cwd /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/wgd/240515_RNADNA_WGD_t18
#BSUB -u RMulqueen@mdanderson.org
#BSUB -J WGD_RNA_WGS

#load modules
module load nextflow/23.04.3

#set up environment variables 
export SCRATCH="/rsrch4/scratch/genetics/rmulqueen"
export projDir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/wgd"
export srcDir="${projDir}/src"
export refDir="${projDir}/ref"
export runDir="${projDir}/rundir/20240720_spDNA_BCIS41T_CHip1_TArgetDR_test2_55_61degree_Ryan_ArcDR_2nd_ECIS61T_BCIS120T_cDNA/" #change this
export icell8_data="${projDir}/240515_RNADNA_WGD_t18/240515_WDR_Fixed_WGDt18d/2024.05.15.14.59-141436/141436_secondscan_WellList.TXT"

export outDir="${projDir}/240515_RNADNA_WGD_t18"
export sample_layout="${outDir}/sample_layout.txt"

mkdir -p ${SCRATCH}/wgd_work 

#call nextflow
nextflow ${srcDir}/WGD_RNA_WGS.groovy \
--projectdir $projDir \
--sequencing_dir $runDir \
--icell8_data $icell8_data \
--out_dir $outDir \
--sample_layout $sample_layout \
--dna_chip_primer MD_N705_out \
--rna_chip_primer DDR_PCR_p5_UDI4_v2 \
-w ${SCRATCH}/wgd_work \
-resume
```

```bash
bsub < 240515_WGD_RNAWGS.nextflow.submit.lsf
```
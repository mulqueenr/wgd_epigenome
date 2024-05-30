nextflow.enable.dsl=2
// Reference files
params.projectdir="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/wgd"
params.src_dir="${params.projectdir}/src/"
params.ref_dir="${params.projectdir}/ref/"

// Input parameters, user specified defaults
params.sequencing_dir = "${params.projectdir}/rundir/240425_VH00219_582_AACFVYNHV"
params.icell8_data="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/wgd/240505_RNADNA_WGD_t0/240322_ARCDR_FIXED_WGDt0_141014C_scan2/141014-manual3_reprocess_WellList.TXT"
params.out_dir="${params.projectdir}/240505_RNADNA_WGD_t0"
params.sample_layout="${params.out_dir}/sample_layout.txt"
params.outname="WGD"
params.dna_chip_primer="MD_N708_out"
params.rna_chip_primer="DDR_PCR_p5_UDI8_V2"
params.readcountfilter="100000"

log.info """

		================================================================
		    WGD PROJECT RNA/WGS ANALYSIS v0.1
		    Based on Kaile Wang and Rui Ye's work and original pipeline
		================================================================
		==Input==
		Input Sequencing Dir : ${params.sequencing_dir}
		Input icell8 Data : ${params.icell8_data}
		Input icell8 Sample Names : ${params.sample_layout}
		Input icell8 DNA Chip Primer : ${params.dna_chip_primer}
		Input icell8 RNA Chip Primer : ${params.rna_chip_primer}

		==Params==
		NF Working Dir : ${workflow.workDir}
		Src directory : ${params.src_dir}
		Read Count Filter (DNA): ${params.readcountfilter}

		==Output==
		Output Project Directory : ${params.projectdir}
		================================================

""".stripIndent()

//PREPROCESSING BLOCK
process BCL_TO_FASTQ  { 
	//Generate Undetermined Fastq Files for Processing
	containerOptions "--bind ${params.ref_dir}:/ref/"
	cpus 50
	label 'bcl2fastq'

	input:
		path(runDir)
	output:
		path("*.fastq.gz")
	script:
	"""
	bcl2fastq \
	--use-bases-mask=Y46,I34,I8,Y50 \
	--create-fastq-for-index-reads \
	-r 20 -p 10 -w 20 \
	--no-lane-splitting \
	-R ${runDir} \
	-o .
	"""
}


process DEMUX_FASTQ { 
	//Assign fastq files to determined index pairs.
	containerOptions "--bind ${params.ref_dir}:/ref/,${params.src_dir}:/src/"
	publishDir "${params.out_dir}", pattern:'metadata.csv', mode: 'copy', overwrite: true
	cpus 50
	label 'bcl2fastq'

	input:
		tuple path(read1), path(read2), path(idx1), path(idx2)
		path(sample)
		path(sample_layout)
	output:
		tuple path("${params.outname}.R1.DNA.fastq.gz"), path("${params.outname}.R2.DNA.fastq.gz")
		tuple path("${params.outname}_S1_L001_R1_001.fastq.gz"), path("${params.outname}_S1_L001_R2_001.fastq.gz")

	script:
		"""
		#chunk fastq to speed up processing
		seqkit split2 ${read1} -p 50 -j ${task.cpus} -O . --by-part-prefix ${read1.simpleName}.chunk -e .gz 
		seqkit split2 ${read2} -p 50 -j ${task.cpus} -O . --by-part-prefix ${read2.simpleName}.chunk -e .gz 
		seqkit split2 ${idx1} -p 50 -j ${task.cpus} -O . --by-part-prefix ${idx1.simpleName}.chunk -e .gz 
		seqkit split2 ${idx2} -p 50 -j ${task.cpus} -O . --by-part-prefix ${idx2.simpleName}.chunk -e .gz 

		#run demultiplexing
		demux() { 
		python /src/fastqsplitter.wdrfixed.py \\
		--fq1 Undetermined_S0_R1_001.chunk0\${1}.fastq.gz \\
		--fq2 Undetermined_S0_R2_001.chunk0\${1}.fastq.gz  \\
		--idx3 Undetermined_S0_I1_001.chunk0\${1}.fastq.gz  \\
		--idx4 Undetermined_S0_I2_001.chunk0\${1}.fastq.gz \\
		--samples ${sample} \\
		--sample_layout ${sample_layout} \\
		--dna_chip_primer ${params.dna_chip_primer} \\
		--rna_chip_primer ${params.rna_chip_primer}

		gzip *chunk0\${1}*fastq
		}

		export -f demux
		parallel -j ${task.cpus} demux ::: \$(echo {01..50})

		zcat *R1*DNA.barc.fastq.gz > ${params.outname}.R1.DNA.fastq.gz
		zcat *R2*DNA.barc.fastq.gz > ${params.outname}.R2.DNA.fastq.gz
		zcat *R1*RNA.barc.fastq.gz > ${params.outname}_S1_L001_R1_001.fastq.gz #have to follow naming convention for cellranger
		zcat *R2*RNA.barc.fastq.gz > ${params.outname}_S1_L001_R2_001.fastq.gz 
		"""
}

// DNA ALIGNMENT AND SPLITTING CELLS 
process DNA_BWA_ALIGN {
	//Map reads with BWA mem
	containerOptions "--bind ${params.ref_dir}:/ref/,${params.src_dir}:/src/"
	cpus 50
	label 'bcl2fastq'

	input:
		tuple path(dna_fq1), path(dna_fq2)
		path bwa_index
		val outname
	output:
		tuple val("${outname}"), path("${outname}.DNA.bam")
	script:
		def idxbase = bwa_index[0].baseName
		"""
		bwa mem \\
		-t ${task.cpus} \\
		${idxbase} \\
		${dna_fq1} \\
		${dna_fq2} \\
		| samtools view -@ ${task.cpus} -b - > ${outname}.DNA.bam
		"""
}

process DNA_SPLIT_BAM { 
	// Generate a count per grouped bam and pass list.
	//Simplified version of the scalebio met postprocessing function (accounts for only one bam)
	cpus 50
	label 'cnv'

	input:
		tuple val(outname), path(bam)
	output:
		path("*.sorted.bam"), optional: true
	script:
	"""
	samtools view -@ ${task.cpus} ${bam} \\
	| awk -v b=${bam} '{split(\$1,a,":"); print a[1],b}' \\
	| sort \\
	| uniq -c \\
	| sort -k1,1n \\
	| awk '\$1>${params.readcountfilter} {print \$0}' > readcount.tsv


	split_bam() {
	test=\$1
	idx=\$(echo \$test | cut -d \' \' -f 2 )
	bam=\$(echo \$test | cut -d \' \' -f 3)

	((samtools view -H \$bam) && (samtools view \$bam | awk -v i=\$idx \'{split(\$1,a,\":\"); if(a[1]==i) print \$0}\')) \\
	| samtools view -bS - \\
	| samtools sort -T . -O BAM -o \${idx}.sorted.bam - 
	}
    
    export -f split_bam
	parallel -j ${task.cpus} -a readcount.tsv split_bam 
    """
}

process DNA_BAM_MARKDUP {
	//Fix mates, sort and mark duplicates in bam
	publishDir "${params.out_dir}/dna_data/cells", mode: 'copy'
	label 'cnv'

	input:
		path bam
	output:
		path("*bbrd.bam")
	script:
	"""
	samtools sort -n -o - ${bam} \\
	| samtools fixmate -m - - \\
	| samtools sort -T . -o - - \\
	| samtools markdup -s - ${bam.baseName}.bbrd.bam 2> ${bam.baseName}.rmdup.stats.txt
	"""
}


// CNV PROFILING 
process DNA_CNV_CLONES {
	//Run CopyKit and output list of bam files by clones
	containerOptions "--bind ${params.src_dir}:/src/"
	publishDir "${params.out_dir}/dna_data", mode: 'copy', pattern: "*"
	cpus 50
	label 'cnv'

	input:
		path sc_sorted_bam
	output:
		path("*rds")
	script:
		"""
		Rscript /src/copykit_run.R \\
		-i "." \\
		-c ${task.cpus}
		"""
}

// RNA PROCESSING
process RNA_CELLRANGER {
	//Run Cellranger on flex rna output
	containerOptions "--bind ${params.ref_dir}:/ref/,${params.src_dir}:/src/"
	publishDir "${params.out_dir}/rna_data", mode: 'copy', pattern: "*"
	cpus 50
	label 'cellranger'

	input:
		tuple path(rna_fq1),path(rna_fq2)
	output:
		path("*")
	script:
	"""
	echo \"""[gene-expression]
	reference,/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
	probe-set,/ref/Chromium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv
	create-bam,true
	r1-length,28
	chemistry,SFRP

	[libraries]
	fastq_id,fastqs,feature_types
	WGD,\$PWD,Gene Expression\""" > sample.csv

	cellranger multi \\
	--id=WGD \\
	--csv=./sample.csv \\
	--localcores ${task.cpus}
	"""
}

workflow {
	// SETTING UP VARIABLES
		bwa_index = file("${params.ref_dir}/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa{,.amb,.ann,.bwt,.pac,.sa}" )
		def fasta_ref = Channel.value(params.ref_dir)
		def outname = Channel.value(params.outname)
		def icell8_data= Channel.fromPath(params.icell8_data)
		def sample_layout= Channel.fromPath(params.sample_layout)

	// BCL TO FASTQ PIPELINE 
		fq = Channel.fromPath(params.sequencing_dir) | BCL_TO_FASTQ
		(dna_fq, rna_fq) = DEMUX_FASTQ(fq,icell8_data,sample_layout)

	// DNA ALIGNMENT AND SPLITTING CELLS AND CNV CALLING
		clone_lists  = DNA_BWA_ALIGN(dna_fq, bwa_index, outname) \
		| DNA_SPLIT_BAM \
		| flatten
		| DNA_BAM_MARKDUP \
		| flatten \
		| collect \
		| DNA_CNV_CLONES

	//RNA PROCESSING VIA CELLRANGER
		RNA_CELLRANGER(rna_fq)
	
}

/*
bsub -Is -W 36:00 -q long -n 10 -M 100 -R rusage[mem=100] /bin/bash

sif="/rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/wgd/src/copykit.sif"

module load singularity

singularity shell \
--bind /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/wgd/ \
--bind /rsrch4/scratch/genetics/rmulqueen \
--bind /rsrch4/scratch/genetics/rmulqueen/wgd_work \
--bind /rsrch5/home/genetics/NAVIN_LAB/Ryan/projects/wgd/src:/src/ \
$sif 

cd /rsrch4/scratch/genetics/rmulqueen/wgd_work/a6/7e935b


*/


/*
DEPROCATED FUNCTIONS
process DNA_SPLIT_BAM_BY_READNAME {
	//Split bam file by read names
	label 'cnv'

	input:
		tuple val(outname), path(bam)
	output:
		path("*sam")

	script:
	"""
	samtools view ${bam} \\
	| awk 'OFS="\\t" {split(\$1,a,":"); print \$0,"XM:Z:"a[1] > a[1]".${outname}.DNA.sam"}'
	"""
}
// CONVERT TO BAMS AND RUN QC 
process DNA_SAM_TO_BAM_CONVERT {
	//Convert sam to bam and headers on single cell bams
	//Hard filters to cells with greater than at least 100k reads.
	label 'cnv'
	containerOptions "--bind ${params.ref_dir}:/ref/"

	input:
		path(sc_sams)
		val fasta_ref
	output:
		path("*.sc.bam"), optional: true

	script:
	"""
	if [[ \$(wc -l < ${sc_sams}) -ge 100000 ]]; then
	samtools view -bT \\
	/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa \\
	${sc_sams} \\
	| samtools sort -o ${sc_sams.baseName}.sc.bam - 
	fi
	"""
}
*/




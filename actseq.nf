//Nextflow pipeline for processing Navin lab standard ACT seq libraries//

// Declare syntax version
//https://github.com/danrlu/Nextflow_cheatsheet/blob/main/nextflow_cheatsheet.pdf
//Runs entirely in conda environment singlecell_hic_pipeline_env
nextflow.enable.dsl=2

// Script parameters
params.flowcellDir = "/volumes/seq/flowcells/MDA/nextseq2000/2023/231228_RM_WGD_ACT/231229_VH00219_560_AAFF2CWM5"
params.outname = "231228_RM_WGD_ACT"
params.outdir = "/volumes/seq/projects/wgd/231228_RM_WGD_ACT"
params.ref = "/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"
params.plate="1,2,3,4"
params.platform="echo"

log.info """

		================================================
		             WGD-ACTseq NF PIPELINE v1.0
		================================================
		Flowcell Dir : ${params.flowcellDir}
		Output Prefix : ${params.outname}
		NF Working Dir : ${workflow.launchDir}
		Output Directory : ${params.outdir}
		Experiment Platform : ${params.platform}
		Experiment Plate : ${params.plate}
		================================================

""".stripIndent()

// BCL TO FASTQ PIPELINE 
process BCL_TO_FASTQ { 
	//Generate Undetermined Fastq Files from BCL Files.
	cpus 100

	input:
		path flowcellDir
	output:
		tuple path("Undetermined_S0_L001_R1_001.fastq.gz"),
		path("Undetermined_S0_L001_R2_001.fastq.gz"),
		path("Undetermined_S0_L001_I1_001.fastq.gz"),
		path("Undetermined_S0_L001_I2_001.fastq.gz")
	script:
		"""
		bcl2fastq -R $flowcellDir \\
		-o . \\
		-r 10 \\
		-p ${task.cpus} \\
		-w 10 \\
		--ignore-missing-bcls \\
		--ignore-missing-filter \\
		--ignore-missing-positions \\
		--ignore-missing-controls \\
		--create-fastq-for-index-reads
		"""
}

process CHUNK_FASTQ {
	//Chunk fastq files to 50 chunks reads to parallelize demultiplexing
	cpus 150

	input:
		tuple path(read1), path(read2), path(ind1), path(ind2)
	output:
		tuple path("${read1.simpleName}.chunk*gz"), path("${read2.simpleName}.chunk*gz"),
			path("${ind1.simpleName}.chunk*gz"), path("${ind2.simpleName}.chunk*gz")
	script:
	"""
	seqkit split2 ${read1} -p 50 -j ${task.cpus} -O . --by-part-prefix ${read1.simpleName}.chunk -e .gz &
	seqkit split2 ${read2} -p 50 -j ${task.cpus} -O . --by-part-prefix ${read2.simpleName}.chunk -e .gz &
	seqkit split2 ${ind1} -p 50 -j ${task.cpus} -O . --by-part-prefix ${ind1.simpleName}.chunk -e .gz &
	seqkit split2 ${ind2} -p 50 -j ${task.cpus} -O . --by-part-prefix ${ind2.simpleName}.chunk -e .gz &
	"""
}

process DEMUX_FASTQ { 
	//Assign fastq files to determined index pairs.
	cpus 50

	input:
		tuple path(read1), path(read2), path(idx1), path(idx2)
	output:
		path("*R1*.barc.fastq")
		path("*R2*.barc.fastq")
	script:
		"""
		python /volumes/seq/projects/wgd/src/fastqsplitter.nf.py \\
		--fq1 $read1 \\
		--fq2 $read2 \\
		--idx3 $idx1 \\
		--idx4 $idx2 \\
		--platform ${params.platform} \\
		--plate ${params.plate}
		"""
}

process UNCHUNK_FASTQ {
	//Chunk fastq files to 5M reads to parallelize demultiplexing
	
	input:
		path(read1)
		path(read2)
	output:
		tuple path("test.R1.fq.gz"), path("test.R2.fq.gz")
	script:
	"""
	cat ${read1} | gzip > test.R1.fq.gz &
	cat ${read2} | gzip > test.R2.fq.gz &
	"""
}

// ALIGNMENT AND SPLITTING CELLS 

process BWA_ALIGN {
	//Map reads with BWA mem
	cpus 50

	input:
		tuple path(read1), path(read2)
		path bwa_index
		val outname
	output:
		tuple val("${outname}"), path("${outname}.bam")
	script:
		def idxbase = bwa_index[0].baseName
		"""
		bwa mem \\
		-t ${task.cpus} \\
		${idxbase} \\
		${read1} \\
		${read2} \\
		| samtools view -@ ${task.cpus} -b - > ${outname}.bam
		"""
}

process SPLIT_BAM_BY_READNAME {
	//Split bam file by read names

	input:
		tuple val(outname), path(bam)
	output:
		path("*sam")

	script:
	"""
	samtools view ${bam} \\
	| awk 'OFS="\\t" {split(\$1,a,":"); print \$0,"XM:Z:"a[1] > a[1]".${bam.simpleName}.sam"}'
	"""
}

// CONVERT TO BAMS AND RUN QC 
process SAM_TO_BAM_CONVERT {
	//Convert sam to bam and headers on single cell bams
	//Hard filters to cells with greater than at least 100k reads.

	input:
		path sc_sams
		val fasta_ref
	output:
		path("*.sc.bam"), optional: true

	script:
	"""
	if [[ \$(wc -l < ${sc_sams}) -ge 100000 ]]; then
	samtools view -bT ${fasta_ref} ${sc_sams} | samtools sort -o ${sc_sams.baseName}.sc.bam - 
	fi
	"""
}

process BAM_MARKDUP {
	//Fix mates, sort and mark duplicates in bam
	publishDir "${params.outdir}/data/cells", mode: 'copy'

	input:
		path bam
	output:
		path("*.sc.bbrd.bam")
	script:
	"""
	samtools sort -n -o - ${bam} \\
	| samtools fixmate -m - - \\
	| samtools sort -T . -o - - \\
	| samtools markdup -s - ${bam.baseName}.bbrd.bam 2> ${bam.baseName}.rmdup.stats.txt
	"""
}

process FASTQC {
	//Generate FastQC per bam
	//Based on https://github.com/nf-core/rnaseq/blob/37f260d360e59df7166cfd60e2b3c9a3999adf75/main.nf#L473
	publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

	input:
		path bam
	output:
    	path("*_fastqc.{zip,html}")
	script:
	"""
	fastqc ${bam}
	"""
}

process MULTIQC {
	//Run MultiQC
	publishDir "${params.outdir}/multiqc", mode: 'copy'

	input:
		path fastqc
	output:
		path("*")
	script:
	"""
	multiqc --outdir . -f .
	"""
}

// CNV PROFILING 
process CNV_CLONES {
	//Run CopyKit and output list of bam files by clones
	publishDir "${params.outdir}/data", mode: 'copy', pattern: "*"
	cpus 50

	input:
		path sc_dedup_bams
	output:
		path("*bam_list*txt")
	script:
		"""
		Rscript /volumes/seq/projects/wgd/src/copykit_cnv_clones.nf.R \\
		. \\
		${params.outdir} \\
		${task.cpus}
		"""
}


workflow {
	// SETTING UP VARIABLES
	bwa_index = file("/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa{,.amb,.ann,.bwt,.pac,.sa}" )
	def fasta_ref = Channel.value(params.ref)
	def outname = Channel.value(params.outname)

	// BCL TO FASTQ PIPELINE 
		(fq1, fq2) = Channel.fromPath(params.flowcellDir) \
		| BCL_TO_FASTQ \
		| CHUNK_FASTQ \
		| transpose \
		| DEMUX_FASTQ

		fq1 = fq1 | collect
		fq2 = fq2 | collect

		unchunked_fqs = UNCHUNK_FASTQ(fq1, fq2)

	// ALIGNMENT AND SPLITTING CELLS 
		sc_sams = BWA_ALIGN(unchunked_fqs, bwa_index, outname) \
		| SPLIT_BAM_BY_READNAME \
		| flatten

		sc_dedup_bams = SAM_TO_BAM_CONVERT(sc_sams, fasta_ref) \
		| flatten \
		| BAM_MARKDUP \
		| flatten

		bam_in = sc_dedup_bams \
		| collect

	// QC ON CELLS
		FASTQC(sc_dedup_bams) \
		| collect \
		| MULTIQC
		
	// CALL CNVs AND GENERATE CLONE LISTS
		clone_lists = CNV_CLONES(bam_in) \
		| flatten
}

/*
default run:
nextflow actseq.nf \
-resume \
-with-report \
--flowcellDir /volumes/seq/flowcells/MDA/nextseq2000/2023/231228_RM_WGD_ACT/231229_VH00219_560_AAFF2CWM5 \
--outname 231228_RM_WGD_ACT \
--plate 1,2,3,4 \
--platform echo \
--outdir /volumes/seq/projects/wgd/231228_RM_WGD_ACT
*/



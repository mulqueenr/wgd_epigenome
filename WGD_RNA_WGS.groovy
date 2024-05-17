nextflow.enable.dsl=2

// Make samplesheet per cell, run bcl2fastq, process RNA via cellranger, process DNA via copykit

//#RunParameters.xml #actual sequencing input, I'm guessing
//    <Read1>46</Read1>
//    <Index1>34</Index1>
//    <Index2>8</Index2>
//    <Read2>50</Read2>

//#RunInfo_backup.xml #RNA?
//<Reads>
//<Read Number="1" NumCycles="46" IsIndexedRead="N" IsReverseComplement="N"/>
//<Read Number="2" NumCycles="34" IsIndexedRead="Y" IsReverseComplement="N"/> #pool barcode+icell8columnidx
//<Read Number="3" NumCycles="8" IsIndexedRead="Y" IsReverseComplement="Y"/> #icell8rowidx
//<Read Number="4" NumCycles="50" IsIndexedRead="N" IsReverseComplement="N"/>

//#RunInfo.xml #DNA?
//<Read Number="1" NumCycles="46" IsIndexedRead="N" IsReverseComplement="N"/>
//<Read Number="2" NumCycles="26" IsIndexedRead="Y" IsReverseComplement="N"/>
//<Read Number="3" NumCycles="16" IsIndexedRead="Y" IsReverseComplement="Y"/>
//<Read Number="4" NumCycles="50" IsIndexedRead="N" IsReverseComplement="N"/>

// Reference files
params.projectdir="/volumes/seq/projects/wgd/"
params.src_dir="${params.projectdir}/src/"
params.wafar_cell_bc = "${params.src_dir}/Wafar_cell.csv"
params.wdr_fix_cell_bc = "${params.src_dir}/WDR_FIX_cellBC.txt"

// Input parameters, user specified defaults
params.sequencing_dir = "/volumes/seq/flowcells/MDA/nextseq2000/2024/20240426_Ryan_231_wgdT0_DNA_RNA_ESO_P17_Chip1_DNA_RNA_reseq_P18_Chip2_DNA/240425_VH00219_582_AACFVYNHV"
params.icell8_filter="/volumes/lab/users/wet_lab/records/Navin_lab_projects/WGD/icell8/RNAWGS/240322_ARCDR_FIXED_WGDt0_141014C_scan2/141014-manual3_FilterFile.csv"
params.icell8_data="/volumes/lab/users/wet_lab/records/Navin_lab_projects/WGD/icell8/RNAWGS/240322_ARCDR_FIXED_WGDt0_141014C_scan2/141014-2_WellList.TXT"
params.out_dir="${params.projectdir}/240322_ARCDR_FIXED_WGDt0"
params.sample_layout="${params.out_dir}/sample_layout.csv"
params.bcl_dir="/volumes/seq/flowcells/MDA/nextseq2000/2024/20240426_Ryan_231_wgdT0_DNA_RNA_ESO_P17_Chip1_DNA_RNA_reseq_P18_Chip2_DNA/240425_VH00219_582_AACFVYNHV"
log.info """

		================================================================
		    WGD PROJECT RNA/WGS ANALYSIS v0.1
		    Based on Kaile Wang and Rui Ye's work and original pipeline
		================================================================
		==Input==
		Input Sequencing Dir : ${params.sequencing_dir}
		Input icell8 Filter : ${params.icell8_filter}
		Input icell8 Data : ${params.icell8_filter}
		Input icell8 Sample Names : ${params.icell8_samples}

		==Params==
		NF Working Dir : ${workflow.workDir}
		Src directory : ${params.src_dir}

		==Output==
		Output Project Directory : ${params.projectdir}
		================================================

""".stripIndent()


//PREPROCESSING BLOCK
process MAKE_SAMPLESHEET { 
	publishDir "${params.out_dir}/preprocessing", mode: 'copy', overwrite: true
	// Generate Samplesheet from simple listed samples and filter file.
	// Use barcode whitelist from Rui Ye, original locations /volumes/USR1/ruiye/WDR_10XFIX/

	input:
	path wafar_cell_bc
	path wdr_fix_cell_bc
	path icell8_filter
	path icell8_data
	path sample_layout

	output:
		path("SampleSheet.csv")
	script:
	"""
	echo "[Header],,,,,,,,
	Investigator Name,Ryan,,,,,,,
	Date,\$(date +'%d/%m/%Y'),,,,,,,
	Workflow,GenerateFASTQ,,,,,,,
	[Data],,,,,,,,
	Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,Sample_Project,Description" > SampleSheet.csv
	//g" >> SampleSheet.csvMPLESHEET | sed 's/\./_/g' | awk '{print ","$1","$1",plate1,A"NR","substr($3,1,8)","substr($3,10,16)","$1","$1}' | sed -e "s/

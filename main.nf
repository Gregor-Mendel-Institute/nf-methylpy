#!/usr/bin/env nextflow
/*
========================================================================================
             B S - S E Q   M E T H Y L A T I O N   B E S T - P R A C T I C E
========================================================================================
 New Methylation (BS-Seq) Best Practice Analysis Pipeline. Started June 2016.
 #### Homepage / Documentation
 https://github.com/nf-core/methylseq
 #### Authors
 Phil Ewels <phil.ewels@scilifelab.se>
----------------------------------------------------------------------------------------
*/


/*

Simply run this

nextflow run ~/mygit/methylseq/main.nf --reads "*bam" --file_ext bam --fasta ~/TAiR10_ARABIDOPSIS/TAIR10_wholeGenome.fasta --outdir output_folder

*/

/*
 * SET UP CONFIGURATION VARIABLES
 */
params.project = "cegs"
build_index = false
params.outdir = './methylpy'
params.umeth = "ChrC:"
params.fasta = false
params.file_ext = "bam"  // please change this accordingly..
params.snpcall = false

params.name = false
params.clusterOptions = false
params.email = false
params.plaintext_email = false
params.genome = false
params.bismark_index = params.genome ? params.genomes[ params.genome ].bismark ?: false : false
params.bwa_meth_index = params.genome ? params.genomes[ params.genome ].bwa_meth ?: false : false
params.fasta_index = params.genome ? params.genomes[ params.genome ].fasta_index ?: false : false

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
try {
  if( ! nextflow.version.matches(">= $params.nf_required_version") ){
    throw GroovyException('Nextflow version too old')
  }
} catch (all) {
  log.error "====================================================\n" +
            "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
            "  Pipeline execution will continue, but things may break.\n" +
            "  Please run `nextflow self-update` to update Nextflow.\n" +
            "============================================================"
}
// Show a big error message if we're running on the base config and an uppmax cluster
if( workflow.profile == 'standard'){
    if ( "hostname".execute().text.contains('.uppmax.uu.se') ) {
        log.error "====================================================\n" +
                  "  WARNING! You are running with the default 'standard'\n" +
                  "  pipeline config profile, which runs on the head node\n" +
                  "  and assumes all software is on the PATH.\n" +
                  "  ALL JOBS ARE RUNNING LOCALLY and stuff will probably break.\n" +
                  "  Please use `-profile uppmax` to run on UPPMAX clusters.\n" +
                  "============================================================"
    }
}

// Validate inputs
if (params.aligner == 'methylpy'){
  "Running methylpy"
} else if (params.aligner != 'bismark' && params.aligner != 'bwameth'){
    exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'bismark', 'bwameth'"
}
if( params.bismark_index && params.aligner == 'bismark' ){
    bismark_index = Channel
        .fromPath(params.bismark_index)
        .ifEmpty { exit 1, "Bismark index not found: ${params.bismark_index}" }
}
else if( params.bwa_meth_index && params.aligner == 'bwameth' ){
    bwa_meth_indices = Channel
        .fromPath( "${params.bwa_meth_index}*" )
        .ifEmpty { exit 1, "bwa-meth index not found: ${params.bwa_meth_index}" }
}
else if( params.fasta_index && params.aligner == 'bwameth' ){
    fasta_index = file(params.fasta_index)
    if( !fasta_index.exists() ) exit 1, "Fasta index file not found: ${params.fasta_index}"
}
else if( !params.fasta ) {
    exit 1, "No reference genome index specified!"
}
if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}
else if( params.aligner == 'bwameth') {
    exit 1, "No Fasta reference specified! This is required by MethylDackel."
}
multiqc_config = file(params.multiqc_config)

if ( params.fasta && params.aligner == 'methylpy' ){
    genome = file(params.fasta)
    reffol = genome.parent
    refid = genome.baseName
    if( !genome.exists() ) exit 1, "Reference fasta file not found: ${params.fasta}"
    methylpy_indices = Channel
      .fromPath( "$reffol/${refid}_methylpy/${refid}*" )
      .ifEmpty { build_index = true }
      .subscribe onComplete: { checked_genome_index = true }
}


// Validate inputs
if( workflow.profile == 'uppmax' || workflow.profile == 'uppmax_devel' ){
    if ( !params.project ) exit 1, "No UPPMAX project ID found! Use --project"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Library prep presets
params.illumina = true
params.rrbs = false
params.pbat = false
params.single_cell = false
params.epignome = false
params.accel = false
params.zymo = false
params.cegx = false
if(params.pbat){
    params.clip_r1 = 6
    params.clip_r2 = 9
    params.three_prime_clip_r1 = 6
    params.three_prime_clip_r2 = 9
} else if(params.single_cell){
    params.clip_r1 = 6
    params.clip_r2 = 6
    params.three_prime_clip_r1 = 6
    params.three_prime_clip_r2 = 6
} else if(params.epignome){
    params.clip_r1 = 8
    params.clip_r2 = 8
    params.three_prime_clip_r1 = 8
    params.three_prime_clip_r2 = 8
} else if(params.accel || params.zymo){
    params.clip_r1 = 10
    params.clip_r2 = 15
    params.three_prime_clip_r1 = 10
    params.three_prime_clip_r2 = 10
} else if(params.cegx){
    params.clip_r1 = 6
    params.clip_r2 = 6
    params.three_prime_clip_r1 = 2
    params.three_prime_clip_r2 = 2
} else {
    params.clip_r1 = 0
    params.clip_r2 = 0
    params.three_prime_clip_r1 = 0
    params.three_prime_clip_r2 = 0
}

/*
 * Create a channel for input read files
 */
num_files = 1
if ( params.file_ext == 'fastq' ){
  num_files = params.singleEnd ? 1 : 2
}
read_files_processing = Channel
    .fromFilePairs( params.reads, size: num_files )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }

log.info "=================================================="
log.info " nf-core/methylseq : Bisulfite-Seq Best Practice v${params.version}"
log.info "=================================================="
def summary = [:]
summary['Run Name']       = custom_runName ?: workflow.runName
summary['Reads']          = params.reads
summary['Aligner']        = params.aligner
summary['Data Type']      = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Genome']         = params.genome ?: genome
if(params.bismark_index) summary['Bismark Index'] = params.bismark_index
if(params.bwa_meth_index) summary['BWA-Meth Index'] = "${params.bwa_meth_index}*"
else if(params.fasta)    summary['Fasta Ref'] = params.fasta
if(params.rrbs) summary['RRBS Mode'] = 'On'
if(params.relaxMismatches) summary['Mismatch Func'] = "L,0,-${params.numMismatches} (Bismark default = L,0,-0.2)"
if(params.notrim)       summary['Trimming Step'] = 'Skipped'
if(params.pbat)         summary['Trim Profile'] = 'PBAT'
if(params.single_cell)  summary['Trim Profile'] = 'Single Cell'
if(params.epignome)     summary['Trim Profile'] = 'TruSeq (EpiGnome)'
if(params.accel)        summary['Trim Profile'] = 'Accel-NGS (Swift)'
if(params.zymo)         summary['Trim Profile'] = 'Zymo Pico-Methyl'
if(params.cegx)         summary['Trim Profile'] = 'CEGX'
summary['Trim R1'] = params.clip_r1
summary['Trim R2'] = params.clip_r2
summary["Trim 3' R1"] = params.three_prime_clip_r1
summary["Trim 3' R2"] = params.three_prime_clip_r2
summary['Deduplication']  = params.nodedup || params.rrbs ? 'No' : 'Yes'
summary['Directional Mode'] = params.single_cell || params.zymo || params.non_directional ? 'No' : 'Yes'
summary['All C Contexts'] = params.comprehensive ? 'Yes' : 'No'
if(params.mindepth) summary['Minimum Depth'] = params.mindepth
if(params.ignoreFlags) summary['MethylDackel'] = 'Ignoring SAM Flags'
summary['Save Reference'] = params.saveReference ? 'Yes' : 'No'
summary['Save Trimmed']   = params.saveTrimmed ? 'Yes' : 'No'
summary['Save Unmapped']  = params.unmapped ? 'Yes' : 'No'
summary['Save Intermeds'] = params.saveAlignedIntermediates ? 'Yes' : 'No'
summary['Max Memory']     = params.max_memory
summary['Max CPUs']       = params.max_cpus
summary['Max Time']       = params.max_time
summary['Output dir']     = params.outdir
summary['Working dir']    = workflow.workDir
summary['Container']      = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(params.project) summary['UPPMAX Project'] = params.project
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "========================================="


// PREPROCESSING - Build methylpy genome index

if (build_index == true && params.aligner == 'methylpy'){
  process makeMethylpyIndex {
    publishDir path: "${reffol}", mode: 'copy',
            saveAs: {filename -> filename.indexOf(".fai") > 0 ? "$filename" : "${refid}_methylpy/$filename"}

    input:
    file genome
    checked_genome_index

    output:
    file "${genome}.fai" into genome_index
    file "${refid}_*" into built_methylpy_index

    script:
    """
    samtools faidx ${genome}
    methylpy build-reference --input-files ${genome} --output-prefix ${refid} --bowtie2 True
    """
  }
} else if (checked_genome_index == true){
  genome_index = Channel.fromPath("${reffol}/${refid}.fasta.fai")
  built_methylpy_index = Channel.fromPath("${reffol}/${refid}_methylpy/${refid}*")
}

// Step 0, preprocessing input read files

if (params.file_ext == "fastq"){
  read_files_processing.into { read_files_fastqc; read_files_trimming }
} else {
  process reads_preprocess {
    tag "$name"
    storeDir "${params.tmpdir}/rawreads"
    label 'env_picard_small'

    input:
    set val(name), file(reads) from read_files_processing

    output:
    set val(name), file("${name}*fastq") into read_files_fastqc
    set val(name), file("${name}*fastq") into read_files_trimming

    script:
    if (params.singleEnd) {
      if (reads.getExtension() == "sra") {
        """
        fastq-dump $reads
        """
      } else if (reads.getExtension() == "bam") {
        """
        picard SamToFastq I=$reads FASTQ=${name}.fastq VALIDATION_STRINGENCY=LENIENT
        """
      }
    } else {
      if (reads[0].getExtension() == "sra") {
        """
        fastq-dump --split-files $reads
        """
      } else if (reads.getExtension() == "bam") {
        """
        picard SamToFastq I=$reads FASTQ=${name}_1.fastq SECOND_END_FASTQ=${name}_2.fastq VALIDATION_STRINGENCY=LENIENT
        """
      }
    }
  }
}

/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    label 'env_quality'
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file '*_fastqc.{zip,html}' into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}

/*
 * STEP 2 - Trim Galore!
 */
if(params.notrim){
    trimmed_reads = read_files_trimming
    trimgalore_results = Channel.from(false)
} else {
    process trim_galore {
        tag "$name"
        label 'env_quality'
        publishDir "${params.outdir}/trim_galore", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
                else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
                else params.saveTrimmed ? filename : null
            }

        input:
        set val(name), file(reads) from read_files_trimming

        output:
        set val(name), file('*fq.gz') into trimmed_reads
        file "*trimming_report.txt" into trimgalore_results
        file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports

        script:
        c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
        c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
        tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
        tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
        rrbs = params.rrbs ? "--rrbs" : ''
        illumina = params.illumina ? "--illumina" : ''
        non_directional = params.rrbs && params.non_directional ? "--non_directional" : ''
        if (params.singleEnd) {
            """
            trim_galore --fastqc --gzip $rrbs $c_r1 $tpc_r1 $reads
            """
        } else {
            """
            trim_galore --paired --fastqc --gzip $rrbs $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
            """
        }
    }
}

// 3 - align with methylpy

if(params.aligner == 'methylpy'){
  process methylpy_align {
    tag "$name"
    publishDir "${params.outdir}", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".bam") > 0) "alignedBams/$filename"
            else if (filename =~ '^allc' ) "allc/$filename"
            else if (filename =~ '^conversion' ) "info/$filename"
            else if (filename =~ '^log' ) "info/log.${name}.txt"
          }
    label 'env_methylpy'

    input:
    set val(name), file(reads) from trimmed_reads
    file fasta
    file(meth_index) from built_methylpy_index.collect()
    file(meth_genome_index) from genome_index.collect()

    output:
    set val(name), file("*processed_reads_no_clonal.bam") into bam_aligned
    set val(name), file("allc_*tsv.gz") into allc
    set val(name), file("conversion_rate_${name}.txt") into conv_rate
    set val(name), file("log.txt") into log_file

    script:
    if (params.singleEnd) {
        """
        export TMPDIR="${params.tmpdir}"
        methylpy single-end-pipeline --read-files ${reads} --sample $name --forward-ref ${refid}_f  --reverse-ref ${refid}_r  --ref-fasta  $fasta   --num-procs ${task.cpus}  --remove-clonal True --binom-test True  --unmethylated-control ${params.umeth} --java-options="-Djava.io.tmpdir=${params.tmpdir}" > log.txt 2>&1
        cat log.txt | grep "non-conversion rate" > conversion_rate_${name}.txt
        """
    } else {
        """
        export TMPDIR="${params.tmpdir}"
        methylpy paired-end-pipeline --read1-files ${reads[0]}  --read2-files ${reads[1]}  --sample ${name}  --forward-ref ${refid}_f  --reverse-ref ${refid}_r  --ref-fasta  $fasta  --num-procs ${task.cpus} --remove-clonal True --binom-test True  --unmethylated-control ${params.umeth} --java-options="-Djava.io.tmpdir=${params.tmpdir}" > log.txt 2>&1
        cat log.txt | grep "non-conversion rate" > conversion_rate_${name}.txt
        """
    }
  }

  if (params.snpcall) {
    bam_aligned.into { bams_index ; bams_snpcall }
  } else {
    bams_index = bam_aligned
  }

  process bam_index {
    tag "$name"
    publishDir "${params.outdir}/alignedBams", mode: 'copy'
    label 'env_picard_small'

    input:
    set val(name), file(bam) from bams_index

    output:
    set val(name), file("${bam}.bai") into aligned_bam_index

    script:
    """
    samtools index $bam
    """
  }

}

/*
Make hdf5 files for the allc
*/
process make_hdf5 {
  tag { "$name" }
  publishDir "${params.outdir}/hdf5", mode: 'copy'
  label 'env_pybshap'

  input:
  set val(name), file(allc) from allc

  output:
  set val(name), file("*hdf5") into hdf5_out

  script:
  """
  bshap methylation_percentage -i $allc -a new -b ChrC,1,100 -o temp -v
  """
}

/*
SNP calling from the methylpy
*/
if (params.snpcall){

  process add_readgroups {
    tag "$name"
    label 'env_picard'

    input:
    set val(name), file(bam) from bams_snpcall

    output:
    set val(name), file("${name}.modified.bam"), file("${name}.modified.bam.bai") into modifiedbam

    script:
    """
    picard AddOrReplaceReadGroups INPUT=$bam OUTPUT=${name}.modified.bam ID=${name} LB=${name} PL=illumina PU=none SM=${name}
    samtools index ${name}.modified.bam
    """
  }


  process do_realignindel {
    tag "$name"
    label 'env_snpcall_small'

    input:
    set val(name), file(bam), file(bam_index) from modifiedbam

    output:
    set val(name), file("${name}.realignedBam.bam"), file("${name}.realignedBam.bam.bai") into realignedbam

    script:
    """
    java -Djava.io.tmpdir=${params.tmpdir} -jar \$EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reffol/${refid}.fasta -I $bam -o ${name}.forIndelRealigner.intervals
    java -Djava.io.tmpdir=${params.tmpdir} -jar \$EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R $reffol/${refid}.fasta -I $bam -targetIntervals ${name}.forIndelRealigner.intervals -o ${name}.realignedBam.bam
    samtools index ${name}.realignedBam.bam
    """
  }

  process do_snpcall {
    tag "$name"
    label 'env_snpcall'

    input:
    set val(name), file(bam), file(bam_index) from realignedbam

    output:
    set val(name), file("${name}.vcf") into vcffile

    script:
    """
    java -Djava.io.tmpdir=${params.tmpdir} -jar \$EBROOTGATK/GenomeAnalysisTK.jar \
    -T HaplotypeCaller -R $reffol/${refid}.fasta \
    -I $bam \
    -o ${name}.vcf \
    -nct ${task.cpus} \
    --output_mode EMIT_ALL_SITES \
    --genotyping_mode DISCOVERY \
    """
  }

  process get_snps_from_vcf {
    tag "$name"
    publishDir "${params.outdir}/vcfbed", mode: 'copy'
    label 'env_snpcall_small'

    input:
    set val(name), file(vcf) from vcffile

    output:
    set val(name), file("${name}.snpvcf.bed") into snpbed

    script:
    """
    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT[\t%GT][\t%AD]\n" $vcf | awk '\$3 !~ /C|G/ && length(\$3) == 1 && length(\$4) == 1 && \$4 !~ /T/ {print \$1 "\t" \$2 "\t" \$5 "\t" \$6}'  > ${name}.snpvcf.bed
    """
  }
}

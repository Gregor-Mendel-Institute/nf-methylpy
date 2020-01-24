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


RUN:

nextflow run ~/mygit/methylseq/methylpyDiffGenomes.nf --inputfiles input_csv --outdir methylpy


 * SET UP CONFIGURATION VARIABLES
 */

params.inputfiles = false
params.outdir = './methylpy'
params.file_ext = false

params.project = "cegs"
params.aligner = "methylpy"
if( params.aligner != "methylpy" ) exit 1, "This pipeline has been written only for methylpy, please choose methylpy as aligner"
params.umeth = "ChrC:"

params.name = false
params.clusterOptions = false
params.email = false
params.plaintext_email = false
params.genome = false

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
input_genomes = Channel
    .fromPath( params.inputfiles )
    .ifEmpty { exit 1, "Provide a tab separated table indicating the read pairs and reference fasta."}
    .splitCsv(sep: '\t')
    .unique{ row -> [row[0], row[1]] }

input_reads = Channel
  .fromPath( params.inputfiles )
  .splitCsv(sep: '\t')

log.info "=================================================="
log.info " nf-core/methylseq : Bisulfite-Seq Best Practice v${params.version}"
log.info "=================================================="
def summary = [:]
summary['Run Name']       = custom_runName ?: workflow.runName
summary['Input file']     = params.inputfiles
summary['Aligner']        = params.aligner
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
process makeMethylpyIndex {
  tag { "$accID" }
  storeDir "$reffol/${refid}_methylpy"

  input:
  set val(accID), genome_file, reads from input_genomes

  output:
  set val(accID), genome, file("${refid}*") into methylpy_indices

  script:
  genome = file(genome_file)
  reffol = genome.parent
  refid = genome.baseName
  f_ref = "$reffol/${refid}_methylpy/${refid}_f"
  r_ref = "$reffol/${refid}_methylpy/${refid}_r"
  """
  samtools faidx ${genome}
  methylpy build-reference --input-files ${genome} --output-prefix ${refid}
  """
}


// Step 0, preprocessing input read files
if (params.file_ext == "fastq"){
  input_reads.into { read_files_fastqc; read_files_trimming }
} else {
  process identify_libraries{
    label 'env_picard_small'
    tag { "${accID}_${read_files.baseName}" }

    input:
    set val(accID), genome, reads from input_reads

    output:
    set val(accID), reads, stdout into read_files_processing

    script:
    read_files = file(reads)
    file_ext = read_files.getExtension()
    if (file_ext == "sra"){
      """
      fastq-dump -I -X 1 -Z --split-spot $read_files  | awk ' NR % 2 == 1 {print substr(\$1,length(\$1),1) } ' | uniq | wc -l
      """
    } else if (file_ext == "bam"){
      """
      (samtools view -h $read_files | head -n 100 | samtools view -f 0x1 | wc -l | awk '{print \$0 + 1 }' ) || true
      """
    }
  }

  process reads_preprocess {
    tag { "${accID}_${read_files.baseName}" }
    storeDir "${params.tmpdir}/rawreads"
    label 'env_picard_small'

    input:
    set val(accID), reads, library_id from read_files_processing

    output:
    set val(accID), file("${prefix}*fastq"), val(library) into read_files_fastqc
    set val(accID), file("${prefix}*fastq"), val(library) into read_files_trimming

    script:
    read_files = file(reads)
    file_ext = read_files.getExtension()
    prefix = read_files.baseName.toString() - ~/(\.sra)?(\.bam)?$/
    library = library_id.toInteger()
    if (library == 1) {
      if (file_ext == "sra") {
        """
        fastq-dump $read_files
        """
      } else if (file_ext == "bam") {
        """
        java -jar \${EBROOTPICARD}/picard.jar SamToFastq I=$read_files FASTQ=${prefix}.fastq VALIDATION_STRINGENCY=LENIENT
        """
      }
    } else {
      if (file_ext == "sra") {
        """
        fastq-dump --split-files $read_files
        """
      } else if (file_ext == "bam") {
        """
        java -jar \${EBROOTPICARD}/picard.jar SamToFastq I=$read_files FASTQ=${prefix}_1.fastq SECOND_END_FASTQ=${prefix}_2.fastq VALIDATION_STRINGENCY=LENIENT
        """
      }
    }
  }
}


/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag { "${accID}_$reads" }
    label 'env_quality'
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(accID), file(reads), val(library_id) from read_files_fastqc

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
        tag { "${accID}_$reads" }
        label 'env_quality'
        publishDir "${params.outdir}/trim_galore", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
                else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
                else params.saveTrimmed ? filename : null
            }

        input:
        set val(accID), file(reads), val(library_id) from read_files_trimming

        output:
        set val(accID), file('*fq.gz'), val(library_id) into trimmed_reads
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
        if (library_id == 1) {
            """
            trim_galore --fastqc --gzip $illumina $rrbs $c_r1 $tpc_r1 $reads
            """
        } else {
            """
            trim_galore --paired --fastqc --gzip $illumina $rrbs $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
            """
        }
    }
}

// 3 - align with methylpy
input_reads_methylpy = trimmed_reads
    .combine( methylpy_indices , by: 0)

process methylpy_align {
  tag { "${accID}_$reads" }
  publishDir "${params.outdir}", mode: 'copy',
      saveAs: {filename ->
          if (filename.indexOf(".bam") > 0) "alignedBams/$filename"
          else if (filename =~ '^allc' ) "allc/$filename"
          else if (filename =~ '^conversion' ) "info/$filename"
          else if (filename =~ '^log' ) "info/log.${name}.txt"
        }
  label 'env_methylpy'

  input:
  set val(accID), file(reads), val(library_id), file(genome), file(meth_index) from input_reads_methylpy

  output:
  set val(prefix), file("*processed_reads_no_clonal.bam") into bam_aligned
  set val(prefix), file("allc_*tsv.gz*") into allc
  set val(prefix), file("conversion_rate_${prefix}.txt") into conv_rate

  script:
  reffol = genome.parent
  refid = genome.baseName
  if (library_id == 1) {
      prefix = reads.toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
      """
      export TMPDIR="${params.tmpdir}"
      methylpy single-end-pipeline --read-files ${reads} --sample $prefix --forward-ref ${refid}_f  --reverse-ref ${refid}_r  --ref-fasta  $genome   --num-procs ${task.cpus}  --remove-clonal True   --binom-test True  --unmethylated-control ${params.umeth} --java-options="-Djava.io.tmpdir=${params.tmpdir}" > log.txt 2>&1
      cat log.txt | grep "non-conversion rate" > conversion_rate_${prefix}.txt
      """
  } else {
      prefix = reads[0].toString() - ~/(_1)?(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
      """
      export TMPDIR="${params.tmpdir}"
      methylpy paired-end-pipeline --read1-files ${reads[0]}  --read2-files ${reads[1]}  --sample $prefix  --forward-ref ${refid}_f  --reverse-ref ${refid}_r  --ref-fasta  $genome  --num-procs ${task.cpus}  --remove-clonal True  --binom-test True  --unmethylated-control ${params.umeth} --java-options="-Djava.io.tmpdir=${params.tmpdir}" > log.txt 2>&1
      cat log.txt | grep "non-conversion rate" > conversion_rate_${prefix}.txt
      """
  }
}

process bam_index {
  tag { "${prefix}" }
  publishDir "${params.outdir}/alignedBams", mode: 'copy'
  label 'env_picard_small'

  input:
  set val(prefix), file(bam) from bam_aligned

  output:
  set val(prefix), file("*processed_reads_no_clonal.bam.bai") into aligned_bam_index

  script:
  """
  samtools index $bam
  """
}

process make_hdf5 {
  tag { "${prefix}" }
  publishDir "${params.outdir}/hdf5", mode: 'copy'
  label 'env_pybshap'

  input:
  set val(prefix), file(allc) from allc

  output:
  set val(prefix), file("*hdf5") into hdf5_out

  script:
  """
  bshap methylation_percentage -i $allc -a new -b Chr1,1,100 -o temp -v
  """
}

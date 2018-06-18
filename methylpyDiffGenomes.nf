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
 * SET UP CONFIGURATION VARIABLES
 */

params.inputfiles = false
params.outdir = './methylpy'

params.project = "cegs"
params.aligner = "methylpy"
if( params.aligner != "methylpy" ) exit 1, "This pipeline has been written only for methylpy, please choose methylpy as aligner"
params.umeth = "ChrC:"
params.tmpdir = "/lustre/scratch/users/rahul.pisupati/tempFiles/"

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
else if( params.aligner == 'bwameth') {
    exit 1, "No Fasta reference specified! This is required by MethylDackel."
}
multiqc_config = file(params.multiqc_config)

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
  methylpy build-reference --input-files ${genome} --output-prefix ${refid} --bowtie2 True
  """
}


// Step 0, preprocessing input read files
if (params.file_ext == "fastq"){
  input_reads.into { read_files_fastqc; read_files_trimming }
} else {
  process identify_libraries{
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

    input:
    set val(accID), reads, val(library_id) from read_files_processing

    output:
    set val(accID), file("${prefix}*fastq"), val(library_id) into read_files_fastqc
    set val(accID), file("${prefix}*fastq"), val(library_id) into read_files_trimming

    script:
    read_files = file(reads)
    file_ext = read_files.getExtension()
    prefix = read_files.baseName.toString() - ~/(\.sra)?(\.bam)?$/
    if (library_id == 1) {
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
        if (params.singleEnd) {
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
        }

  input:
  set val(accID), file(reads), val(library_id), file(genome), file(meth_index) from input_reads_methylpy

  output:
  set val(accID), file("*processed_reads_no_clonal.bam") into bam_aligned
  set val(accID), file("allc_*tsv.gz*") into allc
  set val(accID), file("conversion_rate_${prefix}.txt") into conv_rate

  script:
  reffol = genome.parent
  refid = genome.baseName
  if (library_id == 1) {
      prefix = reads.toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
      """
      export TMPDIR="${params.tmpdir}"
      methylpy single-end-pipeline --read-files ${reads} --sample $prefix --forward-ref ${refid}_f  --reverse-ref ${refid}_r  --ref-fasta  $genome   --num-procs ${task.cpus}  --remove-clonal True  --path-to-picard \$EBROOTPICARD  --binom-test True  --unmethylated-control ${params.umeth} --java-options="-Djava.io.tmpdir=${params.tmpdir}" > log.txt 2>&1
      cat log.txt | grep "non-conversion rate" > conversion_rate_${prefix}.txt
      """
  } else {
      prefix = reads[0].toString() - ~/(_1)?(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
      """
      export TMPDIR="${params.tmpdir}"
      methylpy paired-end-pipeline --read1-files ${reads[0]}  --read2-files ${reads[1]}  --sample ${prefix}  --forward-ref ${refid}_f  --reverse-ref ${refid}_r  --ref-fasta  $genome  --num-procs ${task.cpus}  --remove-clonal True  --path-to-picard \${EBROOTPICARD}  --binom-test True  --unmethylated-control ${params.umeth} --java-options="-Djava.io.tmpdir=${params.tmpdir}" > log.txt 2>&1
      cat log.txt | grep "non-conversion rate" > conversion_rate_${prefix}.txt
      """
  }
}

process bam_index {
  tag { "${accID}_$bam" }
  publishDir "${params.outdir}/alignedBams", mode: 'copy'

  input:
  set val(accID), file(bam) from bam_aligned

  output:
  set val(accID), file("*processed_reads_no_clonal.bam.bai") into aligned_bam_index

  script:
  """
  samtools index $bam
  """
}

/*
 * STEP 8 - Qualimap
process qualimap {
    tag "${bam.baseName}"
    publishDir "${params.outdir}/qualimap", mode: 'copy'

    input:
    file bam from bam_dedup_qualimap

    output:
    file "${bam.baseName}_qualimap" into qualimap_results

    script:
    gcref = params.genome == 'GRCh37' ? '-gd HUMAN' : ''
    gcref = params.genome == 'GRCm38' ? '-gd MOUSE' : ''
    """
    samtools sort $bam -o ${bam.baseName}.sorted.bam
    qualimap bamqc $gcref \\
        -bam ${bam.baseName}.sorted.bam \\
        -outdir ${bam.baseName}_qualimap \\
        --collect-overlap-pairs \\
        --java-mem-size=${task.memory.toGiga()}G \\
        -nt ${task.cpus}
    """
}
*/

/*
 * Parse software version numbers
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo "$params.version" &> v_ngi_methylseq.txt
    echo "$workflow.nextflow.version" &> v_nextflow.txt
    bismark_genome_preparation --version &> v_bismark_genome_preparation.txt
    fastqc --version &> v_fastqc.txt
    cutadapt --version &> v_cutadapt.txt
    trim_galore --version &> v_trim_galore.txt
    bismark --version &> v_bismark.txt
    deduplicate_bismark --version &> v_deduplicate_bismark.txt
    bismark_methylation_extractor --version &> v_bismark_methylation_extractor.txt
    bismark2report --version &> v_bismark2report.txt
    bismark2summary --version &> v_bismark2summary.txt
    samtools --version &> v_samtools.txt
    bwa &> v_bwa.txt 2>&1 || true
    bwameth.py --version &> v_bwameth.txt
    picard MarkDuplicates --version &> v_picard_markdups.txt 2>&1 || true
    MethylDackel --version &> v_methyldackel.txt
    qualimap --version &> v_qualimap.txt
    multiqc --version &> v_multiqc.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}
*/


/*
 * STEP 9 - MultiQC
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file ('fastqc/*') from fastqc_results.toList()
    file ('trimgalore/*') from trimgalore_results.toList()
    file ('bismark/*') from bismark_align_log_3.toList()
    file ('bismark/*') from bismark_dedup_log_3.toList()
    file ('bismark/*') from bismark_splitting_report_3.toList()
    file ('bismark/*') from bismark_mbias_3.toList()
    file ('bismark/*') from bismark_reports_results.toList()
    file ('bismark/*') from bismark_summary_results.toList()
    file ('samtools/*') from flagstat_results.flatten().toList()
    file ('samtools/*') from samtools_stats_results.flatten().toList()
    file ('picard/*') from picard_results.flatten().toList()
    file ('methyldackel/*') from methyldackel_results.flatten().toList()
    file ('qualimap/*') from qualimap_results.toList()
    file ('software_versions/*') from software_versions_yaml.toList()

    output:
    file "*_report.html" into multiqc_report
    file "*_data"
    file '.command.err' into multiqc_stderr

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}
*/

/*
 * Completion e-mail notification
 *
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/methylseq] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/methylseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = params.version
    email_fields['runName'] = workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    email_fields['summary']['Container'] = workflow.container
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/methylseq] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/methylseq] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Switch the embedded MIME images with base64 encoded src
    ngimethylseqlogo = new File("$baseDir/assets/methylseq_logo.png").bytes.encodeBase64().toString()
    scilifelablogo = new File("$baseDir/assets/SciLifeLab_logo.png").bytes.encodeBase64().toString()
    ngilogo = new File("$baseDir/assets/NGI_logo.png").bytes.encodeBase64().toString()
    email_html = email_html.replaceAll(~/cid:ngimethylseqlogo/, "data:image/png;base64,$ngimethylseqlogo")
    email_html = email_html.replaceAll(~/cid:scilifelablogo/, "data:image/png;base64,$scilifelablogo")
    email_html = email_html.replaceAll(~/cid:ngilogo/, "data:image/png;base64,$ngilogo")

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[nf-core/methylseq] Pipeline Complete"

    if(!workflow.success){
        if( workflow.profile == 'standard'){
            if ( "hostname".execute().text.contains('.uppmax.uu.se') ) {
                log.error "====================================================\n" +
                        "  WARNING! You are running with the default 'standard'\n" +
                        "  pipeline config profile, which runs on the head node\n" +
                        "  and assumes all software is on the PATH.\n" +
                        "  This is probably why everything broke.\n" +
                        "  Please use `-profile uppmax` to run on UPPMAX clusters.\n" +
                        "============================================================"
            }
        }
    }

}
*/

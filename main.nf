#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/qtlmap
========================================================================================
 nf-core/qtlmap Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/qtlmap
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
    log.info"""
    =======================================================
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\'
        |\\ | |__  __ /  ` /  \\ |__) |__         }  {
        | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                              `._,._,\'

     nf-core/qtlmap v${workflow.manifest.version}
    =======================================================

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/qtlmap --reads '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --genome                      Name of iGenomes reference
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.

    Options:
      --singleEnd                   Specifies that the input is single end reads

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// TODO nf-core: Add any reference files that are needed
// Configurable reference genomes
// fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
// if ( params.fasta ){
//     fasta = file(params.fasta)
//     if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
// }
//
// NOTE - THIS IS NOT USED IN THIS PIPELINE, EXAMPLE ONLY
// If you want to use the above in a process, define the following:
//   input:
//   file fasta from fasta
//


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


if( workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
  // Check workDir/outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}

// Stage config files
// ch_multiqc_config = Channel.fromPath(params.multiqc_config)
// ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")

/*
 * Create a channel for input read files
 */ 
 if(params.readPaths){
     if(params.singleEnd){
         Channel
             .from(params.readPaths)
             .map { row -> [ row[0], [file(row[1][0])]] }
             .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
             .into { read_files_fastqc; read_files_trimming }
     } else {
         Channel
             .from(params.readPaths)
             .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
             .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
             .into { read_files_fastqc; read_files_trimming }
     }
 } else {
     Channel
         .fromPath( params.expression_matrix)
         .ifEmpty { exit 1, "Cannot find any expression_matrix file: ${params.expression_matrix}\nNB: Path needs to be enclosed in quotes!" }
         .set { expression_matrix_create_QTLTools_input}
     Channel
         .fromPath( params.sample_metadata )
         .ifEmpty { exit 1, "Cannot find any sample metadata file: ${params.sample_metadata}\nNB: Path needs to be enclosed in quotes!" }
         .set { sample_metadata_create_QTLTools_input}
     Channel
         .from( params.conditions )
         .set { conditions_create_QTLTools_input}
         
 }


// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

nf-core/qtlmap v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']        = 'nf-core/qtlmap'
summary['Pipeline Version']     = workflow.manifest.version
summary['Run Name']             = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here

summary['Expression Matrix']    = params.expression_matrix
summary['Gene Metadata']        = params.gene_metadata
summary['Sample Metadata']      = params.sample_metadata
summary['Genotype file']        = params.vcf_file
summary['Variant Info']         = params.variant_info
summary['Cis distance']         = params.cisdistance
summary['Minimum Cis variants'] = params.mincisvariant
summary['Max Memory']           = params.max_memory
summary['Max CPUs']             = params.max_cpus
summary['Max Time']             = params.max_time
summary['Output dir']           = params.outdir
summary['Working dir']          = workflow.workDir
summary['Container Engine']     = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']         = "$HOME"
summary['Current user']         = "$USER"
summary['Current path']         = "$PWD"
summary['Working dir']          = workflow.workDir
summary['Output dir']           = params.outdir
summary['Script dir']           = workflow.projectDir
summary['Config Profile']       = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']        = params.awsregion
   summary['AWS Queue']         = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "========================================="


def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-qtlmap-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/qtlmap Workflow Summary'
    section_href: 'https://github.com/nf-core/qtlmap'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}


/*
 * Parse software version numbers
 */
// process get_software_versions {

//     output:
//     file 'software_versions_mqc.yaml' into software_versions_yaml

//     script:
//     // TODO nf-core: Get all tools to print their version number here
//     """
//     echo $workflow.manifest.version > v_pipeline.txt
//     echo $workflow.nextflow.version > v_nextflow.txt
//     fastqc --version > v_fastqc.txt
//     multiqc --version > v_multiqc.txt
//     scrape_software_versions.py > software_versions_mqc.yaml
//     """
// }



/*
 * STEP 1 - Generate QTLTools input files
 */
process create_QTLTools_input {
    tag "${expression_matrix.baseName}"
    publishDir "${params.outdir}/qtl_input", mode: 'copy'

    input:
    file expression_matrix from expression_matrix_create_QTLTools_input
    file sample_metadata from sample_metadata_create_QTLTools_input

    output:
    file "${params.quantification_method}/*.bed" into condition_beds
    file "${params.quantification_method}/*.sample_names.txt" into condition_samplenames 

    script:
    """
    array_data_to_QTLtools_input.R \\
        -g "${params.gene_metadata}" \\
        -s "${sample_metadata}" \\
        -e "${expression_matrix}" \\
        -v "${params.variant_info}" \\
        --qtlutils "${params.eqtl_utils}" \\
        -o "${params.quantification_method}" \\
        -c ${params.cisdistance} \\
        -m ${params.mincisvariant}
    """
}

/*
 * STEP 2 - Compres and index input bed file
 */ 
process compress_bed {
    tag "${bed_file.baseName}"
    publishDir "${params.outdir}/compressed_bed", mode: 'copy'

    input:
    file bed_file from condition_beds.flatMap()

    output:
    file "${bed_file}.gz" into compressed_beds_perform_pca, compressed_beds_run_nominal, compressed_beds_run_permutation
    file "${bed_file}.gz.tbi" into compressed_bed_indexes_perform_pca, compressed_bed_indexes_run_nominal, compressed_bed_indexes_run_permutation

    script:
    """
    bgzip $bed_file && tabix -p bed ${bed_file}.gz
    """
}

/*
 * STEP 3 - Extract samples from vcf
 */
process extract_samples {
    tag "${sample_names.simpleName}"
    publishDir "${params.outdir}/vcf", mode: 'copy'

    input:
    file sample_names from condition_samplenames.flatMap()

    output:
    file "${sample_names.simpleName}.vcf.gz" into vcfs_extract_variant_info, vcfs_perform_pca, vcfs_run_nominal, vcfs_run_permutation
    file "${sample_names.simpleName}.vcf.gz.csi" into vcf_indexes_perform_pca

    script:
    """
    bcftools view -S $sample_names ${params.vcf_file} -Oz -o ${sample_names.simpleName}.vcf.gz
    bcftools index ${sample_names.simpleName}.vcf.gz
    """
}

/*
 * STEP 4 - Extract variant information from VCF
 */
process extract_variant_info {
    tag "${vcf.simpleName}"
    publishDir "${params.outdir}/final", mode: 'copy'

    input:
    file vcf from vcfs_extract_variant_info
        // vcf = "processed/{study}/qtltools/input/{annot_type}/vcf/{condition}.vcf.gz",
    
    output:
    file "${vcf.simpleName}.variant_information.txt.gz"
        // var_info = "processed/{study}/qtltools/output/{annot_type}/final/{condition}.variant_information.txt.gz"

    script:
    if (params.is_imputed) {
        """
        set +o pipefail; bcftools +fill-tags $vcf | bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%TYPE\\t%AC\\t%AN\\t%MAF\\t%R2\\n' | gzip > ${vcf.simpleName}.variant_information.txt.gz
        """
    } else {
        """
        set +o pipefail; bcftools +fill-tags $vcf | bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%TYPE\\t%AC\\t%AN\\t%MAF\\tNA\\n' | gzip > ${vcf.simpleName}.variant_information.txt.gz
        """
    }
}

// TODO: try to use each input repeater for permutations

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/qtlmap] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/qtlmap] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
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
          log.info "[nf-core/qtlmap] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/qtlmap] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[nf-core/qtlmap] Pipeline Complete"

}

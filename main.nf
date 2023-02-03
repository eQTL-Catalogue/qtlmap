#!/usr/bin/env nextflow
/*
========================================================================================
                          eQTL-Catalogue/qtlmap
========================================================================================
 eQTL-Catalogue/qtlmap Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/eQTL-Catalogue/qtlmap
----------------------------------------------------------------------------------------
*/
nextflow.enable.dsl=2

def helpMessage() {
    log.info"""
    =======================================================
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\'
        |\\ | |__  __ /  ` /  \\ |__) |__         }  {
        | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                              `._,._,\'

     eQTL-Catalogue/qtlmap v${workflow.manifest.version}
    =======================================================

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf\
     -profile tartu_hpc\
     --studyFile testdata/multi_test.tsv\
     --vcf_has_R2_field FALSE\
     --varid_rsid_map_file testdata/varid_rsid_map.tsv.gz\
     --n_batches 25

    Mandatory arguments:
      --studyFile                   Path to the TSV file containing pipeline inputs (VCF, expression matrix, metadata)

    Executions:
      --n_batches                   Number of parallel batches to run QTL Mapping per sample, must exceed the number of chromosomes (default: 400)
      --vcf_has_R2_field            Does the genotype VCF file contain R2 value in the INFO field? (default: true)
      --run_permutation             Calculate permuation p-values for each phenotype group (group_id in the phenotype metadata file) (default: false)
      --run_nominal                 Calculate nominal p-values for each phenotype group (group_id in the phenotype metadata file) (default: true)
      --n_permutations              Number of permutations to be performed per gene when run_permutation = true (default: 1000)

    QTL mapping:
      --cis_window                  The window where to search for associated variants around the phenotype (default: 1000000)
      --n_geno_pcs                  Number of genotype matrix principal components included as covariates in QTL analysis (default: 6).
      --n_pheno_pcs                 Number of phenotype matrix principal components included as covariates in QTL analysis (default: 6).
      --mincisvariant               Minimal numner of variants needed to be found in cis_window of each phenotype (default: 5)
      --covariates                  Comma-separated list of additional covariates included in the analysis (e.g. sex, age, batch). Columns with the exact same names should exist in the sample metadata file. 

    Fine mapping (SuSiE)
      --run_susie                   Perform eQTL fine mapping with SuSiE
      --vcf_genotype_field          Field in the VCF file that is used to construct the dosage matrix. Valid options are GT and DS (default: GT). 
      --write_full_susie            If 'true' then full SuSiE output will not be written to disk (default: true). 
                                    Setting this to 'false' will apply credible set connected component based filtering to SuSiE results. 
                                    This helps to reduce the size of SuSiE output for molecular traits with many correlated sub-phenotypes (e.g. Leafcutter splice-junctions).

    Format results:
      --reformat_sumstats          Add rsid and median TPM columns to the nominal summary statistics files and perform additional formatting to make the files compatible with the eQTL Catalogue (default: true)
      --varid_rsid_map_file         TSV file mapping variant ids in CHR_POS_REF_ALT format to rsids from dbSNP.

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

/*
 * Create a channel for input files
 */ 

//molecular trait data input data
Channel.fromPath(params.studyFile)
    .ifEmpty { error "Cannot find studyFile file in: ${params.studyFile}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.qtl_subset, file(row.count_matrix), file(row.pheno_meta), file(row.sample_meta)]}
    .set { study_file_ch }

// Separate channel for the VCF file
Channel.fromPath(params.studyFile)
    .ifEmpty { error "Cannot find studyFile file in: ${params.studyFile}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.qtl_subset, file(row.vcf) ]}
    .set { vcf_file_ch }

//Another one for the TPM file that is only needed in the end
Channel.fromPath(params.studyFile)
    .ifEmpty { error "Cannot find studyFile file in: ${params.studyFile}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.qtl_subset, file(row.tpm_file)]}
    .set { tpm_file_ch }

Channel.fromPath(params.varid_rsid_map_file)
    .ifEmpty { error "Cannot find varid_rsid_map_file file in: ${params.varid_rsid_map_file}" }
    .set { rsid_map_ch }

// Batch channel
batch_ch = Channel.of(1..params.n_batches)

// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

eQTL-Catalogue/qtlmap v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']        = 'eQTL-Catalogue/qtlmap'
summary['Pipeline Version']     = workflow.manifest.version
summary['Run Name']             = custom_runName ?: workflow.runName
summary['Study file']           = params.studyFile
summary['cis window']           = params.cis_window
summary['Min # cis variants']   = params.mincisvariant
summary['VCF has R2 field']     = params.vcf_has_R2_field
summary['Permutation run']      = params.run_permutation
summary['# of permutations']    = params.n_permutations
summary['Nominal run']          = params.run_nominal
summary['# of batches']         = params.n_batches
summary['# of phenotype PCs']   = params.n_pheno_pcs
summary['# of genotype PCs']    = params.n_geno_pcs
summary['Additonal covariates'] = params.covariates
summary["Run SuSiE"]            = params.run_susie
summary["Write full SuSiE"]     = params.write_full_susie
summary["VCF genotype field"]   = params.vcf_genotype_field
summary['Max Memory']           = params.max_memory
summary['Max CPUs']             = params.max_cpus
summary['Max Time']             = params.max_time
summary['Output dir']           = params.outdir
summary['Working dir']          = workflow.workDir
summary['Container Engine']     = workflow.containerEngine
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

include { vcf_set_variant_ids } from './modules/vcf_set_variant_ids'
include { extract_lead_cc_signal } from './modules/extract_cc_signal'
include { extract_variant_info } from './modules/extract_variant_info'
include { extract_variant_info as extract_variant_info2 } from './modules/extract_variant_info'
include { prepare_molecular_traits; compress_bed; make_pca_covariates } from './modules/prepare_molecular_traits'
include { extract_samples_from_vcf } from './modules/extract_samples_from_vcf'
include { run_permutation; merge_permutation_batches; run_nominal; merge_nominal_batches; sort_qtltools_output} from './modules/map_qtls'
include { join_rsids_var_info; reformat_sumstats; tabix_index} from './modules/reformat_sumstats'
include { vcf_to_dosage } from './modules/vcf_to_dosage'
include { run_susie; merge_susie; sort_susie; extract_cs_variants; merge_cs_sumstats } from './modules/susie'

workflow {

    // Prepare input data for QTL mapping
    vcf_set_variant_ids(vcf_file_ch)
    extract_variant_info(vcf_set_variant_ids.out)
    prepare_molecular_traits(study_file_ch.join(extract_variant_info.out))
    compress_bed(prepare_molecular_traits.out.bed_file)
    extract_samples_from_vcf(vcf_set_variant_ids.out.join(prepare_molecular_traits.out.sample_names))
    make_pca_covariates(prepare_molecular_traits.out.pheno_cov.join(extract_samples_from_vcf.out.vcf))

    //Run nominal and permutation passes
    qtlmap_input_ch = compress_bed.out
      .join(extract_samples_from_vcf.out.vcf)
      .join(extract_samples_from_vcf.out.index)
      .join(make_pca_covariates.out)
    
    //Permutation pass
    if( params.run_permutation ){
      run_permutation(batch_ch, qtlmap_input_ch)
      merge_permutation_batches( run_permutation.out.groupTuple(size: params.n_batches, sort: true) )
    }
    //Nominal pass
    if( params.run_nominal ){
      run_nominal(batch_ch, qtlmap_input_ch)
      merge_nominal_batches( run_nominal.out.groupTuple(size: params.n_batches, sort: true) )
      sort_qtltools_output( merge_nominal_batches.out )

      //Reformat sumstats
      if( params.reformat_sumstats ){
        extract_variant_info2(extract_samples_from_vcf.out.vcf)
        join_rsids_var_info( extract_variant_info2.out, rsid_map_ch.collect() )
        
        reformat_input_ch = sort_qtltools_output.out
          .join(join_rsids_var_info.out)
          .join(prepare_molecular_traits.out.pheno_meta)
          .join(tpm_file_ch)

        reformat_sumstats( reformat_input_ch )
        tabix_index(reformat_sumstats.out)
      }
    }

    //Run SuSiE
    if( params.run_permutation & params.run_susie ){
      vcf_to_dosage(extract_samples_from_vcf.out.vcf)
      susie_ch = study_file_ch
        .join(merge_permutation_batches.out)
        .join(make_pca_covariates.out)
        .join(vcf_to_dosage.out)
      run_susie(susie_ch, batch_ch)
      merge_susie( run_susie.out.groupTuple(size: params.n_batches, sort: true) )
      sort_susie( merge_susie.out )
    }

    //Extract credible set variants from the full summary statistics
    if (params.run_nominal & params.run_permutation & params.run_susie & params.reformat_sumstats){
      extract_lead_cc_signal(sort_susie.out.join(tabix_index.out))
      extract_cs_variants( sort_susie.out.join(tabix_index.out) )
      merge_cs_sumstats( extract_cs_variants.out )
    }

    
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[eQTL-Catalogue/qtlmap] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[eQTL-Catalogue/qtlmap] FAILED: $workflow.runName"
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
          log.info "[eQTL-Catalogue/qtlmap] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[eQTL-Catalogue/qtlmap] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/${params.name}" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[eQTL-Catalogue/qtlmap] Pipeline Complete"

}

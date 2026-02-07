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

Channel.fromPath(params.studyFile)
    .ifEmpty { error "Cannot find studyFile file in: ${params.studyFile}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map { row -> 
        def tpm = row.containsKey('tpm_file') && row.tpm_file ? file(row.tpm_file) : null
        def is_missing = (tpm == null) ? true : false

        if (is_missing) {
            def dummy_file = file("dummy_tpm.parquet")
            dummy_file.text = ""  // Create an empty file
            tpm = dummy_file
        }

        [ row.qtl_subset, tpm, is_missing ]
    }
    .set { tpm_file_ch }

Channel.fromPath(params.rsid_map_file)
    .ifEmpty { error "Cannot find rsid file in: ${params.rsid_map_file}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.chr, file(row.rsid_file)]}
    .set { chr_rsid_map_ch }

// Batch channel
batch_ch = Channel.of(1..params.n_batches)

// Fetch the pipeline version from Git tags
def pipelineVersion = "v0.0.0" // Default version in case git command fails

// Try to fetch the version from Git
try {
    pipelineVersion = "git describe --tags".execute().text.trim()
} catch (Exception e) {
    log.warn "Could not retrieve the pipeline version from Git. Using default version $pipelineVersion."
}

// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

eQTL-Catalogue/qtlmap ${pipelineVersion}"
======================================================="""
def summary = [:]
summary['Pipeline Name']        = workflow.manifest.name
summary['Pipeline Version']     = pipelineVersion
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
include { extract_variant_info } from './modules/extract_variant_info'
include { extract_variant_info as extract_variant_info2 } from './modules/extract_variant_info'
include { prepare_molecular_traits; compress_bed; make_pca_covariates } from './modules/prepare_molecular_traits'
include { extract_samples_from_vcf } from './modules/extract_samples_from_vcf'
include { run_permutation; merge_permutation_batches; run_nominal;convert_merged_permutation_txt_to_pq} from './modules/map_qtls'
include { vcf_to_dosage } from './modules/vcf_to_dosage'
include { run_susie } from './modules/susie'
include { concatenate_pq_files; merge_cs_sumstats } from './modules/concat_pq'
include { concatenate_pq_files as concatenate_pq_files_credible_sets } from './modules/concat_pq'
include { concatenate_pq_files as concatenate_pq_files_cc } from './modules/concat_pq'
include { concatenate_pq_files as concat_pq_all } from './modules/concat_pq'
include { concatenate_pqs_wo_sorting; sort_pq_file } from './modules/concat_pq'
include { generate_sumstat_batches; convert_extracted_variant_info; convert_tpm; convert_pheno_meta} from './modules/generate_sumstat_batches'
include { extract_unique_molecular_trait_id; extract_lead_cc_signal } from './modules/extract_cc_signal'


workflow {

    // Prepare input data for QTL mapping
    if( params.vcf_set_variant_ids ){
      vcf_set_variant_ids(vcf_file_ch)
      vcf_input_ch = vcf_set_variant_ids.out
    }
    else {
      vcf_input_ch = vcf_file_ch
    }
    extract_variant_info(vcf_input_ch)
    prepare_molecular_traits(study_file_ch.join(extract_variant_info.out))
    compress_bed(prepare_molecular_traits.out.bed_file)
    extract_samples_from_vcf(vcf_input_ch.join(prepare_molecular_traits.out.sample_names))
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
      convert_merged_permutation_txt_to_pq(merge_permutation_batches.out)
    }
    //Nominal pass
    if( params.run_nominal ){
      run_nominal(batch_ch, qtlmap_input_ch)
        extract_variant_info2(extract_samples_from_vcf.out.vcf) 
        run_nominal_output= run_nominal.out.map{qtl_group,nominal_file, chromosome,start_pos,end_pos ->[qtl_group,[nominal_file,chromosome,start_pos,end_pos]]}
        nominal_qtl_subset_grouped = run_nominal_output.groupTuple(size: params.n_batches)
        all_nominal_qtl_subset_grouped = nominal_qtl_subset_grouped.map{qtl_group,nominal_run_data ->[qtl_group,nominal_run_data.flatten()]}
        convert_extracted_variant_info(extract_variant_info2.out)
        convert_tpm(tpm_file_ch)
        convert_pheno_meta(prepare_molecular_traits.out.pheno_meta)
        all_nominal_qtl_subset_info = all_nominal_qtl_subset_grouped
          .join(convert_extracted_variant_info.out) 
          .join(convert_pheno_meta.out) 
          .join(convert_tpm.out)
        nominal_qtl_subset_info_correct_format_ch = all_nominal_qtl_subset_info
            .flatMap { qtl_group, nominal_run_files_regions, extracted_variant_info, pheno_meta, tpm_file, tpm_missing ->
              nominal_run_files_regions.collate(4)
            .collect { nominal_run_file -> [qtl_group, nominal_run_file, extracted_variant_info, pheno_meta, tpm_file, tpm_missing] }}
            .map { qtl_group, nominal_run_file_region_list, extracted_variant_info, pheno_meta, tpm_file, tpm_missing ->
              def nominal_file = nominal_run_file_region_list[0]
              def chr = nominal_run_file_region_list[1]
              def start = nominal_run_file_region_list[2]
              def end = nominal_run_file_region_list[3]
              [chr, start, end, qtl_group, nominal_file, extracted_variant_info, pheno_meta, tpm_file, tpm_missing]}        
        generate_sumstat_batches_input_ch = chr_rsid_map_ch.cross(nominal_qtl_subset_info_correct_format_ch).map { rsid_data, nominal_data -> 
          def chromosome = rsid_data[0]  
          def rsid_map = rsid_data[1]    
          def start = nominal_data[1]    
          def end = nominal_data[2]      
          def qtl_group = nominal_data[3] 
          def nominal_file = nominal_data[4]  
          def extracted_variant_info = nominal_data[5]  
          def pheno_meta = nominal_data[6]  
          def tpm_file = nominal_data[7]
          def tpm_missing = nominal_data[8]
          [qtl_group, rsid_map, chromosome, start, end, nominal_file, extracted_variant_info, pheno_meta, tpm_file, tpm_missing]}
        generate_sumstat_batches(generate_sumstat_batches_input_ch)
        //Concat nominal batches
        grouped_sumstats_batches = generate_sumstat_batches.out.map { qtl_subset, file, chr, start, end -> [qtl_subset, file] }.groupTuple(size: params.n_batches)
        concat_pq_all(grouped_sumstats_batches, "all")
    }
    //Run SuSiE
    if( params.run_permutation & params.run_susie ){
      vcf_to_dosage(extract_samples_from_vcf.out.vcf)
      susie_ch = study_file_ch
        .join(merge_permutation_batches.out)
        .join(make_pca_covariates.out)
        .join(vcf_to_dosage.out)
      run_susie(susie_ch, batch_ch)
    }
    if( params.run_permutation & params.run_susie & params.run_nominal ){
      grouped_susie_cs = run_susie.out.in_cs_variant_batch.groupTuple( size: params.n_batches)
      concatenate_pq_files(grouped_susie_cs, "merged_susie")
      extract_unique_molecular_trait_id(concatenate_pq_files.out)
      sumstat_batches_crossed_uniq_mol_trait_ids = extract_unique_molecular_trait_id.out.cross(generate_sumstat_batches.out)
      extract_lead_cc_signal_ch = sumstat_batches_crossed_uniq_mol_trait_ids.map { nested_item ->
          def (unique_trait, batch) = nested_item
          def (qtl_subset, unique_molecular_trait_ids) = unique_trait
          def (batch_qtl_subset, batch_file, chr, start_pos, end_pos) = batch
          tuple(qtl_subset, unique_molecular_trait_ids, batch_file, chr, start_pos, end_pos)
      }
      extract_lead_cc_signal(extract_lead_cc_signal_ch)
      concatenate_pq_files.out
      .combine(generate_sumstat_batches.out, by: 0)  
      .map { qtl_subset, merged_parquet_file, sumstat_batch, chrom, start_pos, end_pos ->
        return tuple(qtl_subset, sumstat_batch, chrom, start_pos, end_pos, merged_parquet_file)
          }
      .set { merged_susie_sumstat_ch }
      merge_cs_sumstats(merged_susie_sumstat_ch)
      grouped_merge_cs_sumstats = merge_cs_sumstats.out.groupTuple(size: params.n_batches)
      concatenate_pq_files_credible_sets(grouped_merge_cs_sumstats, "credible_sets")
      grouped_susie_lbf = run_susie.out.lbf_variable_batch.groupTuple( size: params.n_batches)
      extract_lead_cc_signal_grouped_output = extract_lead_cc_signal.out.groupTuple(size: params.n_batches)
      concatenate_pq_files_cc(extract_lead_cc_signal_grouped_output,"cc")
      if( params.run_merge_lbf){
        concatenate_pqs_wo_sorting(grouped_susie_lbf, "lbf_variable")
        sort_pq_file(concatenate_pqs_wo_sorting.out)
      }
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

#!/usr/bin/env nextflow
/*
========================================================================================
                         kerimoff/qtlmap
========================================================================================
 kerimoff/qtlmap Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/kerimoff/qtlmap
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =======================================================
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\'
        |\\ | |__  __ /  ` /  \\ |__) |__         }  {
        | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                              `._,._,\'

     kerimoff/qtlmap v${workflow.manifest.version}
    =======================================================

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf\
     -profile tartu_hpc\
     --expression_matrix testdata/GEUVADIS_cqn.tsv\
     --phenotype_metadata testdata/GEUVADIS_phenotype_metadata.tsv\
     --sample_metadata testdata/GEUVADIS_sample_metadata.tsv\
     --genotype_vcf testdata/GEUVADIS_genotypes.vcf.gz\
     --is_imputed FALSE

    Mandatory arguments:
      --expression_matrix           Path to phenotype matrix data (.tsv - normalized and quality controlled)
      --phenotype_metadata          Path to phenotype metadata (.tsv)
      --sample_metadata             Path to sample metadata (.tsv)
      --genotype_vcf                Path to genotype (VCF) file 

    Options:
      --cis_window                  The window where to search for associated variants around the phenotype (default: 1000000)
      --mincisvariant               Minimum variants needed to be found in cis_window (default: 56)   
      --n_batches                   Number of parallel batches to run QTL Mapping per sample (default: 400)
      --is_imputed                  Does the genotype VCF file contain R2 value in the INFO field? (default: true)
      --run_permutation             Calculate permuation p-values for each phenotype group (group_id in the phenotype metadata file) (default: false)

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
         .ifEmpty { exit 1, "Cannot find the expression matrix file: ${params.expression_matrix}\n" }
         .set { expression_matrix_create_QTLTools_input}
     Channel
         .fromPath( params.sample_metadata )
         .ifEmpty { exit 1, "Cannot find the sample metadata file: ${params.sample_metadata}\n" }
         .set { sample_metadata_create_QTLTools_input}   
    Channel
         .fromPath( params.phenotype_metadata )
         .ifEmpty { exit 1, "Cannot find the phenotype metadata file: ${params.phenotype_metadata}\n" }
         .set { phenotype_metadata_create_QTLTools_input}
    Channel
         .fromPath( params.genotype_vcf )
         .ifEmpty { exit 1, "Cannot find the genotype vcf file: ${params.genotype_vcf}\n" }
         .into { genotype_vcf_extract_variant_info; genotype_vcf_extract_samples }
         
 }


// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

kerimoff/qtlmap v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']        = 'kerimoff/qtlmap'
summary['Pipeline Version']     = workflow.manifest.version
summary['Run Name']             = custom_runName ?: workflow.runName
summary['Expression Matrix']    = params.expression_matrix
summary['Phenotype Metadata']   = params.phenotype_metadata
summary['Sample Metadata']      = params.sample_metadata
summary['Genotype VCF file']    = params.genotype_vcf
summary['Cis window']           = params.cis_window
summary['Minimum Cis variants'] = params.mincisvariant
summary['Is imputed']           = params.is_imputed
summary['Permutation run']      = params.run_permutation
summary['# of batches']         = params.n_batches
summary['# of phenotype pcs']   = params.n_pheno_pcs
summary['# of genotype pcs']    = params.n_geno_pcs
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


/*
 * STEP 0 - Extract variant information from VCF
 */
process extract_all_variant_info {
    tag "${vcf.simpleName}"
    // publishDir "${params.outdir}/final", mode: 'copy'

    input:
    file vcf from genotype_vcf_extract_variant_info
    
    output:
    file "${vcf.simpleName}.variant_information.txt.gz" into variant_info_create_QTLTools_input

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

/*
 * STEP 1 - Generate QTLTools input files
 */
process create_QTLTools_input {
    tag "${expression_matrix.baseName}"
    publishDir "${params.outdir}/", mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".phenoPCA.tsv") > 0) "PCA/$filename" else null
    }

    input:
    file expression_matrix from expression_matrix_create_QTLTools_input.collect()
    file sample_metadata from sample_metadata_create_QTLTools_input.collect()
    file phenotype_metadata from phenotype_metadata_create_QTLTools_input.collect()
    file study_variant_info from variant_info_create_QTLTools_input.collect()
    
    output: // set can be used to pass condition val and file as tuple to the channel 
    file "*.bed" into condition_beds mode flatten
    file "*.sample_names.txt" into condition_samplenames mode flatten
    file "*.phenoPCA.tsv" into pheno_PCAs mode flatten

    script:
    """
    Rscript $baseDir/bin/group_by_qtlgroup.R \\
        -p "$phenotype_metadata" \\
        -s "$sample_metadata" \\
        -e "$expression_matrix" \\
        -v "$study_variant_info" \\
        -o "." \\
        -c ${params.cis_window} \\
        -m ${params.mincisvariant}
    """
}

/*
 * STEP 2 - Compres and index input bed file
 */ 
process compress_bed {
    tag "${bed_file.baseName}"
    // publishDir "${params.outdir}/compressed_bed", mode: 'copy'

    input:
    file bed_file from condition_beds

    output:
    set val(bed_file.simpleName), file("${bed_file}.gz") into compressed_beds //compressed_beds_perform_pca, compressed_beds_run_nominal, compressed_beds_run_permutation
    set val(bed_file.simpleName), file("${bed_file}.gz.tbi") into compressed_bed_indexes //compressed_bed_indexes_perform_pca, compressed_bed_indexes_run_nominal, compressed_bed_indexes_run_permutation

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
    // publishDir "${params.outdir}/vcf", mode: 'copy'

    input:
    file sample_names from condition_samplenames
    file genotype_vcf from genotype_vcf_extract_samples.collect()

    output:
    set val(sample_names.simpleName), file("${sample_names.simpleName}.vcf.gz") into vcfs_extract_variant_info, vcfs, vcfs_perform_pca 
    set val(sample_names.simpleName), file("${sample_names.simpleName}.vcf.gz.csi") into vcf_indexes 

    script:
    """
    bcftools view -S $sample_names $genotype_vcf -Oz -o ${sample_names.simpleName}.vcf.gz
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
    set condition, file(vcf) from vcfs_extract_variant_info
    
    output:
    file "${condition}.variant_information.txt.gz"

    script:
    if (params.is_imputed) {
        """
        set +o pipefail; bcftools +fill-tags $vcf | bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%TYPE\\t%AC\\t%AN\\t%MAF\\t%R2\\n' | gzip > ${condition}.variant_information.txt.gz
        """
    } else {
        """
        set +o pipefail; bcftools +fill-tags $vcf | bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%TYPE\\t%AC\\t%AN\\t%MAF\\tNA\\n' | gzip > ${condition}.variant_information.txt.gz
        """
    }
}

compressed_beds.join(compressed_bed_indexes).join(vcfs).join(vcf_indexes).into{ tuple_run_nominal; tuple_run_permutation }
pheno_PCAs.map { file -> [ file.simpleName, file] }.join(vcfs_perform_pca).set{ tuple_perform_pca }

/*
 * STEP 5 - Perform PCA on the genotype and phenotype data
 */
process perform_pca {
    tag "${condition}"
    publishDir "${params.outdir}/PCA", mode: 'copy'

    input:
    set condition, file(phenotype_pca), file(vcf) from tuple_perform_pca

    output:
    file "${condition}.geno.pca*"
    set val(condition), file("${condition}.covariates.txt") into covariates_run_nominal, covariates_run_permutation

    script:
    """
    plink2 --vcf $vcf --vcf-half-call h --indep-pairwise 50000 200 0.05 --out ${vcf.simpleName}_pruned_variants --threads ${task.cpus} --memory 12000
    plink2 --vcf $vcf --vcf-half-call h --extract ${vcf.simpleName}_pruned_variants.prune.in --make-bed --out ${vcf.simpleName}_pruned
    plink2 -bfile ${vcf.simpleName}_pruned --pca ${params.n_geno_pcs} header tabs
    cat plink.eigenvec \\
        | sed '1s/IID/genotype_id/' \\
        | sed '1s/PC/geno_PC/g' \\
        | csvtk cut -t -f -"FID" \\
        | csvtk transpose -t > ${condition}.geno.pca
    head -n ${params.n_pheno_pcs + 1} $phenotype_pca > ${condition}.covariates.txt    
    set +o pipefail; tail -n+2 ${condition}.geno.pca | head -n ${params.n_geno_pcs} >> ${condition}.covariates.txt
    """
}

tuple_run_nominal
        .join(covariates_run_nominal)
        .set{tuple_run_nominal}

tuple_run_permutation
        .join(covariates_run_permutation)
        .set{tuple_run_permutation}

/*
 * STEP 6 - Run QTLtools in permutation mode
 */
process run_permutation {
    tag "${condition} - ${batch_index}/${params.n_batches}"
    // publishDir "${params.outdir}/temp_batches", mode: 'copy'
    
    when:
    params.run_permutation

    input:
    each batch_index from 1..params.n_batches
    set condition, file(bed), file(bed_index), file(vcf), file(vcf_index), file(covariate) from tuple_run_permutation

    output:
    set val(condition), file("${condition}.permutation.batch.${batch_index}.${params.n_batches}.txt") into batch_files_merge_permutation_batches

    script:
    """
    QTLtools cis --vcf $vcf --bed $bed --cov $covariate --chunk $batch_index ${params.n_batches} --out ${condition}.permutation.batch.${batch_index}.${params.n_batches}.txt --window ${params.cis_window} --permute 10000 --grp-best
    """
}

/*
 * STEP 7 - Merge permutation batches from QTLtools
 */
process merge_permutation_batches {
    tag "${condition}"
    publishDir "${params.outdir}/final", mode: 'copy'
    
    when:
    params.run_permutation

    input:
    set condition, batch_file_names from batch_files_merge_permutation_batches.groupTuple(size: params.n_batches, sort: true)  

    output:
    file "${condition}.permuted.txt.gz"

    script:
    """
    cat ${batch_file_names.join(' ')} | bgzip > ${condition}.permuted.txt.gz
    """
}


/*
 * STEP 8 - Run QTLtools in nominal mode
 */
process run_nominal {
    tag "${condition} - ${batch_index}/${params.n_batches}"
    // publishDir "${params.outdir}/temp_batches", mode: 'copy'
    
    input:
    each batch_index from 1..params.n_batches
    set condition, file(bed), file(bed_index), file(vcf), file(vcf_index), file(covariate) from tuple_run_nominal

    output:
    set val(condition), file("${condition}.nominal.batch.${batch_index}.${params.n_batches}.txt") into batch_files_merge_nominal_batches

    script:
    """
	QTLtools cis --vcf $vcf --bed $bed --cov $covariate --chunk $batch_index ${params.n_batches} --out ${condition}.nominal.batch.${batch_index}.${params.n_batches}.txt --window ${params.cis_window} --nominal 1
    """
}

/*
 * STEP 9 - Merge nominal batches from QTLtools
 */
process merge_nominal_batches {
    tag "${condition}"
    // publishDir "${params.outdir}/Nominal_merged", mode: 'copy'

    input:
    set condition, batch_file_names from batch_files_merge_nominal_batches.groupTuple(size: params.n_batches, sort: true)  

    output:
    set val(condition), file("${condition}.nominal.txt.gz") into nominal_merged_files_replace_space_tabs

    script:
    """
    cat ${batch_file_names.join(' ')} | bgzip > ${condition}.nominal.txt.gz
    """
}

/*
 * STEP 10 - Replace spaces with tabs 
 */
process replace_space_tabs {
    tag "${condition}"
    // publishDir "${params.outdir}/Nominal_merged", mode: 'copy'
	
    input:
    set condition, file(nominal_merged) from nominal_merged_files_replace_space_tabs

    output:
    set val(condition), file("${condition}.nominal.tab.txt.gz") into nominal_merged_tab_sort_qtltools_output
    
    script:
    """
    gzip -dc $nominal_merged | awk -v OFS='\\t' '{{\$1=\$1; print \$0}}' | gzip > ${condition}.nominal.tab.txt.gz
    """
}

/*
 * STEP 11 - Add SNP coordinates to QTLTools output file
 */
process sort_qtltools_output {
    tag "${condition}"
    publishDir "${params.outdir}/final", mode: 'copy'

    input:
    set condition, file(nominal_merged) from nominal_merged_tab_sort_qtltools_output

    output:
    set val(condition), file("${condition}.nominal.sorted.txt.gz") into sorted_merged_nominal_index_qtltools_output

    script:
    """
    gzip -dc $nominal_merged | LANG=C sort -k9,9 -k10,10n -k11,11n -S11G --parallel=8 | bgzip > ${condition}.nominal.sorted.txt.gz
    """
}

/*
 * STEP 12 - Tabix-index QTLtools output files
 */
process index_qtltools_output {
    tag "${condition}"
    publishDir "${params.outdir}/final", mode: 'copy'

    input:
    set condition, file(sorted_merged_nominal) from sorted_merged_nominal_index_qtltools_output

    output:
    file "${condition}.nominal.sorted.txt.gz.tbi"

    script:
    """
    tabix -s9 -b10 -e11 -f $sorted_merged_nominal
    """
}

// TODO: try to use each input repeater for permutations

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[kerimoff/qtlmap] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[kerimoff/qtlmap] FAILED: $workflow.runName"
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
          log.info "[kerimoff/qtlmap] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[kerimoff/qtlmap] Sent summary e-mail to $params.email (mail)"
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

    log.info "[kerimoff/qtlmap] Pipeline Complete"

}

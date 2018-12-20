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
fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}
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
ch_multiqc_config = Channel.fromPath(params.multiqc_config)
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")

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
         .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
         .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
         .into { read_files_fastqc; read_files_trimming }
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
summary['Pipeline Name']  = 'nf-core/qtlmap'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Reads']        = params.reads
summary['Fasta Ref']    = params.fasta
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
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
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    // TODO nf-core: Get all tools to print their version number here
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}



/*
 * STEP 1 - Generate QTLTools input files
 */
process create_QTLTools_input {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    Rscript array_data_to_QTLtools_input.R\
     -g "../metadata/gene_metadata/HumanHT-12_V4_gene_metadata.txt.gz"\
     -s "../metadata/cleaned/Fairfax_2014.tsv"\
     -e "/gpfs/hpc/home/a72094/projects/RNAseq_pipeline/results/expression_matrices/HumanHT-12_V4/Fairfax_2014.tsv.gz"\
     -v "/gpfs/hpchome/a72094/hpc/datasets/controlled_access/Fairfax_2014/genotypes/Michigan_GRCh37_1KGPhase3_061118/GRCh38/Fairfax_2014_GRCh38.variant_information.txt.gz"\
     --qtlutils "../../eQTLUtils"\
     -o "../processed/Fairfax_2014/qtltools/input/array/"\
     -c 1000001\
     -m 6
    """
}



/*
 * STEP 2 - Compres and index input bed file
 */
process compress_bed {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'
    // 	threads: 1
	// resources:
	// 	mem = 100
    input:
        // bed = "processed/{study}/qtltools/input/{annot_type}/{condition}.bed"

    file multiqc_config from ch_multiqc_config
    // TODO nf-core: Add in log files from your new processes for MultiQC to find!
    file ('fastqc/*') from fastqc_results.collect().ifEmpty([])
    file ('software_versions/*') from software_versions_yaml
    file workflow_summary from create_workflow_summary(summary)

    output:
    	// 	bed = protected("processed/{study}/qtltools/input/{annot_type}/{condition}.bed.gz"),
		// bed_index = protected("processed/{study}/qtltools/input/{annot_type}/{condition}.bed.gz.tbi")
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    // TODO nf-core: Specify which MultiQC modules to use with -m for a faster run time
    """
    module load samtools-1.6
    bgzip {input.bed} && tabix -p bed {output.bed}
    """
}



/*
 * STEP 3 - Extract samples from vcf
 */
process extract_samples {
    publishDir "${params.outdir}/Documentation", mode: 'copy'
	// threads: 1
	// resources:
	// 	mem = 100
    input:
        // samples = "processed/{study}/qtltools/input/{annot_type}/{condition}.sample_names.txt"
    file output_docs from ch_output_docs

    output:
        // vcf = "processed/{study}/qtltools/input/{annot_type}/vcf/{condition}.vcf.gz",
		// vcf_index = "processed/{study}/qtltools/input/{annot_type}/vcf/{condition}.vcf.gz.csi"
    file "results_description.html"

    script:
    """
    module load bcftools-1.8
    bcftools view -S {input.samples} {config[vcf_file]} -Oz -o {output.vcf}
    bcftools index {output.vcf}
    """
}

/*
 * STEP 4 - Extract variant information from VCF
 */
process extract_samples {
    publishDir "${params.outdir}/Documentation", mode: 'copy'
	// threads: 1
	// resources:
	// 	mem = 100

    input:
        // vcf = "processed/{study}/qtltools/input/{annot_type}/vcf/{condition}.vcf.gz",
    file output_docs from ch_output_docs

    output:
        // var_info = "processed/{study}/qtltools/output/{annot_type}/final/{condition}.variant_information.txt.gz"
    file "results_description.html"

    script:
    """
    module load bcftools-1.8
    bcftools view -S {input.samples} {config[vcf_file]} -Oz -o {output.vcf}
    bcftools index {output.vcf}
    """
}


/*
 * STEP 5 - Perform PCA on the genotype and phenotype data
 */
process perform_pca {
    publishDir "${params.outdir}/Documentation", mode: 'copy'
	// threads: 1
	// resources:
	// 	mem = 2000

    input:
        // bed = "processed/{study}/qtltools/input/{annot_type}/{condition}.bed.gz",
		// bed_index = "processed/{study}/qtltools/input/{annot_type}/{condition}.bed.gz.tbi",
		// vcf = "processed/{study}/qtltools/input/{annot_type}/vcf/{condition}.vcf.gz",
		// vcf_index = "processed/{study}/qtltools/input/{annot_type}/vcf/{condition}.vcf.gz.csi"
    file output_docs from ch_output_docs

    output:
        // covariates = "processed/{study}/qtltools/input/{annot_type}/{condition}.covariates.txt"
    file "results_description.html"

        // params:
        // 	pheno_pca = "processed/{study}/qtltools/input/{annot_type}/{condition}.pheno",
        // 	geno_pca = "processed/{study}/qtltools/input/{annot_type}/vcf/{condition}.geno"

    script:
    """
		QTLtools pca --bed {input.bed} --center --scale --out {params.pheno_pca}
		QTLtools pca --vcf {input.vcf} --maf 0.05 --center --scale --distance 50000 --out {params.geno_pca}
		head -n 7 {params.pheno_pca}.pca > {output.covariates}
		set +o pipefail; tail -n+2 {params.geno_pca}.pca | head -n 3 >> {output.covariates}
    """
}


/*
 * STEP 6 - Run QTLtools in permutation mode
 */
process run_permutation {
    publishDir "${params.outdir}/Documentation", mode: 'copy'
	// threads: 1
	// resources:
	// 	mem = 5000

    input:
		// bed = "processed/{study}/qtltools/input/{annot_type}/{condition}.bed.gz",
		// bed_index = "processed/{study}/qtltools/input/{annot_type}/{condition}.bed.gz.tbi",
		// covariates = "processed/{study}/qtltools/input/{annot_type}/{condition}.covariates.txt",
		// vcf = "processed/{study}/qtltools/input/{annot_type}/vcf/{condition}.vcf.gz"
    file output_docs from ch_output_docs

    output:
        // temp("processed/{study}/qtltools/output/{annot_type}/batches/{condition}.permutation.batch.{batch}.{n_batches}.txt")
    file "results_description.html"

	// params:
	// 	chunk = "{batch} {n_batches}"

    script:
    """
    QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.covariates} --chunk {params.chunk} --out {output} --window {config[cis_window]} --permute 10000 --grp-best
    """
}


/*
 * STEP 7 - Merge all batches from QTLtools
 */
process merge_permutation_batches {
    publishDir "${params.outdir}/Documentation", mode: 'copy'
	// threads: 1
	// resources:
	// 	mem = 100
    input:
        // expand("processed/{{study}}/qtltools/output/{{annot_type}}/batches/{{condition}}.permutation.batch.{batch}.{n_batches}.txt", 
		// 	batch=[i for i in range(1, config["n_batches"] + 1)],
		// 	n_batches = config["n_batches"])
    file output_docs from ch_output_docs

    output:
        // "processed/{study}/qtltools/output/{annot_type}/{condition}.permuted.txt.gz"
    file "results_description.html"

    script:
    """
		module load samtools-1.6
		cat {input} | bgzip > {output}
    """
}

/*
 * STEP 8 - Run QTLtools in nominal mode
 */
process run_nominal {
    publishDir "${params.outdir}/Documentation", mode: 'copy'
	// threads: 1
	// resources:
	// 	mem = 5000
    input:
		// bed = "processed/{study}/qtltools/input/{annot_type}/{condition}.bed.gz",
		// bed_index = "processed/{study}/qtltools/input/{annot_type}/{condition}.bed.gz.tbi",
		// covariates = "processed/{study}/qtltools/input/{annot_type}/{condition}.covariates.txt",
		// vcf = "processed/{study}/qtltools/input/{annot_type}/vcf/{condition}.vcf.gz"
    file output_docs from ch_output_docs

    output:
        // temp("processed/{study}/qtltools/output/{annot_type}/{condition}.nominal.txt.gz")
    file "results_description.html"

	// params:
	// 	chunk = "{batch} {n_batches}"
    script:
    """
		module load samtools-1.6
		cat {input} | bgzip > {output}
    """
}

/*
 * STEP 9 - Merge all batches from QTLtools
 */
process merge_nominal_batches {
    publishDir "${params.outdir}/Documentation", mode: 'copy'
	// threads: 1
	// resources:
	// 	mem = 100
    input:
		// 		expand("processed/{{study}}/qtltools/output/{{annot_type}}/nominal_batches/{{condition}}.nominal.batch.{batch}.{n_batches}.txt", 
			// batch=[i for i in range(1, config["n_batches"] + 1)],
			// n_batches = config["n_batches"])
    file output_docs from ch_output_docs

    output:
        // temp("processed/{study}/qtltools/output/{annot_type}/{condition}.nominal.txt.gz")
    file "results_description.html"

    script:
    """
		module load samtools-1.6
		cat {input} | bgzip > {output}
    """
}

/*
 * STEP 10 - Replace tabs 
 */
process replace_space_tabs {
    publishDir "${params.outdir}/Documentation", mode: 'copy'
	// threads: 2
	// resources:
	// 	mem = 1000
    input:
		// "processed/{study}/qtltools/output/{annot_type}/{condition}.nominal.txt.gz"
    file output_docs from ch_output_docs

    output:
        // "processed/{study}/qtltools/output/{annot_type}/tab/{condition}.nominal.txt.gz"
    file "results_description.html"

    script:
    """
    gzip -dc {input} | awk -v OFS='\\t' '{{$1=$1; print $0}}' | gzip > {output}
    """
}


/*
 * STEP 11 - Add SNP coordinates to QTLTools output file
 */
process sort_qtltools_output {
    publishDir "${params.outdir}/Documentation", mode: 'copy'
	// threads: 10
	// resources:
	// 	mem = 12000
    input:
		// "processed/{study}/qtltools/output/{annot_type}/tab/{condition}.nominal.txt.gz"
    file output_docs from ch_output_docs

    output:
        // protected("processed/{study}/qtltools/output/{annot_type}/final/{condition}.nominal.sorted.txt.gz")
    file "results_description.html"

	// params:
	// 	chunk = "{batch} {n_batches}"
    script:
    """
    	module load samtools-1.6
		gzip -dc {input} | LANG=C sort -k9,9 -k10,10n -k11,11n -S11G --parallel=8 | bgzip > {output}
    """
}

/*
 * STEP 12 - Tabix-index QTLtools output files
 */
process index_qtltools_output {
    publishDir "${params.outdir}/Documentation", mode: 'copy'
	// threads: 1
	// resources:
	// 	mem = 1000
    input:
		// "processed/{study}/qtltools/output/{annot_type}/final/{condition}.nominal.sorted.txt.gz"
    file output_docs from ch_output_docs

    output:
        // "processed/{study}/qtltools/output/{annot_type}/final/{condition}.nominal.sorted.txt.gz.tbi"
    file "results_description.html"

	// params:
	// 	chunk = "{batch} {n_batches}"
    script:
    """
		module load samtools-1.6
		tabix -s9 -b10 -e11 -f {input}
    """
}





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

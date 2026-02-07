process run_susie{
    container = 'quay.io/kfkf33/susier:v24.01.1'
    publishDir "${params.outdir}/susie_batches/${qtl_subset}/cs/", mode: 'copy', pattern: "${qtl_subset}.${batch_index}_${params.n_batches}.parquet"
    publishDir "${params.outdir}/susie_batches/${qtl_subset}/lbf/", mode: 'copy', pattern: "${qtl_subset}.${batch_index}_${params.n_batches}.lbf_variable.parquet"
    publishDir "${params.outdir}/susie_batches/${qtl_subset}/full/", mode: 'copy', pattern: "${qtl_subset}.${batch_index}_${params.n_batches}.full_susie.parquet"

    input:
    tuple val(qtl_subset), file(expression_matrix), file(phenotype_meta), file(sample_meta), file(phenotype_list), file(covariates), file(genotype_matrix), file(genotype_matrix_index)
    each batch_index

    output:
    tuple val(qtl_subset), path("${qtl_subset}.${batch_index}_${params.n_batches}.parquet"), emit: in_cs_variant_batch 
    tuple val(qtl_subset), path("${qtl_subset}.${batch_index}_${params.n_batches}.lbf_variable.parquet"), emit: lbf_variable_batch
    tuple val(qtl_subset), path("${qtl_subset}.${batch_index}_${params.n_batches}.full_susie.parquet"), emit: full_susie_batch

    script:
    """
    Rscript $baseDir/bin/run_susie.R --expression_matrix ${expression_matrix}\
     --phenotype_meta ${phenotype_meta}\
     --sample_meta ${sample_meta}\
     --phenotype_list ${phenotype_list}\
     --covariates ${covariates}\
     --genotype_matrix ${genotype_matrix}\
     --chunk '${batch_index} ${params.n_batches}'\
     --cisdistance ${params.cis_window}\
     --out_prefix '${qtl_subset}.${batch_index}_${params.n_batches}'\
     --eqtlutils null\
     --write_full_susie ${params.write_full_susie}\
     --finemap_by_group_id ${params.finemap_by_group_id}
    """
}
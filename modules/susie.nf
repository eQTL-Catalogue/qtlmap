process run_susie{
    container = 'quay.io/kfkf33/susier'
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
     --write_full_susie ${params.write_full_susie}
    """
}


process concatenate_pqs_wo_sorting {
    tag "${qtl_subset}"
    container = 'quay.io/kfkf33/duckdb_env'


    input:
    tuple val(qtl_subset), val(files)
    val(output_postfix)

    output:
    tuple val(qtl_subset), path("${qtl_subset}_${output_postfix}.parquet")

    script:
    """
    concatenate_pqs_without_sorting.py -f ${files.join(' ')} -o ${qtl_subset}_${output_postfix}.parquet -m ${task.memory.toMega() / 1024}
    """
}

process sort_pq_file {
    tag "${qtl_subset}"
    container = 'quay.io/kfkf33/duckdb_env'
    publishDir "${params.outdir}/susie/", mode: 'copy', pattern: "*merged.parquet"

    input:
    tuple val(qtl_subset), path(pq_file)

    output:
    tuple val(qtl_subset), path("*merged.parquet")

    script:
    """
    sort_concatenated_pq.py -i ${pq_file} -m ${task.memory.toMega() / 1024} -n ${pq_file.simpleName}
    """
}

process concatenate_pq_files {
    tag "${qtl_subset}"
    container = 'quay.io/kfkf33/duckdb_env'
    publishDir "${params.outdir}/susie/", mode: 'copy', pattern: "*credible_sets.parquet"


    input:
    tuple val(qtl_subset), val(files)
    val(output_postfix)

    output:
    tuple val(qtl_subset), path("${qtl_subset}_${output_postfix}.parquet")

    script:
    """
    concatenate_pq_files.py -f ${files.join(' ')} -o ${qtl_subset}_${output_postfix}.parquet
    """
}

process merge_cs_sumstats{
    tag "${qtl_subset}"
    container = 'quay.io/kfkf33/duckdb_env'

    input:
    tuple val(qtl_subset), path(sumstat_batch), val(chrom), val(start_pos), val(end_pos),path(merged_susie_file)

    output:
    tuple val(qtl_subset), path("merged_cs_sumstat_${qtl_subset}_${chrom}_${start_pos}_${end_pos}.parquet")

    script:
    """
    merge_cs_sumstats.py \
        --merged_susie_file ${merged_susie_file} \
        --sumstat_btach_file ${sumstat_batch} \
        --chrom ${chrom} \
        --start_pos ${start_pos} \
        --end_pos ${end_pos} \
        --output_file merged_cs_sumstat_${qtl_subset}_${chrom}_${start_pos}_${end_pos}.parquet
    """
}
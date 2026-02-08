process concatenate_pqs_wo_sorting {
    tag "${qtl_subset}"
    container = 'quay.io/kfkf33/duckdb_env:v24.01.1'


    input:
    tuple val(qtl_subset), val(files)
    val(output_postfix)

    output:
    tuple val(qtl_subset), val(output_postfix), path("${qtl_subset}_${output_postfix}.parquet")

    script:
    """
    concatenate_pqs_without_sorting.py -f ${files.join(' ')} -o ${qtl_subset}_${output_postfix}.parquet -m ${task.memory.toMega() / 1024}
    """
}

process sort_pq_file {
    tag "${qtl_subset}"
    container = 'quay.io/kfkf33/duckdb_env:v24.01.1'
    publishDir "${params.outdir}/susie/${qtl_subset}/", mode: 'copy', pattern: "${qtl_subset}.${output_postfix}.parquet"

    input:
    tuple val(qtl_subset),val(output_postfix), path(pq_file)

    output:
    tuple val(qtl_subset), path("${qtl_subset}.${output_postfix}.parquet")

    script:
    """
    sort_concatenated_pq.py -i ${pq_file} -m ${task.memory.toMega() / 1024} -n ${qtl_subset}.${output_postfix}.parquet
    """
}

process concatenate_pq_files {
    tag "${qtl_subset}"
    container = 'quay.io/kfkf33/duckdb_env:v24.01.1'
    publishDir "${params.outdir}/susie/${qtl_subset}/", mode: 'copy', pattern: "*credible_sets.parquet"
    publishDir "${params.outdir}/sumstats/${qtl_subset}/", mode: 'copy', pattern: "*cc.parquet"
    publishDir "${params.outdir}/sumstats/${qtl_subset}/", mode: 'copy', pattern: "*all.parquet"

    input:
    tuple val(qtl_subset), val(files)
    val(output_postfix)

    output:
    tuple val(qtl_subset), path("${qtl_subset}.${output_postfix}.parquet")

    script:
    """
    concatenate_pq_files.py -f ${files.join(' ')} -o ${qtl_subset}.${output_postfix}.parquet -m ${task.memory.toMega() / 1024}
    """
}

process merge_cs_sumstats{
    tag "${qtl_subset}"
    container = 'quay.io/kfkf33/duckdb_env:v24.01.1'

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
        --memory_limit ${task.memory.toMega() / 1024} \
        --output_file merged_cs_sumstat_${qtl_subset}_${chrom}_${start_pos}_${end_pos}.parquet
    """
}
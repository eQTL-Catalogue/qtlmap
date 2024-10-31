process generate_sumstat_batches {
    tag "${qtl_subset}"
    publishDir "${params.outdir}/nominal_sumstats_batches/${qtl_subset}", mode: 'copy', pattern: "${qtl_subset}_chr*.parquet"
    container = 'quay.io/kfkf33/duckdb_env'

    input:
    tuple val(qtl_subset), path(rsid_map), val(chr),val(start_pos),val(end_pos), path(summ_stats_batch),path(var_info), path(phenotype_metadata), path(median_tpm)

    output:
    tuple val(qtl_subset), path("${qtl_subset}_chr_${region}.parquet"), val(chr), val(start_pos), val(end_pos)

    script:
    region = chr + "_" + start_pos + "_" + end_pos

    """
    $baseDir/bin/generate_sumstats.py \
        -s $summ_stats_batch \
        -v $var_info \
        -r $rsid_map \
        -p $phenotype_metadata \
        -m $median_tpm \
        -a $start_pos \
        -b $end_pos \
        -o ${qtl_subset}_chr_${region}.parquet
    """
}

process convert_extracted_variant_info {
    tag "${qtl_subset}"
    container = 'quay.io/kfkf33/duckdb_env'

    input:
    tuple val(qtl_subset), path(variant_info)

    output:
    tuple val(qtl_subset), path("${qtl_subset}_${variant_info.simpleName}.parquet")

    script:
    """
    $baseDir/bin/convert_txt_to_pq.py \
        -i $variant_info \
        -c chromosome,position,variant,ref,alt,type,ac,an,r2 \
        -s '{"chromosome":"VARCHAR","position":"INTEGER","variant":"VARCHAR","ref":"VARCHAR","alt":"VARCHAR","type":"VARCHAR","ac":"INTEGER","an":"INTEGER","maf":"DOUBLE","r2":"VARCHAR"}' \
        -o ${qtl_subset}_${variant_info.simpleName}.parquet
    """
}

process convert_tmp {
    tag "${qtl_subset}"
    container = 'quay.io/kfkf33/duckdb_env'

    input:
    tuple val(qtl_subset), path(tmp_file)

    output:
    tuple val(qtl_subset), path("${qtl_subset}_${tmp_file.simpleName}.parquet")

    script:
    """
     $baseDir/bin/convert_txt_to_pq.py \
        -i $tmp_file \
        -c phenotype_id,median_tpm \
        -o ${qtl_subset}_${tmp_file.simpleName}.parquet
    """
}

process convert_pheno_meta {
    tag "${qtl_subset}"
    container = 'quay.io/kfkf33/duckdb_env'

    input:
    tuple val(qtl_subset), path(phenotype_metadata)

    output:
    tuple val(qtl_subset), path("${qtl_subset}_${phenotype_metadata.simpleName}.parquet")

    script:
    """
    $baseDir/bin/convert_txt_to_pq.py \
        -i $phenotype_metadata \
        -c phenotype_id,group_id,gene_id \
        -o ${qtl_subset}_${phenotype_metadata.simpleName}.parquet
    """
}
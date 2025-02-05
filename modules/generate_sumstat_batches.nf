process generate_sumstat_batches {
    tag "${qtl_subset}"
    publishDir "${params.outdir}/sumstats/${qtl_subset}/all/", mode: 'copy', pattern: "${qtl_subset}_chr*.parquet"
    container = 'quay.io/kfkf33/duckdb_env:v24.01.1'

    input:
    tuple val(qtl_subset), path(rsid_map), val(chr),val(start_pos),val(end_pos), path(summ_stats_batch),path(var_info), path(phenotype_metadata), path(median_tpm), val(tpm_missing)

    output:
    tuple val(qtl_subset), path("${qtl_subset}_chr_${region}.parquet"), val(chr), val(start_pos), val(end_pos)

    script:
    region = chr + "_" + start_pos + "_" + end_pos
    missing_tpm_arg = tpm_missing ? "-t 1" : "-t 0"

    """
    $baseDir/bin/generate_sumstats.py \
        -s $summ_stats_batch \
        -v $var_info \
        -r $rsid_map \
        -p $phenotype_metadata \
        -o ${qtl_subset}_chr_${region}.parquet \
        -a $start_pos \
        -b $end_pos \
        -m $median_tpm \
        $missing_tpm_arg
    """
}

process convert_extracted_variant_info {
    tag "${qtl_subset}"
    container = 'quay.io/kfkf33/duckdb_env:v24.01.1'

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

process convert_tpm {
    tag "${qtl_subset}"
    container = 'quay.io/kfkf33/duckdb_env:v24.01.1'

    input:
    tuple val(qtl_subset), path(tpm_file), val(tpm_missing)

    output:
    tuple val(qtl_subset), path("${qtl_subset}_${tpm_file.simpleName}.parquet"), val(tpm_missing)

    script:

    if (!tpm_missing) {
        """
        $baseDir/bin/convert_txt_to_pq.py \
            -o ${qtl_subset}_${tpm_file.simpleName}.parquet \
            -c phenotype_id,median_tpm \
            -i ${tpm_file}
        """
    } else {
        """
        touch ${qtl_subset}_${tpm_file.simpleName}.parquet  # Create an dummy parquet file
        """
    }
}

process convert_pheno_meta {
    tag "${qtl_subset}"
    container = 'quay.io/kfkf33/duckdb_env:v24.01.1'

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
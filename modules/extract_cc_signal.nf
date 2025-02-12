process extract_unique_molecular_trait_id {
    tag "${qtl_subset}"
    container = 'quay.io/kfkf33/duckdb_env:v24.01.1'

    input:
    tuple val(qtl_subset), path(concatenated_susie_output)


    output:
    tuple val(qtl_subset), path("${qtl_subset}_extracted_unique_molecular_trait_ids.parquet")

    script:
    """
    extract_unique_molecular_trait_ids.py -f $concatenated_susie_output -o ${qtl_subset}_extracted_unique_molecular_trait_ids.parquet -m ${task.memory.toMega() / 1024}

    """
}

process extract_lead_cc_signal {
    tag "${qtl_subset}"
    container = 'quay.io/kfkf33/duckdb_env:v24.01.1'
    //publishDir "${params.outdir}/lead_cc_signal/${qtl_subset}", mode: 'copy', pattern: "${qtl_subset}_chr${region}_cc.parquet"


    input:
    tuple val(qtl_subset), path(unique_molecular_trait_ids), path(qtl_ss) , val(chr), val(start_pos), val(end_pos)

    output:
    tuple val(qtl_subset), file("${qtl_subset}_chr${region}_cc.parquet")

    script:
    region = chr + "_" + start_pos + "_" + end_pos

    """
    extract_lead_cc_signal.py -s $qtl_ss -u $unique_molecular_trait_ids -o ${qtl_subset}_chr${region}_cc.parquet -m ${task.memory.toMega() / 1024}

    """
} 
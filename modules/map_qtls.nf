/*
 * STEP 6 - Run QTLtools in permutation mode
 */
process run_permutation {
    tag "${qtl_subset} - ${batch_index}/${params.n_batches}"
    container = 'quay.io/eqtlcatalogue/qtltools:v22.03.1'

    input:
    each batch_index
    tuple val(qtl_subset), file(bed), file(bed_index), file(fastqtl_bed), file(fastqtl_bed_index), file(vcf), file(vcf_index), file(covariate)

    output:
    tuple val(qtl_subset), file("${qtl_subset}.permutation.batch.${batch_index}.${params.n_batches}.txt")

    script:
    """
    QTLtools cis --vcf $vcf --bed $bed --cov $covariate --chunk $batch_index ${params.n_batches} --out ${qtl_subset}.permutation.batch.${batch_index}.${params.n_batches}.txt --window ${params.cis_window} --permute ${params.n_permutations} --grp-best
    """
}

/*
 * STEP 7 - Merge permutation batches from QTLtools
 */
process merge_permutation_batches {
    tag "${qtl_subset}"
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'

    input:
    tuple val(qtl_subset), file(batch_file_names)

    output:
    tuple val(qtl_subset), file("${qtl_subset}_permuted.tsv.gz")

    script:
    """
    cat ${batch_file_names.join(' ')} | csvtk space2tab | sort -k11n -k12n > merged.txt
    cut -f 1,6,7,8,10,11,12,18,20,21,22 merged.txt | csvtk add-header -t -n molecular_trait_object_id,molecular_trait_id,n_traits,n_variants,variant,chromosome,position,pvalue,beta,p_perm,p_beta | gzip > ${qtl_subset}_permuted.tsv.gz
    """
}

process convert_merged_permutation_txt_to_pq {
    tag "${qtl_subset}"
    publishDir "${params.outdir}/sumstats/${qtl_subset}", mode: 'copy'
    container = 'quay.io/kfkf33/duckdb_env:v24.01.1'

    input:
    tuple val(qtl_subset), path(input_file)

    output:
    tuple val(qtl_subset), path("${qtl_subset}.permuted.parquet")

    script:
    """
    $baseDir/bin/convert_txt_to_pq.py \
        -i $input_file \
        -m ${task.memory.toMega() / 1024} \
        -c molecular_trait_object_id,molecular_trait_id,n_traits,n_variants,variant,chromosome,position,pvalue,beta,p_perm,p_beta \
        -o ${qtl_subset}.permuted.parquet
    """
}

/*
 * STEP 8 - Run QTLtools in nominal mode
 */
process run_nominal {
    tag "${qtl_subset} - ${batch_index}/${params.n_batches}"
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'
    
    input:
    each batch_index
    tuple val(qtl_subset), file(bed), file(bed_index), file(fastqtl_bed), file(fastqtl_bed_index), file(vcf), file(vcf_index), file(covariate)

    output:
    tuple val(qtl_subset), file("${qtl_subset}.nominal.batch.${batch_index}.${params.n_batches}.txt"), env(chrom), env(start_pos), env(end_pos)

    script:
    """
	fastQTL --vcf $vcf --bed $fastqtl_bed --cov $covariate \\
        --chunk $batch_index ${params.n_batches} \\
        --out ${qtl_subset}.nominal.batch.${batch_index}.${params.n_batches}.txt \\
        --window ${params.cis_window} \\
        --ma-sample-threshold 1

    retries=5
    while [[ ! -s .command.out && \$retries -gt 0 ]]; do
        echo "Waiting for .command.out to be available..."
        sleep 2
        ((\$retries--))
    done

    if [[ ! -s .command.out ]]; then
        echo "Error: .command.out file not found or empty. Exiting."
        exit 1
    fi

    analyzed_region=\$(grep -A1 "Reading genotype data" .command.out | grep "region =" | head -n 1 | awk -F'=' '{print \$2}' | sed 's/[:\\-]/_/g')

    if [[ -z "\$analyzed_region" ]]; then
        echo "Error: analyzed_region is empty. Exiting."
        exit 1
    fi

    # Extract the genotype data region from the .command.out file
    read chrom start_pos end_pos <<< \$(echo "\$analyzed_region" | awk -F'_' '{print \$1, \$2, \$3}')
    """
}
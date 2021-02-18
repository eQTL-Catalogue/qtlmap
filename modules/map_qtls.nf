/*
 * STEP 6 - Run QTLtools in permutation mode
 */
process run_permutation {
    tag "${qtl_subset} - ${batch_index}/${params.n_batches}"
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'

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
    publishDir "${params.outdir}/sumstats", mode: 'copy'
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'

    input:
    tuple val(qtl_subset), file(batch_file_names)

    output:
    tuple val(qtl_subset), file("${qtl_subset}.permuted.txt.gz")

    script:
    """
    cat ${batch_file_names.join(' ')} | csvtk space2tab | sort -k11n -k12n > merged.txt
    cut -f 1,6,7,8,10,11,12,18,19,20,21 merged.txt | csvtk add-header -t -n molecular_trait_object_id,molecular_trait_id,n_traits,n_variants,variant,chromosome,position,pvalue,beta,p_perm,p_beta | bgzip > ${qtl_subset}.permuted.txt.gz
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
    tuple val(qtl_subset), file("${qtl_subset}.nominal.batch.${batch_index}.${params.n_batches}.txt")

    script:
    """
	fastQTL --vcf $vcf --bed $fastqtl_bed --cov $covariate \\
        --chunk $batch_index ${params.n_batches} \\
        --out ${qtl_subset}.nominal.batch.${batch_index}.${params.n_batches}.txt \\
        --window ${params.cis_window} \\
        --ma-sample-threshold 1
    """
}

/*
 * STEP 9 - Merge nominal batches from QTLtools
 */
process merge_nominal_batches {
    tag "${qtl_subset}"
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'

    input:
    tuple val(qtl_subset), file(batch_file_names)  

    output:
    tuple val(qtl_subset), file("${qtl_subset}.nominal.tab.txt.gz")

    script:
    """
    cat ${batch_file_names.join(' ')} | \\
        csvtk space2tab -T | \\
        csvtk sep -H -t -f 2 -s "_" | \\
        csvtk replace -t -H -f 10 -p ^chr | \\
        csvtk cut -t -f1,10,11,12,13,2,4,5,6,7,8,9 | \\
        bgzip > ${qtl_subset}.nominal.tab.txt.gz
    """
}

/*
 * STEP 11 - Sort fastQTL nominal pass outout and add header
 */
process sort_qtltools_output {
    tag "${qtl_subset}"
    publishDir path: { !params.reformat_sumstats ? "${params.outdir}/sumstats" : params.outdir },
            saveAs: { !params.reformat_sumstats ? it : null }, mode: 'copy'
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'

    input:
    tuple val(qtl_subset), file(nominal_merged)

    output:
    tuple val(qtl_subset), file("${qtl_subset}.nominal.sorted.norsid.tsv.gz")

    script:
    """
    gzip -dc $nominal_merged | LANG=C sort -k2,2 -k3,3n -S11G --parallel=8 | uniq | \\
        bgzip > ${qtl_subset}.nominal.sorted.norsid.tsv.gz
    """
}
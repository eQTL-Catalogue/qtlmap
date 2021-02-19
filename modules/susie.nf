process run_susie{
    container = 'quay.io/eqtlcatalogue/susie-finemapping:v20.08.1'

    input:
    tuple val(qtl_subset), file(expression_matrix), file(phenotype_meta), file(sample_meta), file(phenotype_list), file(covariates), file(genotype_matrix), file(genotype_matrix_index)
    each batch_index

    output:
    tuple val(qtl_subset), file("${qtl_subset}.${batch_index}_${params.n_batches}.txt"), file("${qtl_subset}.${batch_index}_${params.n_batches}.cred.txt"), file("${qtl_subset}.${batch_index}_${params.n_batches}.snp.txt")

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
     --permuted true
    """
}

process merge_susie{
    container = 'quay.io/eqtlcatalogue/susie-finemapping:v20.08.1'

    publishDir {params.save_susie_full ? {"${params.outdir}/susie_full/", mode: 'copy', pattern: "*.cred.txt.gz"} : false }
    publishDir {params.save_susie_full ? {"${params.outdir}/susie_full/", mode: 'copy', pattern: "*.snp.txt.gz"} : false }

    input:
    tuple val(qtl_subset), file(in_cs_variant_batch_names), file(credible_set_batch_names), file(variant_batch_names)
    
    output:
    tuple val(qtl_subset), file("${qtl_subset}.txt.gz"), file("${qtl_subset}.cred.txt.gz"), file("${qtl_subset}.snp.txt.gz")

    script:
    """
    awk 'NR == 1 || FNR > 1{print}' ${in_cs_variant_batch_names.join(' ')} | gzip -c > ${qtl_subset}.txt.gz
    awk 'NR == 1 || FNR > 1{print}' ${credible_set_batch_names.join(' ')} | gzip -c > ${qtl_subset}.cred.txt.gz
    awk 'NR == 1 || FNR > 1{print}' ${variant_batch_names.join(' ')} | gzip -c > ${qtl_subset}.snp.txt.gz
    """
}

process sort_susie{
    container = 'quay.io/eqtlcatalogue/susie-finemapping:v20.08.1'

    publishDir "${params.outdir}/susie/", mode: 'copy', pattern: "*.purity_filtered.txt.gz"

    input:
    tuple val(qtl_subset), file(merged_susie_output), file(susie_cred_output), file(susie_snp_output)

    output:
    tuple val(qtl_subset), file("${qtl_subset}.purity_filtered.txt.gz")

    script:
    """
    gunzip -c ${merged_susie_output} > susie_merged.txt
    (head -n 1 susie_merged.txt && tail -n +2 susie_merged.txt | sort -k3 -k4n ) | gzip > ${qtl_subset}.purity_filtered.txt.gz
    """
}
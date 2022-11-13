process run_susie{
    container = 'quay.io/eqtlcatalogue/susier:v21.10.2'

    input:
    tuple val(qtl_subset), file(expression_matrix), file(phenotype_meta), file(sample_meta), file(phenotype_list), file(covariates), file(genotype_matrix), file(genotype_matrix_index)
    each batch_index

    output:
    tuple val(qtl_subset), file("${qtl_subset}.${batch_index}_${params.n_batches}.txt"), file("${qtl_subset}.${batch_index}_${params.n_batches}.cred.txt"), file("${qtl_subset}.${batch_index}_${params.n_batches}.snp.txt"), file("${qtl_subset}.${batch_index}_${params.n_batches}.lbf_variable.txt")

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
     --permuted true\
     --skip_full ${params.susie_skip_full}
    """
}

process merge_susie{
    container = 'quay.io/eqtlcatalogue/susie-finemapping:v20.08.1'

    publishDir "${params.outdir}/susie_full/", mode: 'copy', pattern: "*.cred.txt.gz"
    publishDir "${params.outdir}/susie_full/", mode: 'copy', pattern: "*.snp.txt.gz"
    publishDir "${params.outdir}/susie_lbf/", mode: 'copy', pattern: "*.lbf_variable.txt.gz"

    input:
    tuple val(qtl_subset), file(in_cs_variant_batch_names), file(credible_set_batch_names), file(variant_batch_names), file(lbf_variable_batch_names)
    
    output:
    tuple val(qtl_subset), file("${qtl_subset}.txt.gz"), file("${qtl_subset}.cred.txt.gz"), file("${qtl_subset}.snp.txt.gz"), file("${qtl_subset}.lbf_variable.txt.gz")

    script:
    """
    awk 'NR == 1 || FNR > 1{print}' ${in_cs_variant_batch_names.join(' ')} | gzip -c > ${qtl_subset}.txt.gz
    awk 'NR == 1 || FNR > 1{print}' ${credible_set_batch_names.join(' ')} | gzip -c > ${qtl_subset}.cred.txt.gz
    awk 'NR == 1 || FNR > 1{print}' ${variant_batch_names.join(' ')} | gzip -c > ${qtl_subset}.snp.txt.gz
    awk 'NR == 1 || FNR > 1{print}' ${lbf_variable_batch_names.join(' ')} | gzip -c > ${qtl_subset}.lbf_variable.txt.gz
    """
}

process sort_susie{
    container = 'quay.io/eqtlcatalogue/susie-finemapping:v20.08.1'

    publishDir "${params.outdir}/susie/", mode: 'copy', pattern: "*.purity_filtered.txt.gz"

    input:
    tuple val(qtl_subset), file(merged_susie_output), file(susie_cred_output), file(susie_snp_output), file(susie_lbf_output)

    output:
    tuple val(qtl_subset), file("${qtl_subset}.purity_filtered.txt.gz")

    script:
    """
    gunzip -c ${merged_susie_output} > susie_merged.txt
    (head -n 1 susie_merged.txt && tail -n +2 susie_merged.txt | sort -k3 -k4n ) | gzip > ${qtl_subset}.purity_filtered.txt.gz
    """
}

process extract_cs_variants{
    container = 'quay.io/eqtlcatalogue/susie-finemapping:v20.08.1'

    input:
    tuple val(qtl_subset), file(credible_sets), file(qtl_ss), file(qtl_ss_index)

    output:
    tuple val(qtl_subset), file(credible_sets), file("${qtl_subset}.extracted_sumstats.tsv.gz")

    script:
    """
    #Extract variant coordinates from the credible set file
    csvtk cut -t -T -f chromosome,position ${credible_sets} | tail -n +2 | sort -k1n -k2n | uniq > selected_regions.tsv

    #Extract variants from the summary stats file
    set +o pipefail; zcat ${qtl_ss} | head -n1 | gzip > header.txt.gz
    set +o pipefail; tabix -R selected_regions.tsv ${qtl_ss} | gzip > filtered_sumstats.tsv.gz
    set +o pipefail; zcat header.txt.gz filtered_sumstats.tsv.gz | gzip > ${qtl_subset}.extracted_sumstats.tsv.gz
    """
}

process merge_cs_sumstats{
    publishDir "${params.outdir}/susie_merged/", mode: 'copy', pattern: "*.credible_sets.tsv.gz"
    container = 'quay.io/eqtlcatalogue/susie-finemapping:v20.08.1'

    input:
    tuple val(qtl_subset), file(credible_sets), file(sumstats)

    output:
    tuple val(qtl_subset), file("${qtl_subset}.credible_sets.tsv.gz")

    script:
    """
    Rscript $baseDir/bin/susie_merge_cs.R --cs_results ${credible_sets}\
     --sumstats ${sumstats}\
     --out ${qtl_subset}.credible_sets.tsv.gz
    """
}
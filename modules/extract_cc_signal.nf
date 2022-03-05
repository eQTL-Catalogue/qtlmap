/*
 * Prepare molecular trait data by:
 * 1. Build connected components where the credible sets sharing at least 1 variant become one connected component
 * 2. Extract nominal summary statistics of lead signal (phenotype id) generared by build_connected_components process
 */
process build_connected_components {
    tag "${qtl_subset}"
    publishDir "${params.outdir}/cc_signal/", mode: 'copy'
    container = 'quay.io/eqtlcatalogue/concomp:latest'

    input:
    tuple val(qtl_subset), file(susie_purity_filtered), file(expression_matrix), file(phenotype_metadata), file(sample_metadata)

    output: 
    tuple val(qtl_subset), file("*.txt")

    script:
    """
    Rscript $baseDir/bin/build_connected_components.R \\
        -s "$susie_purity_filtered" \\
        -p "$phenotype_metadata" \\
        -q "$qtl_subset" \\
        -o "./${qtl_subset}_cc_phenotypes.txt" \\
        -c ${params.cs_size_threshold} 
    """
}

// Extract nominal summary statistics of lead signal (phenotype id) generared by build_connected_components process
process extract_lead_cc_signal {
    tag "${qtl_subset}"
    publishDir "${params.outdir}/cc_signal/", mode: 'copy'
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'

    input:
    tuple val(qtl_subset), file(cc_signals_phenotypes), file(nominal_sumstats), file(nominal_sumstats_index)

    output:
    tuple val(qtl_subset), file("${qtl_subset}.cc.tsv.gz"), file("${qtl_subset}.cc.tsv.gz.tbi")

    script:
    """
    zcat $nominal_sumstats | csvtk -t grep -P $cc_signals_phenotypes | bgzip > ${qtl_subset}.nominal.sorted.bgzip.tsv.gz
    mv ${qtl_subset}.nominal.sorted.bgzip.tsv.gz ${qtl_subset}.cc.tsv.gz
    tabix -s2 -b3 -e3 -S1 -f ${qtl_subset}.cc.tsv.gz
    """
}


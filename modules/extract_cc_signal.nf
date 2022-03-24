

// Extract nominal summary statistics of lead signal (phenotype id) generared by build_connected_components process
process extract_lead_cc_signal {
    tag "${qtl_subset}"
    publishDir "${params.outdir}/sumstats", mode: 'copy'
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'

    input:
    tuple val(qtl_subset), file(credible_sets), file(qtl_ss), file(qtl_ss_index)

    output:
    tuple val(qtl_subset), file("${qtl_subset}.cc.tsv.gz"), file("${qtl_subset}.cc.tsv.gz.tbi")

    script:
    """
    csvtk -t cut -f molecular_trait_id $credible_sets | csvtk -t uniq > cc_signals_phenotypes.txt
    zcat $qtl_ss | csvtk -t grep -P cc_signals_phenotypes.txt | bgzip > ${qtl_subset}.cc.tsv.gz
    tabix -s2 -b3 -e3 -S1 -f ${qtl_subset}.cc.tsv.gz
    """
}


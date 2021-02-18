process vcf_set_variant_ids {
    tag "${qtl_subset}"
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'

    input:
    tuple val(qtl_subset), file(vcf)

    output:
    tuple val(qtl_subset), file("${vcf.simpleName}_renamed.vcf.gz")

    script:
    """
    bcftools annotate --set-id 'chr%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' $vcf -Oz -o ${vcf.simpleName}_renamed.vcf.gz
    """
}
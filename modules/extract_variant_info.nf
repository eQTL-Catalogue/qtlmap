// Extract variant information from VCF
process extract_variant_info {
    tag "${qtl_subset}"
    publishDir "${params.outdir}/varinfo", mode: 'copy'
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'

    input:
    tuple val(qtl_subset), file(vcf)
    
    output:
    tuple val(qtl_subset), file("${qtl_subset}.variant_information.txt.gz")

    script:
    if (params.vcf_has_R2_field) {
        """
        set +o pipefail; bcftools +fill-tags $vcf | bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%TYPE\\t%AC\\t%AN\\t%MAF\\t%R2\\n' | gzip > ${qtl_subset}.variant_information.txt.gz
        """
    } else {
        """
        set +o pipefail; bcftools +fill-tags $vcf | bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%TYPE\\t%AC\\t%AN\\t%MAF\\tNA\\n' | gzip > ${qtl_subset}.variant_information.txt.gz
        """
    }
}
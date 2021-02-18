
process extract_samples_from_vcf {
    tag "${qtl_subset}"
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'

    input:
    tuple val(qtl_subset), file(genotype_vcf), file(sample_names)

    output:
    tuple val(qtl_subset), file("${sample_names.simpleName}.vcf.gz"), emit: vcf 
    tuple val(qtl_subset), file("${sample_names.simpleName}.vcf.gz.tbi"), emit: index

    script:
    """
    bcftools view -S $sample_names $genotype_vcf -Oz -o ${sample_names.simpleName}.vcf.gz
    tabix -p vcf ${sample_names.simpleName}.vcf.gz
    """
}
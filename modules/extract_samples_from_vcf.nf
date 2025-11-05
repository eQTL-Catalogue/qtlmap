process extract_samples_from_vcf {
    tag "${qtl_subset}"
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'

    input:
    tuple val(qtl_subset), file(genotype_vcf), file(sample_names)

    output:
    tuple val(qtl_subset), file("${sample_names.simpleName}.vcf.gz"), emit: vcf 
    tuple val(qtl_subset), file("${sample_names.simpleName}.vcf.gz.tbi"), emit: index

    script:
    if (params.vcf_extract_samples){
        """
        bcftools view --threads ${task.cpus} -S $sample_names $genotype_vcf -Oz -o ${sample_names.simpleName}_extract.vcf.gz
        bcftools +fill-tags --threads ${task.cpus} ${sample_names.simpleName}_extract.vcf.gz -Oz -o ${sample_names.simpleName}_extract_filltags.vcf.gz
        bcftools view --threads ${task.cpus} -i 'AN[0]*MAF[0]>5 & MAF[0]>0.01' ${sample_names.simpleName}_extract_filltags.vcf.gz -Oz -o ${sample_names.simpleName}.vcf.gz
        tabix -p vcf ${sample_names.simpleName}.vcf.gz
        """
    } else {
        """
        cp $genotype_vcf ${sample_names.simpleName}.vcf.gz
        tabix -p vcf ${sample_names.simpleName}.vcf.gz
        """
    }
}
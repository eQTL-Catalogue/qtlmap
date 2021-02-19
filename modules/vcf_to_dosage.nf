process vcf_to_dosage{
    container = 'quay.io/eqtlcatalogue/susie-finemapping:v20.08.1'

    input:
    tuple val(qtl_subset), file(vcf)

    output:
    tuple val(qtl_subset), file("${vcf.simpleName}.dose.tsv.gz"), file("${vcf.simpleName}.dose.tsv.gz.tbi")

    script:
    if(params.vcf_genotype_field == 'DS'){
        """
        #Extract header
        printf 'CHROM\\nPOS\\nREF\\nALT\\n' > 4_columns.tsv
        bcftools query -l ${vcf} > sample_list.tsv
        cat 4_columns.tsv sample_list.tsv > header.tsv
        csvtk transpose header.tsv -T | gzip > header_row.tsv.gz

        #Extract dosage and merge
        bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%DS]\\n" ${vcf} | gzip > dose_matrix.tsv.gz
        zcat header_row.tsv.gz dose_matrix.tsv.gz | bgzip > ${vcf.simpleName}.dose.tsv.gz
        tabix -s1 -b2 -e2 -S1 ${vcf.simpleName}.dose.tsv.gz
        """
    } else if (params.vcf_genotype_field == 'GT'){
        """
        #Extract header
        printf 'CHROM\\nPOS\\nREF\\nALT\\n' > 4_columns.tsv
        bcftools query -l ${vcf} > sample_list.tsv
        cat 4_columns.tsv sample_list.tsv > header.tsv
        csvtk transpose header.tsv -T | gzip > header_row.tsv.gz

        #Extract dosage and merge
        bcftools +dosage ${vcf} -- -t GT | tail -n+2 | gzip > dose_matrix.tsv.gz
        zcat header_row.tsv.gz dose_matrix.tsv.gz | bgzip > ${vcf.simpleName}.dose.tsv.gz
        tabix -s1 -b2 -e2 -S1 ${vcf.simpleName}.dose.tsv.gz
        """
    } 
}
